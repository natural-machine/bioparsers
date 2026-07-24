#!/usr/bin/env python
"""Build Recipe: Capture the Pfam annotation field set + caption-ready text, by Pfam.

A sibling of ``build_pfam_legacy_by_pfam.py``. Where that recipe assembles a
legacy ``[final]text_caption``, this one keeps the structured fields and the
cleaned, caption-ready text **side by side** and assembles **no** caption — so a
downstream trainer can compose captions on the fly (annotation dropout, field
randomization, reordering). It is the Pfam counterpart of
``build_swissprot_caption_fields_by_pfam.py``.

Like the legacy Pfam recipe, an entry is a **join of three sources** performed
by ``bioparsers.builders.pfam.run_pfam_join`` (see that recipe's docstring for
the member / families / uniprot inputs and the ``--uniprot-cache`` flag).

Each output record carries:

  - ``fields``: the raw extracted fields — ``family_name`` / ``family_description``
    (from the Pfam family), then the UniProt-derived fields when the member
    resolves to a UniProt entry (per-CC-topic lists, ``protein_name``,
    ``gene_ontology`` aspect-ordered, ``lineage``, ...). Only populated keys.
  - ``caption_fields``: a parallel dict where each present field is the
    *cleaned, concatenated* string it would contribute to a caption — no
    ``LABEL:`` prefix, no ``"The organism lineage is"`` preamble, no caption
    assembly, and (unlike the legacy caption) no forced-empty fields.

Field extraction matches the legacy recipe (CC topics as lists of all blocks,
GENE ONTOLOGY = GO term names grouped by aspect C -> F -> P, etc.).

Usage:
    python recipes/build_pfam_caption_fields_by_pfam.py data/pfam_fasta.jsonl.gz \\
        --pfam-ids PF00018 \\
        --pfam-families data/pfam.jsonl.gz \\
        --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \\
        --uniprot-cache data/.uniprot_cache/ \\
        -o outputs/pfam_caption_fields_SH3.jsonl
"""

import argparse
import re
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.pfam import filters, helpers, run_pfam_join
from bioparsers.builders.uniprot import helpers as uniprot_helpers


# (caption-field label, output field key) for the UniProt-derived section, in
# the legacy Pfam-caption order (kept for a stable, documented field order; see
# build_pfam_legacy_by_pfam.py for how it was derived from FINAL_SH3_pfam.csv).
SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("FUNCTION", "function"),
    ("CATALYTIC ACTIVITY", "catalytic_activity"),
    ("ACTIVITY REGULATION", "activity_regulation"),
    ("MISCELLANEOUS", "miscellaneous"),
    ("DOMAIN", "domain"),
    ("PATHWAY", "pathway"),
    ("SUBCELLULAR LOCATION", "subcellular_location"),
    ("SUBUNIT", "subunit"),
    ("COFACTOR", "cofactor"),
    ("SIMILARITY", "similarity"),
    ("BIOPHYSICOCHEMICAL PROPERTIES", "biophysicochemical_properties"),
    ("TISSUE SPECIFICITY", "tissue_specificity"),
    ("PTM", "ptm"),
    ("INDUCTION", "induction"),
    ("BIOTECHNOLOGY", "biotechnology"),
    ("GENE ONTOLOGY", "gene_ontology"),
    ("LINEAGE", "lineage"),
]

# Order of keys in the emitted dicts: family fields first, then the UniProt set.
FIELD_ORDER = ["family_name", "family_description"] + [key for _, key in SPEC]

_SIMPLE_TOPICS = {
    "FUNCTION": "function",
    "ACTIVITY REGULATION": "activity_regulation",
    "MISCELLANEOUS": "miscellaneous",
    "DOMAIN": "domain",
    "PATHWAY": "pathway",
    "SUBUNIT": "subunit",
    "SIMILARITY": "similarity",
    "BIOPHYSICOCHEMICAL PROPERTIES": "biophysicochemical_properties",
    "TISSUE SPECIFICITY": "tissue_specificity",
    "PTM": "ptm",
    "INDUCTION": "induction",
    "BIOTECHNOLOGY": "biotechnology",
}

# Cleaning regexes, ported verbatim from biom3.dbio.caption for parity.
_PUBMED_RE = re.compile(r"\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)")
_EVIDENCE_RE = re.compile(r"\s*\{ECO:\d+.*?\}")
_MULTI_DOT_RE = re.compile(r"\.(\s*\.)+")
_MULTI_SPACE_RE = re.compile(r"  +")
_ISOFORM_RE = re.compile(r"\[Isoform[^\]]*\]:\s*")
_ASPECT_ORDER = {"C": 0, "F": 1, "P": 2}


# ===========================================================================
# Builder + CLI
# ===========================================================================

class PfamCaptionFieldsBuilder(Builder):
    """Pfam annotation fields + caption-ready text, by Pfam (no assembled caption).

    Consumes the joined member records produced by
    ``bioparsers.builders.pfam.run_pfam_join`` (each a Pfam member dict with
    ``family`` metadata and the matched ``uniprot`` record, or ``None``).

    Output record::

        {
          "accession": str | None,   # member UniProt accession ("" -> None)
          "sequence":  str,          # domain region sequence
          "region":    str | None,   # aligned region, e.g. "55-110"
          "pfam_ids":  [str],        # [family Pfam accession]
          "fields": {                # raw extracted fields (the inputs)
            "family_name": str, "family_description": str,
            "protein_name": str, "function": [str], ...,   # when a UniProt match
            "gene_ontology": [str], "lineage": [str],
          },
          "caption_fields": {        # cleaned, concatenated text per field
            "family_name": str, "family_description": str,
            "protein_name": str, "function": str, ..., "lineage": str,
          }
        }

    No ``[final]text_caption`` is assembled — ``caption_fields`` holds each
    field's caption-ready string (evidence/PubMed stripped, blocks joined, no
    ``LABEL:`` prefix) so callers can compose captions as they wish.

    Options:
      - ``min_length``: drop members whose domain sequence is shorter than this.
    """

    name = "pfam_caption_fields"

    def __init__(self, *, min_length: int = 0):
        self.min_length = min_length

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for member in records:
            if not keep_len(member):
                continue

            family = member.get("family") or {}
            uniprot = member.get("uniprot")
            fields = {}
            name = helpers.family_name(family)
            if name:
                fields["family_name"] = name
            desc = helpers.family_description(family)
            if desc:
                fields["family_description"] = desc
            if uniprot is not None:
                fields.update(extract_fields(uniprot))

            yield {
                "accession": member.get("accession") or None,
                "sequence": member.get("sequence"),
                "region": member.get("region"),
                "pfam_ids": [member["pfam_accession"]] if member.get("pfam_accession") else [],
                "fields": fields,
                "caption_fields": caption_fields(fields),
            }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=PfamCaptionFieldsBuilder.description)
    p.add_argument("input", help="parsed Pfam member FASTA JSONL (plain or .gz)")
    p.add_argument("--pfam-ids", nargs="+", required=True, metavar="PFAM_ID",
                   dest="pfam_ids", help="one or more Pfam accessions, e.g. PF00018")
    p.add_argument("-o", "--output", required=True,
                   help="output JSONL path; per-ID mode inserts the Pfam ID")
    p.add_argument("--join", action="store_true",
                   help="write a single union file instead of one per Pfam ID")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in each output's build manifest")
    p.add_argument("--pfam-families", required=True, metavar="JSONL",
                   help="parsed Pfam family JSONL (bioparsers pfam) supplying "
                        "family_name / family_description (required)")
    p.add_argument("--uniprot", nargs="+", metavar="JSONL",
                   default=["data/uniprot_sprot.jsonl.gz",
                            "data/uniprot_trembl.jsonl.gz"],
                   help="parsed UniProt JSONL file(s) joined on the member "
                        "accession (default: Swiss-Prot then TrEMBL)")
    p.add_argument("--uniprot-cache", metavar="PATH", default=None,
                   help="cache the resolved UniProt subset for reuse; a plain "
                        "path caches the union of all --pfam-ids, a directory "
                        "(trailing '/') keeps one cache file per Pfam ID")
    p.add_argument("--min-length", type=int, default=0,
                   help="drop members whose domain sequence is shorter than this")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    families = helpers.load_family_metadata(args.pfam_families, args.pfam_ids)
    builder = PfamCaptionFieldsBuilder(min_length=args.min_length)
    run_pfam_join(builder, args.input, args.pfam_ids,
                  families=families, uniprot_paths=args.uniprot,
                  output=args.output, join=args.join, gzip=args.gzip,
                  description=args.description, cache_path=args.uniprot_cache)


# ===========================================================================
# Field extraction & caption-text helpers
# ===========================================================================

def _strip_evidence(text: str) -> str:
    return _EVIDENCE_RE.sub("", text).strip()


def _strip_pubmed(text: str) -> str:
    text = _PUBMED_RE.sub("", text)
    text = _MULTI_DOT_RE.sub(".", text)
    text = _MULTI_SPACE_RE.sub(" ", text).strip()
    return text


def _clean_block(text: str) -> str:
    """A stored comment block: evidence- and PubMed-stripped, with the
    doubled-period / double-space artifacts those removals leave repaired."""
    return _strip_pubmed(_strip_evidence(text))


def _all_comments(rec, topic):
    return [c["text"] for c in rec.get("comments", []) if c.get("topic") == topic]


def _full_name(rec):
    desc = rec.get("description") or {}
    for key in ("rec_name", "sub_name"):
        block = desc.get(key)
        if block and block.get("full"):
            return _strip_evidence(block["full"]).rstrip(";").rstrip(".").strip()
    return None


def _go_terms(rec):
    """GO cross-ref term names, aspect-prefix stripped, ordered C -> F -> P
    (within an aspect, source order is preserved)."""
    items = []
    for line in rec.get("cross_references", {}).get("GO", []):
        parts = line.split(";")
        if len(parts) < 3:
            continue
        desc = parts[2].strip().rstrip(".")
        aspect = None
        if len(desc) > 2 and desc[1] == ":":
            aspect, desc = desc[0], desc[2:]
        items.append((aspect, desc))
    items.sort(key=lambda t: _ASPECT_ORDER.get(t[0], 3))
    return [d for _, d in items]


def extract_fields(rec) -> dict:
    """Map a matched UniProt record to the Pfam annotation field set.

    CC-comment fields are collected as **lists of all blocks** of that topic,
    evidence-stripped. ``protein_name`` is a single string; ``gene_ontology``
    (aspect-ordered) and ``lineage`` are lists.
    """
    fields = {}

    name = _full_name(rec)
    if name:
        fields["protein_name"] = name

    for topic, key in _SIMPLE_TOPICS.items():
        blocks = [b for b in (_clean_block(t) for t in _all_comments(rec, topic)) if b]
        if blocks:
            fields[key] = blocks

    reactions = []
    for text in _all_comments(rec, "CATALYTIC ACTIVITY"):
        m = re.search(r"Reaction=([^;]+)", text)
        if m:
            reaction = _clean_block(m.group(1))
            if reaction:
                reactions.append(reaction)
    if reactions:
        fields["catalytic_activity"] = reactions

    cofactors = uniprot_helpers.cofactor_names(rec)
    if cofactors:
        fields["cofactor"] = cofactors

    sublocs = []
    for text in _all_comments(rec, "SUBCELLULAR LOCATION"):
        body = _ISOFORM_RE.sub("", _strip_evidence(text).split("Note=")[0])
        val = _clean_block(body).rstrip(".")
        if val:
            sublocs.append(val)
    if sublocs:
        fields["subcellular_location"] = sublocs

    go = _go_terms(rec)
    if go:
        fields["gene_ontology"] = go

    lineage = rec.get("lineage")
    if lineage:
        fields["lineage"] = list(lineage)

    return fields


def _field_text(key, value) -> str:
    """Render a field value to its bare, concatenated caption text — no
    ``LABEL:`` prefix and (unlike the legacy caption) no ``lineage`` preamble.
    """
    if key in ("lineage", "gene_ontology", "cofactor"):
        return ", ".join(value)
    if isinstance(value, list):
        return " ".join(value)        # join all blocks of a CC topic
    return value                      # str: family_name/description, protein_name


def caption_fields(fields: dict) -> dict:
    """Build the ``{field: cleaned-concatenated-text}`` dict from the extracted
    fields, in :data:`FIELD_ORDER`. Each value is evidence/PubMed-stripped, with
    blocks of a CC topic concatenated; empty fields are omitted.
    """
    out = {}
    for key in FIELD_ORDER:
        value = fields.get(key)
        if not value:
            continue
        text = _strip_pubmed(_strip_evidence(_field_text(key, value))).rstrip(".").strip()
        if text:
            out[key] = text
    return out


if __name__ == "__main__":
    main()
