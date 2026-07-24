#!/usr/bin/env python
"""Build Recipe: Capture the full Swiss-Prot annotation field set, by Pfam.

A sibling of ``build_swissprot_legacy_by_pfam.py``. Where that recipe trims to
the legacy caption field set and emits an assembled ``[final]text_caption``,
this one keeps the **full** set of fields we extract (including GENE ONTOLOGY)
and does **not** assemble a caption. Instead each record carries:

  - ``fields``: the raw extracted fields (per-CC-topic lists, etc.), exactly
    as captured — the inputs, not a rendered string.
  - ``caption_fields``: a parallel dict where each field is the *cleaned,
    concatenated* string it would contribute to a caption (e.g. all DOMAIN
    blocks joined; lineage / GO as a comma-joined list), with no ``LABEL:``
    prefix and no caption assembly.

This keeps the structured data and the caption-ready text side by side, so
downstream callers can compose captions however they like.

Field extraction (CC-comment fields are lists of all blocks, evidence-stripped):
  - PROTEIN NAME = DE RecName.full (SubName.full fallback) — single string.
  - Each CC topic (FUNCTION, DOMAIN, ...) = list of all its blocks.
  - CATALYTIC ACTIVITY = list of each block's ``Reaction=`` prose.
  - SUBCELLULAR LOCATION = list of each block's text before ``Note=`` (no trailing ".").
  - GENE ONTOLOGY = list of GO cross-ref term text (aspect prefix stripped).
  - LINEAGE = list of OC taxa.

The ``family_names`` caption field requires external Pfam metadata (PF
accession -> full family name, e.g. PF00018 -> "SH3 domain"), which a parsed
UniProt record does not carry. ``--pfam-names`` (a two-column ``PF<TAB>name``
TSV, as written by ``scripts/parse_pfam_names.sh``) is required.

Usage (filter Swiss-Prot JSONL to the SH3 Pfam, one file):
    python recipes/build_swissprot_caption_fields_by_pfam.py \\
        data/uniprot_sprot.jsonl.gz \\
        --pfam-ids PF00018 \\
        --pfam-names data/pfam_names.tsv \\
        -o outputs/swissprot_caption_fields_SH3.jsonl
"""

import argparse
import re
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.uniprot import filters, helpers, run_by_pfam


# (caption-field label, output field key) — the full extracted set, in order.
SPEC = [
    ("PROTEIN NAME", "protein_name"),
    ("FUNCTION", "function"),
    ("CATALYTIC ACTIVITY", "catalytic_activity"),
    ("COFACTOR", "cofactor"),
    ("ACTIVITY REGULATION", "activity_regulation"),
    ("BIOPHYSICOCHEMICAL PROPERTIES", "biophysicochemical_properties"),
    ("PATHWAY", "pathway"),
    ("SUBUNIT", "subunit"),
    ("SUBCELLULAR LOCATION", "subcellular_location"),
    ("TISSUE SPECIFICITY", "tissue_specificity"),
    ("DOMAIN", "domain"),
    ("PTM", "ptm"),
    ("SIMILARITY", "similarity"),
    ("MISCELLANEOUS", "miscellaneous"),
    ("INDUCTION", "induction"),
    ("DEVELOPMENTAL STAGE", "developmental_stage"),
    ("BIOTECHNOLOGY", "biotechnology"),
    ("GENE ONTOLOGY", "gene_ontology"),
    ("LINEAGE", "lineage"),
]

# CC topics collected verbatim — every block of the topic, evidence-stripped,
# stored as a list. (CATALYTIC ACTIVITY, SUBCELLULAR LOCATION and COFACTOR are
# structured rather than prose, and get special per-block handling below.)
_SIMPLE_TOPICS = {
    "FUNCTION": "function",
    "ACTIVITY REGULATION": "activity_regulation",
    "BIOPHYSICOCHEMICAL PROPERTIES": "biophysicochemical_properties",
    "PATHWAY": "pathway",
    "SUBUNIT": "subunit",
    "TISSUE SPECIFICITY": "tissue_specificity",
    "DOMAIN": "domain",
    "PTM": "ptm",
    "SIMILARITY": "similarity",
    "MISCELLANEOUS": "miscellaneous",
    "INDUCTION": "induction",
    "DEVELOPMENTAL STAGE": "developmental_stage",
    "BIOTECHNOLOGY": "biotechnology",
}

# Cleaning regexes, ported verbatim from biom3.dbio.caption for parity.
_PUBMED_RE = re.compile(r"\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)")
_EVIDENCE_RE = re.compile(r"\s*\{ECO:\d+.*?\}")
_MULTI_DOT_RE = re.compile(r"\.(\s*\.)+")
_MULTI_SPACE_RE = re.compile(r"  +")
# "[Isoform N]:" / "[Isoform Name]:" markers prefixing SUBCELLULAR LOCATION text.
_ISOFORM_RE = re.compile(r"\[Isoform[^\]]*\]:\s*")


# ===========================================================================
# Builder + CLI
# ===========================================================================

class SwissProtCaptionFieldsBuilder(Builder):
    """Full Swiss-Prot annotation fields + caption-ready text, by Pfam.

    Output record::

        {
          "accession": str,        # primary accession
          "sequence":  str,        # amino-acid sequence
          "pfam_ids":  [str],      # Pfam accessions from DR cross-refs
          "fields": {              # raw extracted fields (the inputs)
            "protein_name": str,
            "function": [str], "domain": [str], ...,   # ALL blocks per CC topic
            "gene_ontology": [str], "lineage": [str],
          },
          "caption_fields": {      # cleaned, concatenated text per field
            "protein_name": str, "function": str, "domain": str, ...,
            "family_names": str,   # from --pfam-names (required)
          }
        }

    No ``[final]text_caption`` is assembled — ``caption_fields`` holds each
    field's caption-ready string (evidence/PubMed stripped, blocks joined,
    no ``LABEL:`` prefix) so callers can compose captions as they wish.

    Options:
      - ``reviewed_only``: keep only Swiss-Prot (Reviewed) entries (default True).
      - ``min_length``: drop entries shorter than this many residues.
      - ``pfam_family_names``: {PF: name} map enabling the family_names entry.
    """

    name = "swissprot_caption_fields"

    def __init__(self, *, reviewed_only: bool = True, min_length: int = 0,
                 pfam_family_names: dict | None = None):
        self.reviewed_only = reviewed_only
        self.min_length = min_length
        self.pfam_family_names = pfam_family_names

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for rec in records:
            if self.reviewed_only and not filters.is_reviewed(rec):
                continue
            if not keep_len(rec):
                continue

            pfam_ids = helpers.pfam_ids(rec)
            family_names = None
            if self.pfam_family_names is not None:
                family_names = [self.pfam_family_names.get(p, p) for p in pfam_ids]

            fields = extract_fields(rec)
            yield {
                "accession": rec.get("primary_accession"),
                "sequence": rec.get("sequence"),
                "pfam_ids": pfam_ids,
                "fields": fields,
                "caption_fields": caption_fields(fields, family_names=family_names),
            }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SwissProtCaptionFieldsBuilder.description)
    p.add_argument("input", help="parsed Swiss-Prot JSONL (plain or .gz)")
    p.add_argument("--pfam-ids", nargs="+", required=True, metavar="PFAM_ID",
                   dest="pfam_ids", help="one or more Pfam accessions, e.g. PF00018")
    p.add_argument("-o", "--output", required=True,
                   help="output JSONL path; per-ID mode inserts the Pfam ID")
    p.add_argument("--join", action="store_true",
                   help="write a single union file instead of one per Pfam ID")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in each output's build manifest")
    p.add_argument("--pfam-names", required=True, metavar="TSV",
                   help="two-column TSV (PF accession<TAB>family name), as written "
                        "by scripts/parse_pfam_names.sh, to fill the family_names "
                        "caption field (required)")
    p.add_argument("--min-length", type=int, default=0)
    return p.parse_args(argv)


def _load_pfam_names(path: str) -> dict:
    pfam_names = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            accession, _, name = line.partition("\t")
            pfam_names[accession] = name
    return pfam_names


def main(argv=None):
    args = parse_args(argv)
    pfam_names = _load_pfam_names(args.pfam_names)
    builder = SwissProtCaptionFieldsBuilder(min_length=args.min_length,
                                            pfam_family_names=pfam_names)
    run_by_pfam(builder, args.input, args.pfam_ids, args.output,
                join=args.join, gzip=args.gzip, description=args.description)


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
    out = []
    for line in rec.get("cross_references", {}).get("GO", []):
        parts = line.split(";")
        if len(parts) >= 3:
            desc = parts[2].strip().rstrip(".")
            if len(desc) > 2 and desc[1] == ":":  # drop C:/F:/P: aspect prefix
                desc = desc[2:]
            out.append(desc)
    return out


def extract_fields(rec) -> dict:
    """Map a parsed UniProt record to the full annotation field set.

    CC-comment fields are collected as **lists of all blocks** of that topic
    (not just the first), evidence-stripped. ``protein_name`` is a single
    string; ``gene_ontology`` and ``lineage`` are lists.
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

    cofactors = helpers.cofactor_names(rec)
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
    return value


def caption_fields(fields: dict, family_names=None) -> dict:
    """Build the ``{field: cleaned-concatenated-text}`` dict from the extracted
    fields, in SPEC order. Each value is evidence/PubMed-stripped, with blocks
    of a CC topic concatenated. ``family_names`` (from external Pfam metadata)
    is added as ``family_names`` when supplied.
    """
    out = {}
    for _, key in SPEC:
        value = fields.get(key)
        if not value:
            continue
        text = _strip_pubmed(_strip_evidence(_field_text(key, value))).rstrip(".").strip()
        if text:
            out[key] = text
    if family_names:
        out["family_names"] = ", ".join(family_names)
    return out


if __name__ == "__main__":
    main()
