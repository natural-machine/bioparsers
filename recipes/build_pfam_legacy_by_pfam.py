#!/usr/bin/env python
"""Build Recipe: Capture fields used in legacy BioM3 Pfam captions, by Pfam.

The Pfam-section counterpart to ``build_swissprot_legacy_by_pfam.py``. It
reproduces (approximately — see below) the **Pfam** section of the legacy BioM3
SH3 finetuning dataset ``FINAL_SH3_pfam.csv`` (columns: primary_Accession,
protein_sequence, [final]text_caption, pfam_label), where each row pairs a
**Pfam domain region** sequence with a caption.

Unlike the Swiss-Prot recipe (a single streaming pass over UniProt), a Pfam
entry is a **join of three sources**, performed by
``bioparsers.builders.pfam.run_pfam_join``:

  - members  : a parsed Pfam member FASTA JSONL (``bioparsers pfam-fasta``) —
               one redundancy-reduced domain sequence per family member, the
               ``protein_sequence`` of each row, carrying the member accession.
  - families : ``{PF: {name, description}}`` family metadata (``--pfam-families``,
               a parsed Pfam family JSONL) — supplies FAMILY NAME / FAMILY
               DESCRIPTION, the same for every member of a family.
  - uniprot  : one or more parsed UniProt JSONL files (``--uniprot``; Swiss-Prot
               + TrEMBL by default) — joined on the member accession to supply
               PROTEIN NAME, FUNCTION, ..., GENE ONTOLOGY, LINEAGE. Most Pfam
               members are TrEMBL, so TrEMBL is needed for full fidelity.

Caption form (derived empirically from ``FINAL_SH3_pfam.csv``)::

    "FAMILY NAME: <name>. FAMILY DESCRIPTION: <blurb>.PROTEIN NAME: <...>.
     FUNCTION: <...>. ... GENE ONTOLOGY: <a, b, c>. LINEAGE: The organism
     lineage is <x, y, z>."

  - The FAMILY section (NAME + DESCRIPTION) is always present (family metadata).
  - The UniProt section is appended **directly** (no separator) when the member
    resolves to a UniProt entry; its fields are joined by ". " in the order
    below. A member with no accession / no UniProt match gets the FAMILY
    section only.
  - PROTEIN NAME, GENE ONTOLOGY, and LINEAGE are emitted whenever a UniProt
    entry is found, even when empty (matching the legacy "GENE ONTOLOGY: ." for
    entries with no GO terms); the other fields appear only when populated.
  - GENE ONTOLOGY is the GO cross-ref term names, grouped by aspect C -> F -> P
    (the legacy ordering); LINEAGE is the OC taxa list.

Exact row-for-row reproduction is not possible: the published dataset was built
against an older Pfam release, so the current redundancy-reduced member set
differs (release drift). Like the Swiss-Prot recipe, this recreates the Pfam
section's *form and content* approximately.

Usage (reproduce the SH3 Pfam section):
    python recipes/build_pfam_legacy_by_pfam.py \\
        data/pfam_fasta.jsonl.gz \\
        --pfam-ids PF00018 \\
        --pfam-families data/pfam.jsonl.gz \\
        --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \\
        -o outputs/pfam_legacy_SH3.jsonl
"""

import argparse
import re
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.pfam import filters, helpers, run_pfam_join
from bioparsers.builders.uniprot import helpers as uniprot_helpers


# (caption label, output field key) for the UniProt-derived section, in the
# legacy Pfam-caption order. This order was derived empirically from
# FINAL_SH3_pfam.csv (a clean total order over the fields seen in the SH3 rows);
# fields not observed in the SH3 data (TISSUE SPECIFICITY, PTM, ...) are kept in
# the spec for generality, placed before GENE ONTOLOGY in their UniProt order.
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
    # Present in the broader Pfam field set but absent from the SH3 rows; kept
    # for generality (no SH3-row effect since they are never populated there).
    ("BIOPHYSICOCHEMICAL PROPERTIES", "biophysicochemical_properties"),
    ("TISSUE SPECIFICITY", "tissue_specificity"),
    ("PTM", "ptm"),
    ("INDUCTION", "induction"),
    ("BIOTECHNOLOGY", "biotechnology"),
    ("GENE ONTOLOGY", "gene_ontology"),
    ("LINEAGE", "lineage"),
]

# Fields emitted whenever a UniProt entry is found, even when their value is
# empty (the legacy always writes these three once a member resolves).
_FORCED = {"protein_name", "gene_ontology", "lineage"}

# CC topics collected verbatim — every block of the topic, evidence-stripped,
# stored as a list. (CATALYTIC ACTIVITY, SUBCELLULAR LOCATION and COFACTOR are
# structured rather than prose, and get special per-block handling below.)
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

# Caption cleaning regexes, ported verbatim from biom3.dbio.caption for parity.
_PUBMED_RE = re.compile(r"\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)")
_EVIDENCE_RE = re.compile(r"\s*\{ECO:\d+.*?\}")
_MULTI_DOT_RE = re.compile(r"\.(\s*\.)+")
_MULTI_SPACE_RE = re.compile(r"  +")
# "[Isoform N]:" / "[Isoform Name]:" markers prefixing SUBCELLULAR LOCATION text.
_ISOFORM_RE = re.compile(r"\[Isoform[^\]]*\]:\s*")

# GO aspect ordering used by the legacy caption: Component, then Function, then
# Process. Within an aspect, source order is preserved (stable sort).
_ASPECT_ORDER = {"C": 0, "F": 1, "P": 2}


# ===========================================================================
# Builder + CLI
# ===========================================================================

class PfamLegacyBuilder(Builder):
    """Legacy BioM3 Pfam caption records.

    Consumes the joined member records produced by
    ``bioparsers.builders.pfam.run_pfam_join`` — each a Pfam member dict
    (``accession``, ``region``, domain ``sequence``, ``pfam_accession``)
    augmented with ``family`` (``{name, description}``) and ``uniprot`` (the
    matched UniProt record dict, or ``None``).

    Output record::

        {
          "accession": str | None,   # member UniProt accession ("" -> None)
          "sequence":  str,          # domain region sequence
          "region":    str | None,   # aligned region, e.g. "55-110"
          "pfam_ids":  [str],        # [family Pfam accession]
          "caption":   str,          # legacy-style [final]text_caption
          "fields": {                # the structured fields the caption is built from
            "family_name": str, "family_description": str,
            "protein_name": str,             # when a UniProt entry is matched
            "function": [str], "domain": [str], ...,   # ALL blocks per CC topic
            "gene_ontology": [str],          # aspect-ordered GO term names
            "lineage": [str],                # OC taxa list
          }
        }

    ``caption`` reproduces the legacy Pfam caption: a FAMILY NAME / FAMILY
    DESCRIPTION section (always), with the UniProt-derived section appended
    directly when the member resolves to a UniProt entry. See the module
    docstring for the exact form.

    Options:
      - ``min_length``: drop members whose domain sequence is shorter than this.
      - ``include_fields``: also emit the structured ``fields`` dict (default).
    """

    name = "pfam_legacy"

    def __init__(self, *, min_length: int = 0, include_fields: bool = True):
        self.min_length = min_length
        self.include_fields = include_fields

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

            out = {
                "accession": member.get("accession") or None,
                "sequence": member.get("sequence"),
                "region": member.get("region"),
                "pfam_ids": [member["pfam_accession"]] if member.get("pfam_accession") else [],
                "caption": compose_caption(family, fields,
                                           has_uniprot=uniprot is not None),
            }
            if self.include_fields:
                out["fields"] = fields
            yield out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=PfamLegacyBuilder.description)
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
                        "FAMILY NAME / FAMILY DESCRIPTION (required)")
    p.add_argument("--uniprot", nargs="+", metavar="JSONL",
                   default=["data/uniprot_sprot.jsonl.gz",
                            "data/uniprot_trembl.jsonl.gz"],
                   help="parsed UniProt JSONL file(s) joined on the member "
                        "accession (default: Swiss-Prot then TrEMBL)")
    p.add_argument("--uniprot-cache", metavar="PATH", default=None,
                   help="cache the resolved UniProt subset for reuse; a later "
                        "run whose accessions are a subset (same --uniprot "
                        "sources) skips the UniProt scan. A plain path caches "
                        "the union of all --pfam-ids; a directory (or trailing "
                        "'/') keeps one cache file per Pfam ID, e.g. "
                        "data/.uniprot_cache/")
    p.add_argument("--min-length", type=int, default=0,
                   help="drop members whose domain sequence is shorter than this")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    families = helpers.load_family_metadata(args.pfam_families, args.pfam_ids)
    builder = PfamLegacyBuilder(min_length=args.min_length)
    run_pfam_join(builder, args.input, args.pfam_ids,
                    families=families, uniprot_paths=args.uniprot,
                    output=args.output, join=args.join, gzip=args.gzip,
                    description=args.description, cache_path=args.uniprot_cache)


# ===========================================================================
# Field extraction & caption helpers
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
    """GO cross-ref term names, aspect-prefix stripped, ordered C -> F -> P.

    Within an aspect the source order is preserved (stable sort), matching the
    legacy Pfam caption's grouped GENE ONTOLOGY field.
    """
    items = []
    for line in rec.get("cross_references", {}).get("GO", []):
        parts = line.split(";")
        if len(parts) < 3:
            continue
        desc = parts[2].strip().rstrip(".")
        aspect = None
        if len(desc) > 2 and desc[1] == ":":  # drop C:/F:/P: aspect prefix
            aspect, desc = desc[0], desc[2:]
        items.append((aspect, desc))
    items.sort(key=lambda t: _ASPECT_ORDER.get(t[0], 3))
    return [d for _, d in items]


def extract_fields(rec) -> dict:
    """Map a matched UniProt record to the legacy Pfam annotation field set.

    CC-comment fields are collected as **lists of all blocks** of that topic
    (not just the first), evidence-stripped. ``protein_name`` is a single
    string; ``gene_ontology`` (aspect-ordered) and ``lineage`` are lists.
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


def _caption_value(key, value) -> str:
    """Render a UniProt field value to its caption string ("" when absent)."""
    if not value:
        return ""
    if key == "lineage":
        return "The organism lineage is " + ", ".join(value)
    if key in ("gene_ontology", "cofactor"):
        return ", ".join(value)
    if isinstance(value, list):
        return " ".join(value)        # join all blocks of a CC topic
    return value


def _uniprot_section(fields: dict) -> str:
    """The UniProt-derived caption section: ``"LABEL: value"`` parts joined by
    ". ", in SPEC order. ``_FORCED`` fields are emitted even when empty (the
    legacy always writes PROTEIN NAME / GENE ONTOLOGY / LINEAGE once a member
    resolves to a UniProt entry)."""
    parts = []
    for label, key in SPEC:
        val = _strip_pubmed(_strip_evidence(_caption_value(key, fields.get(key)))).rstrip(".")
        if val or key in _FORCED:
            parts.append(f"{label}: {val}")
    section = ". ".join(parts)
    if section and not section.endswith("."):
        section += "."
    return section


def compose_caption(family: dict, fields: dict, *, has_uniprot: bool) -> str:
    """Compose the legacy Pfam caption: a FAMILY NAME / FAMILY DESCRIPTION
    section (always), with the UniProt-derived section appended directly when
    *has_uniprot* (no separator — the description already ends its own sentence,
    mirroring the legacy output)."""
    name = helpers.family_name(family) or ""
    desc = helpers.family_description(family) or ""
    caption = f"FAMILY NAME: {name}. FAMILY DESCRIPTION: {desc}"
    if has_uniprot:
        caption += _uniprot_section(fields)
    return caption


if __name__ == "__main__":
    main()
