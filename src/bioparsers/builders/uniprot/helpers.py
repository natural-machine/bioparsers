"""Field-level helpers shared by dataset builders.

Pure functions operating on a parsed UniProt record dict (the
``as_dict()`` form emitted by the parser's JSONL) or on individual field
values. No I/O. Two artifacts of the source format that callers almost
always want to clean up:

- **Evidence tags** ``{ECO:...}`` ride along inside free-text fields such
  as comment ``text`` and ``keywords`` (the parser keeps them faithfully).
  :func:`strip_evidence` removes them.
- **Wrapping artifacts** — collapsed double spaces left by the parser's
  line joins. :func:`clean_text` normalizes whitespace too.

Note: structured ``description`` name values (e.g. ``rec_name.full``) are
already evidence-free — the parser separates their ``{ECO:...}`` into a
sibling ``evidence`` list — so :func:`full_name` needs no stripping.
"""

from __future__ import annotations

import re

#: One ``{ECO:...}`` evidence tag, with any single space in front of it so
#: ``"Heme {ECO:0000256|ARBA:ARBA00022617}"`` -> ``"Heme"``.
_EVIDENCE_RE = re.compile(r"\s*\{ECO:[^}]*\}")
_WS_RE = re.compile(r"\s+")
_SPACE_BEFORE_PUNCT_RE = re.compile(r"\s+([.,;:)])")
_REPEAT_PERIOD_RE = re.compile(r"\.{2,}")

#: A single inline citation reference: a PubMed ID or a ``Ref.N`` pointer.
_CITATION_TOKEN = r"(?:PubMed:\d+|Ref\.\d+)"
#: A parenthetical that is *entirely* a citation list, e.g.
#: ``"(PubMed:10707082)"`` or ``"(PubMed:1, PubMed:2, Ref.3)"`` — tokens
#: joined by ``,`` or ``and``, tolerating a stray leading comma left by an
#: earlier removal. The leading ``\s*`` consumes the space before it.
_CITATION_GROUP_RE = re.compile(
    rf"\s*\(\s*,?\s*{_CITATION_TOKEN}(?:\s*(?:,|and)\s*{_CITATION_TOKEN})*\s*\)"
)


def strip_evidence(text: str) -> str:
    """Remove ``{ECO:...}`` evidence tags (and the space before them)."""
    return _EVIDENCE_RE.sub("", text)


def strip_citations(text: str) -> str:
    """Remove parenthetical citation groups (and the space before them).

    Targets parentheticals composed *entirely* of PubMed/``Ref.``
    references, e.g. ``"(PubMed:10707082)"`` or
    ``"(PubMed:1, PubMed:2, Ref.3)"``. This is deliberately conservative:
    a PubMed id embedded in descriptive parenthetical prose (e.g.
    ``"(at pH 7.5, in PubMed:15629119)"``) is left intact, since excising
    just the id would leave dangling words.
    """
    return _CITATION_GROUP_RE.sub("", text)


def clean_text(text: str) -> str:
    """Strip evidence tags and citations, then normalize residual artifacts.

    Removes ``{ECO:...}`` evidence tags and parenthetical PubMed/``Ref.``
    citation groups, collapses whitespace, and repairs the punctuation
    artifacts these removals leave behind: UniProt often writes
    ``"... sentence. {ECO:...}."``, so stripping yields a doubled period —
    collapsed here to one — and stray spaces before closing punctuation
    are removed.
    """
    s = strip_citations(strip_evidence(text))
    s = _WS_RE.sub(" ", s).strip()
    s = _SPACE_BEFORE_PUNCT_RE.sub(r"\1", s)
    s = _REPEAT_PERIOD_RE.sub(".", s)
    return s


def full_name(rec: dict) -> str | None:
    """Primary protein name: ``description.rec_name.full``, falling back to
    ``sub_name.full`` (TrEMBL), else ``None``. Already evidence-free.
    """
    desc = rec.get("description") or {}
    for key in ("rec_name", "sub_name"):
        block = desc.get(key)
        if block and block.get("full"):
            return block["full"]
    return None


def comment_texts(rec: dict, topic: str) -> list[str]:
    """The ``text`` of every comment whose ``topic`` equals *topic*.

    A UniProt entry can carry several distinct ``-!-`` comment blocks of
    the same topic (e.g. multiple FUNCTION or DOMAIN annotations); this
    returns them all, in source order, as a list — the lossless primitive
    for callers that want to keep them separate.
    """
    return [c["text"] for c in rec.get("comments", []) if c.get("topic") == topic]


def joined_comment(rec: dict, topic: str, *, sep: str = " ") -> str:
    """All comments of *topic*, each cleaned, joined by *sep* (source order).

    Handles the common case of multiple same-topic comment blocks by
    cleaning each independently (so block-boundary artifacts don't merge)
    and concatenating. Returns ``""`` when the entry has no comment of that
    topic. Use :func:`comment_texts` to keep the blocks as a list.
    """
    parts = [clean_text(p) for p in comment_texts(rec, topic)]
    return sep.join(p for p in parts if p)


#: The ``Name=`` value of one cofactor: everything up to the field-terminating
#: ``;``. A COFACTOR block is structured, not prose —
#: ``"Name=Zn(2+); Xref=ChEBI:CHEBI:29105; Evidence=...;"`` — and one block can
#: declare several cofactors in sequence.
_COFACTOR_NAME_RE = re.compile(r"Name=([^;]+)")


def cofactor_names(rec: dict) -> list[str]:
    """The cofactor names an entry declares, e.g. ``["Cu cation", "Zn(2+)"]``.

    COFACTOR is the one caption-relevant comment topic whose text is a
    structured record rather than prose, so taking it verbatim (as the other
    topics are taken) drags ChEBI cross-references and evidence codes into what
    should read as a cofactor list::

        raw   Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence={ECO:...};
        here  Mg(2+)

    This matches the legacy BioM3 captions, which render ``COFACTOR: Mg(2+).``
    — verified against all 23 COFACTOR-bearing rows of
    ``FINAL_SH3_swissprot.csv``.

    Names are collected across every COFACTOR block and every ``Name=`` within
    a block, in source order, since an entry may bind more than one cofactor
    (Cu/Zn superoxide dismutase declares both). ``Note=`` prose and all other
    subfields are dropped.
    """
    names = []
    for text in comment_texts(rec, "COFACTOR"):
        for match in _COFACTOR_NAME_RE.finditer(text):
            name = clean_text(match.group(1)).strip().rstrip(".")
            if name:
                names.append(name)
    return names


def keywords(rec: dict) -> list[str]:
    """Keyword list with ``{ECO:...}`` evidence tags stripped."""
    return [strip_evidence(k).strip() for k in rec.get("keywords", [])]


def pfam_ids(rec: dict) -> list[str]:
    """Pfam accessions (e.g. ``["PF00199"]``) from ``cross_references``.

    Each Pfam ``DR`` line is stored verbatim, e.g.
    ``"Pfam; PF00199; Catalase; 1."``; the accession is the second
    ``;``-delimited field.
    """
    out = []
    for line in rec.get("cross_references", {}).get("Pfam", []):
        fields = [f.strip() for f in line.split(";")]
        if len(fields) >= 2 and fields[1]:
            out.append(fields[1])
    return out
