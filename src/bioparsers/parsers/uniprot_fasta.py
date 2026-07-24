"""Parser for the UniProtKB FASTA format (Swiss-Prot and TrEMBL).

The FASTA distribution of UniProtKB, and any FASTA derived from it that keeps
its headers intact — BLAST/HMMER hit sets, alignment exports, curated
sub-selections. The header carries the accession, so a file of these headers
is enough to recover the full annotation for every entry: parse here, then
retrieve the flat-file records with ``bioparsers.fetch.uniprot_rest`` (or by
scanning a local ``.dat``) and parse those with
``bioparsers.parsers.uniprot_dat``.

Contract: ``iter_records(path) -> Iterator[UniProtFastaRecord]``. Captures what
the header states, with no reshaping and no lookup against any other source.
The header's ``db`` token (``sp``/``tr``) distinguishes Swiss-Prot from TrEMBL
and is also mapped to the ``Reviewed``/``Unreviewed`` word that the ``.dat``
ID line uses, so records from both parsers key alike.

The header grammar is the one documented at
https://www.uniprot.org/help/fasta-headers::

    >db|UniqueIdentifier|EntryName ProteinName OS=... OX=... [GN=...] PE=... SV=...

``ProteinName`` is free text and may itself contain ``=`` or brackets (e.g.
``Superoxide dismutase [Cu-Zn]``), so it is taken as everything up to the first
recognized ``KEY=`` token rather than by splitting on whitespace. Every
``KEY=value`` field is optional — TrEMBL entries routinely omit ``GN``, and
non-canonical exports may omit more — and each is ``None`` when absent.

UniProtFastaRecord fields (``record_type="uniprot_fasta"``; the annotations
on the class are the executable copy of this list)
--------------------------------------------------------------------
- ``db``                : str, ``"sp"`` (Swiss-Prot) or ``"tr"`` (TrEMBL)
- ``status``            : str, ``"Reviewed"`` / ``"Unreviewed"`` — the word the
                          ``.dat`` ID line carries, derived from ``db``
- ``accession``         : str, the UniProt accession (header field 2)
- ``entry_name``        : str, the ID-line mnemonic (header field 3)
- ``protein_name``      : str, free text between the entry name and the first
                          ``KEY=`` token (``""`` when the header has none)
- ``organism``          [OS] : str | None
- ``taxon_id``          [OX] : int | None, NCBI TaxID
- ``gene_name``         [GN] : str | None
- ``protein_existence`` [PE] : int | None, evidence level 1-5
- ``sequence_version``  [SV] : int | None
- ``sequence``          : str, residues with all whitespace removed
- ``sequence_length``   : int, length of ``sequence``
- ``header``            : str, the raw header verbatim, ``>`` stripped

Fail-loud (raises ``ParseError``)
---------------------------------
- compressed-stream truncation (via ``base.iter_lines``)
- a non-empty file whose first non-blank line is not a ``>`` header
- a header that is not ``db|accession|entry_name`` (fewer than 3 ``|`` fields,
  or an empty accession)
- a ``db`` token other than ``sp`` / ``tr``
- a header with no sequence following it
- a non-integer ``OX`` / ``PE`` / ``SV`` value
"""

from __future__ import annotations

import re
from typing import ClassVar, Iterator

from bioparsers.parsers.base import ParseError, Record, iter_lines

RECORD_TYPE = "uniprot_fasta"

#: ``db`` token -> the STATUS word the ``.dat`` ID line uses for that section,
#: so ``UniProtFastaRecord.status`` and ``UniProtRecord.status`` are comparable.
_DB_STATUS = {"sp": "Reviewed", "tr": "Unreviewed"}


# ===========================================================================
# Public API
# ===========================================================================

class UniProtFastaRecord(Record):
    """One UniProtKB FASTA entry. The annotations below are the single
    executable schema (Pylance-typed); ``Record.__init__`` enforces that
    ``parse_header`` plus the sequence emit exactly these keys.
    """

    record_type: ClassVar[str] = RECORD_TYPE

    db: str
    status: str
    accession: str
    entry_name: str
    protein_name: str
    organism: str | None
    taxon_id: int | None
    gene_name: str | None
    protein_existence: int | None
    sequence_version: int | None
    sequence: str
    sequence_length: int
    header: str


def iter_records(path: str) -> Iterator[UniProtFastaRecord]:
    """Yield one :class:`UniProtFastaRecord` per entry in the FASTA at *path*.

    Reads through :func:`base.iter_lines` (fail-loud on a truncated compressed
    stream), groups each header with the sequence lines that follow it, and
    delegates header parsing to :func:`parse_header`. Raises ``ParseError`` if
    the file does not start at a ``>`` header or any header has no sequence.
    """
    header: str | None = None
    chunks: list[str] = []

    for line in iter_lines(path):
        if line.startswith(">"):
            if header is not None:
                yield _build(header, chunks, path)
            header = line[1:].rstrip("\n").rstrip("\r")
            chunks = []
            continue
        if header is None:
            if line.strip() == "":
                continue
            raise ParseError(
                f"{path}: expected a '>' header at start of file, got {line!r}"
            )
        chunks.append(line)

    if header is not None:
        yield _build(header, chunks, path)


def iter_accessions(path: str) -> Iterator[str]:
    """Yield the accession of every entry in the FASTA at *path*, in file
    order and including duplicates.

    The cheap path when only the identifiers are wanted — feeding an
    annotation lookup, for instance — since it skips sequence assembly::

        accessions = list(dict.fromkeys(iter_accessions("hits.fasta")))
    """
    for line in iter_lines(path):
        if line.startswith(">"):
            yield parse_header(line[1:].rstrip("\n").rstrip("\r"))["accession"]


def parse_header(header: str) -> dict:
    """Parse one UniProt FASTA header (``>`` already stripped) into its
    fields — every key of :class:`UniProtFastaRecord` except ``sequence`` and
    ``sequence_length``.

    Raises ``ParseError`` on a header that is not ``db|accession|entry_name``,
    an unknown ``db`` token, or a non-integer ``OX`` / ``PE`` / ``SV``.
    """
    ident, _, description = header.partition(" ")
    parts = ident.split("|")
    if len(parts) < 3:
        raise ParseError(
            f"malformed UniProt FASTA header (expected 'db|accession|entry_name'): "
            f"{header!r}"
        )
    db, accession, entry_name = parts[0], parts[1].strip(), parts[2]
    if db not in _DB_STATUS:
        raise ParseError(
            f"unknown UniProt db token {db!r} (expected 'sp' or 'tr'): {header!r}"
        )
    if not accession:
        raise ParseError(f"UniProt FASTA header has an empty accession: {header!r}")

    fields = _split_keyed(description)
    return {
        "db": db,
        "status": _DB_STATUS[db],
        "accession": accession,
        "entry_name": entry_name,
        "protein_name": fields.pop("_name"),
        "organism": fields.get("OS"),
        "taxon_id": _as_int(fields.get("OX"), "OX", header),
        "gene_name": fields.get("GN"),
        "protein_existence": _as_int(fields.get("PE"), "PE", header),
        "sequence_version": _as_int(fields.get("SV"), "SV", header),
        "header": header,
    }


# ===========================================================================
# Implementation details
# ===========================================================================

#: The documented header keys. Anchored on whitespace (or string start) so an
#: ``=`` inside free-text protein names cannot be mistaken for a field break.
_KEYED_RE = re.compile(r"(?:(?<=\s)|^)(OS|OX|GN|PE|SV)=")


def _split_keyed(description: str) -> dict:
    """Split a header description into ``{"_name": <free text>, KEY: value}``.

    The value of each key runs to the start of the next key (or end of string),
    which is what lets ``OS=`` hold a multi-word organism containing brackets
    and parentheses, e.g. ``Kluyveromyces lactis (strain ATCC 8585 / ...)``.
    """
    matches = list(_KEYED_RE.finditer(description))
    if not matches:
        return {"_name": description.strip()}

    out = {"_name": description[: matches[0].start()].strip()}
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(description)
        value = description[m.end() : end].strip()
        # A repeated key keeps the first occurrence — the canonical format has
        # no repeats, so this only guards against odd third-party headers.
        out.setdefault(m.group(1), value)
    return out


def _as_int(value, key: str, header: str) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError as exc:
        raise ParseError(
            f"non-integer {key}={value!r} in UniProt FASTA header: {header!r}"
        ) from exc


def _build(header: str, chunks: list[str], path: str) -> UniProtFastaRecord:
    """Assemble one record from a header and its raw sequence lines."""
    sequence = re.sub(r"\s", "", "".join(chunks))
    fields = parse_header(header)
    if not sequence:
        raise ParseError(f"{path}: {fields['accession']}: entry has no sequence")
    fields["sequence"] = sequence
    fields["sequence_length"] = len(sequence)
    return UniProtFastaRecord(**fields)
