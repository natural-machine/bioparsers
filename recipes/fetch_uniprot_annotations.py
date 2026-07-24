#!/usr/bin/env python
"""Recover full UniProtKB annotation for the entries named in a FASTA file.

Reads UniProt-style FASTA headers (``>sp|P00441|SODC_HUMAN ...``), collects the
accessions, retrieves each entry from the UniProt REST API, and writes the
parsed records as JSONL — one complete ``UniProtRecord`` per line, the same
shape ``bioparsers uniprot`` produces from a local ``.dat``.

The point is to turn a bare hit set — a jackhmmer/hmmsearch output, an
alignment export, any curated FASTA that kept its headers — back into
annotated records without scanning a distribution file. Swiss-Prot and TrEMBL
accessions may be mixed; each record's ``status`` field says which section it
came from.

    python recipes/fetch_uniprot_annotations.py sod1s_final.fasta \\
        -o outputs/sod1_annotations.jsonl

Temper expectations for TrEMBL. Unreviewed entries carry ``SubName`` instead of
``RecName``, no curated ``CC`` topics (``FUNCTION``, ``SUBUNIT``, ``DOMAIN``,
...), and mostly automatic cross-references. A hit set that is 97% TrEMBL will
be correspondingly thin in exactly the fields the caption builders consume —
``--summary`` reports the per-field coverage so this is visible rather than
discovered later.

Options
-------
``-o/--output``      JSONL destination (default stdout); ``-z`` to gzip.
``--save-flat PATH`` Also keep the raw UniProt flat file the API returned.
``--summary``        Per-field annotation coverage, split Reviewed/Unreviewed.
``--missing PATH``   Write accessions the API did not return, one per line.
``--batch-size``     Accessions per request (default/maximum 100).

``--save-flat`` is worth using for anything you intend to keep. The JSONL is
this package's *reading* of the entries; the flat file is what UniProt
actually served, and it reparses offline to the same records::

    bioparsers uniprot outputs/sod1s_final.dat.gz -o rebuilt.jsonl

So a parser fix costs a reparse rather than a refetch, and the bytes stay
pinned to the day of the request even though the API tracks the current
release.

Accessions the API does not return — deleted or merged since the FASTA was
made — are reported on stderr and, with ``--missing``, listed to a file. They
are never silently dropped.

With ``-o``, a ``<output>.manifest.json`` sidecar is written alongside, in the
same convention the builder recipes use. It records the two facts a fetch
output cannot carry itself: the **UniProt release** that served the request
(the API tracks the current one, so the same command months later may return
different annotation), and the **accessions that yielded no record**, split by
cause — absent from the JSONL by definition, so invisible unless written down.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import sys

from bioparsers.builders import Stage, write_manifest
from bioparsers.fetch import uniprot_rest
from bioparsers.parsers import ParseError, dump_jsonl
from bioparsers.parsers import uniprot_fasta

#: CC topics the caption builders draw on, checked for coverage by --summary.
_CAPTION_TOPICS = (
    "FUNCTION", "CATALYTIC ACTIVITY", "SUBUNIT", "DOMAIN", "SIMILARITY",
    "SUBCELLULAR LOCATION", "PTM", "TISSUE SPECIFICITY", "COFACTOR", "PATHWAY",
)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="fetch_uniprot_annotations.py",
        description=__doc__.split("\n\n")[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("fasta", help="UniProt-style FASTA (plain or gzipped)")
    parser.add_argument("-o", "--output", default=None,
                        help="JSONL output path (default: stdout)")
    parser.add_argument("-z", "--gzip", action="store_true",
                        help="gzip-compress the output")
    parser.add_argument("--save-flat", default=None, metavar="PATH",
                        help="also write the raw UniProt flat file the API "
                             "returned (gzipped when PATH ends in .gz). "
                             "Reparses offline with 'bioparsers uniprot PATH' "
                             "to the same records, so a later parser fix needs "
                             "no refetch, and pins what the API served today")
    parser.add_argument("--missing", default=None, metavar="PATH",
                        help="write unreturned accessions here, one per line")
    parser.add_argument("--summary", action="store_true",
                        help="report per-field annotation coverage to stderr")
    parser.add_argument("--batch-size", type=int,
                        default=uniprot_rest.MAX_BATCH, metavar="N",
                        help=f"accessions per request "
                             f"(default/max {uniprot_rest.MAX_BATCH})")
    parser.add_argument("--description", default=None,
                        help="free-text note recorded in the manifest")
    return parser.parse_args(argv)


#: Manifest producer identity for this stage — the fetch is not a Builder, but
#: its provenance is exactly as worth recording (arguably more, since the API
#: tracks a moving release while a local .dat is frozen).
_FETCH_STAGE = Stage(
    "uniprot_rest_fetch",
    "UniProtKB entries retrieved by accession from the UniProt REST API, the "
    "accessions read from a UniProt-style FASTA. Output is one parsed "
    "UniProtRecord per line, identical in shape to a local .dat parse.",
)


def _provenance(args, accessions, unique, count, releases, invalid, missing) -> dict:
    """Assemble the fetch-specific manifest keys.

    Two things here are unrecoverable from the output alone and are the reason
    this manifest exists: the **UniProt release** that answered (the API serves
    the current one, so an identical rerun later can return different
    annotation), and the **accessions that produced no record** — absent from
    the output by construction, hence invisible without being written down.
    ``unresolved`` splits those by cause, since malformed and deleted mean
    different things about the input.
    """
    return {
        "source": {
            "fasta": args.fasta,
            "headers": len(accessions),
            "unique_accessions": len(unique),
        },
        "fetch": {
            "endpoint": uniprot_rest.ENDPOINT,
            "batch_size": min(args.batch_size, uniprot_rest.MAX_BATCH),
            "uniprot_releases": releases,
        },
        "outputs": {
            "jsonl": args.output,
            "flat_file": args.save_flat,
            "unresolved_list": args.missing,
        },
        "counts": {
            "requested": len(unique),
            "retrieved": count,
            "unresolved": len(invalid) + len(missing),
        },
        "unresolved": {
            "invalid_format": sorted(invalid),
            "not_returned": sorted(missing),
        },
    }


def _preview(accessions: list[str], limit: int = 10) -> str:
    """First *limit* accessions, comma-joined, with an ellipsis if truncated."""
    shown = ", ".join(accessions[:limit])
    return f"{shown} ..." if len(accessions) > limit else shown


def _open_output(path: str | None, compress: bool):
    if path is None:
        if compress:
            return gzip.open(sys.stdout.buffer, "wt", encoding="utf-8")
        return contextlib.nullcontext(sys.stdout)
    opener = gzip.open if compress else open
    return opener(path, "wt", encoding="utf-8")


def _summarize(records) -> None:
    """Yield through *records*, tallying annotation coverage, and print the
    tally to stderr once the stream is exhausted.

    A generator rather than a second pass so the summary costs no extra
    memory beyond the counters — records stream to disk as they arrive.
    """
    seen = {"Reviewed": 0, "Unreviewed": 0}
    topics = {"Reviewed": {}, "Unreviewed": {}}
    named = {"Reviewed": 0, "Unreviewed": 0}

    for record in records:
        status = record["status"]
        if status not in seen:
            seen[status] = 0
            topics[status] = {}
            named[status] = 0
        seen[status] += 1
        description = record["description"]
        if description["rec_name"] or description["sub_name"]:
            named[status] += 1
        present = {c["topic"] for c in record["comments"]}
        for topic in present & set(_CAPTION_TOPICS):
            topics[status][topic] = topics[status].get(topic, 0) + 1
        yield record

    total = sum(seen.values())
    if not total:
        return
    print("\nannotation coverage", file=sys.stderr)
    for status, count in seen.items():
        if not count:
            continue
        print(f"\n  {status}: {count} entries", file=sys.stderr)
        print(f"    {'PROTEIN NAME':<22} {named[status]:>6}  "
              f"{named[status] / count:6.1%}", file=sys.stderr)
        for topic in _CAPTION_TOPICS:
            n = topics[status].get(topic, 0)
            print(f"    {topic:<22} {n:>6}  {n / count:6.1%}", file=sys.stderr)


def main(argv=None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    try:
        accessions = list(uniprot_fasta.iter_accessions(args.fasta))
    except ParseError as exc:
        print(f"fetch_uniprot_annotations: {exc}", file=sys.stderr)
        return 1

    unique = list(dict.fromkeys(accessions))
    print(f"{len(accessions)} headers, {len(unique)} unique accessions",
          file=sys.stderr)
    if not unique:
        print("fetch_uniprot_annotations: no accessions found", file=sys.stderr)
        return 1

    missing: list[str] = []
    invalid: list[str] = []
    releases: list[dict] = []

    def note_batch(index, count):
        print(f"  batch {index}/{count}", file=sys.stderr)

    try:
        with contextlib.ExitStack() as stack:
            raw_sink = None
            if args.save_flat:
                opener = gzip.open if args.save_flat.endswith(".gz") else open
                raw_sink = stack.enter_context(
                    opener(args.save_flat, "wt", encoding="utf-8")
                ).write

            records = uniprot_rest.iter_records(
                unique,
                batch_size=args.batch_size,
                on_missing=missing.extend,
                on_invalid=invalid.extend,
                on_batch=note_batch,
                on_release=releases.append,
                raw_sink=raw_sink,
            )
            if args.summary:
                records = _summarize(records)

            handle = stack.enter_context(_open_output(args.output, args.gzip))
            count = dump_jsonl(records, handle)
    except (uniprot_rest.FetchError, ParseError) as exc:
        print(f"fetch_uniprot_annotations: {exc}", file=sys.stderr)
        return 1

    print(f"\n{count} records retrieved for {len(unique)} accessions",
          file=sys.stderr)
    if args.save_flat:
        print(f"raw flat file written to {args.save_flat}", file=sys.stderr)
    if invalid:
        print(f"{len(invalid)} identifiers are not valid UniProt accessions "
              f"and were not requested: {_preview(invalid)}", file=sys.stderr)
    if missing:
        print(f"{len(missing)} accessions not returned by the API "
              f"(deleted or never assigned): {_preview(missing)}",
              file=sys.stderr)
    unresolved = invalid + missing
    if unresolved and args.missing:
        with open(args.missing, "w") as handle:
            handle.write("\n".join(unresolved) + "\n")
        print(f"  {len(unresolved)} unresolved accessions written to "
              f"{args.missing}", file=sys.stderr)

    if args.output:
        path = write_manifest(
            _FETCH_STAGE, args.output + ".manifest.json",
            description=args.description, output=args.output,
            record_count=count,
            extra=_provenance(args, accessions, unique, count,
                              releases, invalid, missing),
        )
        print(f"manifest: {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
