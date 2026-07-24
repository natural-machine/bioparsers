"""Command-line interface for bioparsers.

Each subcommand maps to one reference-database parser and streams the
parsed records as JSONL — one compact JSON object per line — to stdout
or to ``-o PATH``. Input is read transparently whether plain or gzipped,
and the tool fails loud: a truncated or corrupt input raises
``ParseError``, which is reported on stderr with a non-zero exit code
rather than producing a silently short result.

    bioparsers uniprot uniprot_sprot.dat.gz > out.jsonl
    bioparsers uniprot in.dat -o out.jsonl
    bioparsers uniprot-fasta hits.fasta > hits.jsonl
    bioparsers pfam Pfam-A.full.gz -o pfam.jsonl
    bioparsers pfam-fasta Pfam-A.fasta.gz --pfam-id PF00018 > sh3_members.jsonl
    bioparsers csv SH3_supplement_data.csv > supplement.jsonl

The ``pfam``/``pfam-fasta`` subcommands add Pfam-specific options (``--pfam-id``,
``--with-member-accessions`` / ``--with-member-sequences``, ``--join``); the
``csv`` subcommand reads a delimited table (one ``CsvRecord`` per row, keyed by
the header) and adds ``--delimiter``. See ``--help``. JSONL is the only *output*
format (mirroring the library's ``dump_jsonl`` helper); CSV/Parquet output is
deliberately out of scope.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import os
import sys
from typing import Callable, Iterator

from bioparsers.parsers import ParseError, Record, dump_jsonl, dump_jsonl_split
from bioparsers.parsers import csv_table, pfam_fasta, pfam_stockholm, uniprot_dat
from bioparsers.parsers import uniprot_fasta

#: Subcommand name -> ``iter_records`` callable for that database. Adding a
#: parser is one entry here; the subcommand and its arguments are generated.
_PARSERS: dict[str, Callable[[str], Iterator[Record]]] = {
    "uniprot": uniprot_dat.iter_records,
    "uniprot-fasta": uniprot_fasta.iter_records,
    "pfam": pfam_stockholm.iter_records,
    "pfam-fasta": pfam_fasta.iter_records,
    "csv": csv_table.iter_records,
}

#: Subcommands that accept a repeatable ``--pfam-id`` family filter.
_PFAM_FILTERABLE = {"pfam", "pfam-fasta"}


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="bioparsers",
        description="Parse biological reference databases to JSONL.",
    )
    sub = parser.add_subparsers(dest="parser", required=True, metavar="PARSER")
    for name in _PARSERS:
        p = sub.add_parser(name, help=f"parse a {name} flat-file to JSONL")
        p.add_argument(
            "input", help="path to the input file (plain or gzipped)"
        )
        p.add_argument(
            "-o", "--output", default=None,
            help="output path for the JSONL (default: stdout)",
        )
        p.add_argument(
            "-z", "--gzip", action="store_true",
            help="gzip-compress the output (stdout or -o file)",
        )
        p.add_argument(
            "--progress", nargs="?", type=int, const=100_000, default=None,
            metavar="N",
            help="print a record-count heartbeat to stderr every N records "
                 "(default 100000 when given with no value)",
        )
        if name == "csv":
            p.add_argument(
                "--delimiter", default=None, metavar="CHAR",
                help="field delimiter (default: tab for .tsv/.tab, else comma)",
            )
        if name in _PFAM_FILTERABLE:
            p.add_argument(
                "--pfam-id", action="append", default=None, metavar="PFXXXXX",
                help="restrict to this Pfam accession (repeatable); scanning "
                     "stops once all requested families are found",
            )
        if name == "pfam":
            p.add_argument(
                "--with-member-accessions", action="store_true",
                help="include the per-member list (accession/name/region); "
                     "omitted by default since num_sequences gives the count",
            )
            p.add_argument(
                "--with-member-sequences", action="store_true",
                help="attach each member's ungapped sequence, derived from the "
                     "alignment and validated against its region span (implies "
                     "--with-member-accessions)",
            )
            p.add_argument(
                "--join", action="store_true",
                help="with --pfam-id, write one unioned stream to -o/stdout; "
                     "default is one file per family, with -o naming the output "
                     "directory (pfam_<accession>.jsonl)",
            )
    return parser.parse_args(argv)


def _parser_kwargs(args) -> dict:
    """Build the parser-specific keyword arguments for ``iter_records`` from the
    parsed CLI args (only the Pfam parsers currently take any)."""
    kwargs: dict = {}
    if args.parser == "csv" and getattr(args, "delimiter", None):
        kwargs["delimiter"] = args.delimiter
    if args.parser in _PFAM_FILTERABLE and getattr(args, "pfam_id", None):
        kwargs["accessions"] = args.pfam_id
    if args.parser == "pfam":
        if getattr(args, "with_member_accessions", False):
            kwargs["with_member_accessions"] = True
        if getattr(args, "with_member_sequences", False):
            kwargs["with_member_sequences"] = True
    return kwargs


def _with_progress(records: Iterator[Record], every: int) -> Iterator[Record]:
    """Yield from *records*, printing a count heartbeat to stderr every
    *every* records. Heartbeats go to stderr so they survive a stdout pipe.
    """
    count = 0
    for count, rec in enumerate(records, 1):
        if count % every == 0:
            print(f"  ... {count} records", file=sys.stderr)
        yield rec


def _open_output(output_path: str | None, compress: bool):
    """Return a context manager yielding the text stream to write JSONL to:
    *output_path* (or stdout when None), gzip-compressed when *compress*.
    stdout is wrapped in a nullcontext so it is never closed.
    """
    if output_path is None:
        if compress:
            return gzip.open(sys.stdout.buffer, "wt", encoding="utf-8")
        return contextlib.nullcontext(sys.stdout)
    opener = gzip.open if compress else open
    return opener(output_path, "wt", encoding="utf-8")


def run(
    parser_name: str,
    input_path: str,
    output_path: str | None,
    progress: int | None = None,
    compress: bool = False,
    parser_kwargs: dict | None = None,
    split: bool = False,
) -> int:
    """Parse *input_path* with the named parser and write JSONL to
    *output_path* (or stdout when None), gzip-compressed when *compress*.
    Returns the record count.

    *parser_kwargs* are forwarded to the parser's ``iter_records`` (e.g. the
    ``pfam`` parser's ``accessions`` / ``with_member_sequences``). With *split*
    (``pfam`` extraction without ``--join``), each family is written to its own
    ``pfam_<accession>.jsonl`` file under the *output_path* directory (default
    cwd) instead of a single stream. When *progress* is a positive int, a
    heartbeat is printed to stderr every *progress* records. Propagates
    ``ParseError`` from the parser; ``main`` turns it into a non-zero exit.
    """
    iter_records = _PARSERS[parser_name]
    records = iter_records(input_path, **(parser_kwargs or {}))
    if progress:
        records = _with_progress(records, progress)
    if split:
        outdir = output_path or "."
        suffix = ".jsonl.gz" if compress else ".jsonl"
        counts = dump_jsonl_split(
            records, outdir, key=lambda r: r["accession"],
            prefix="pfam_", suffix=suffix, compress=compress,
        )
        for acc, n in counts.items():
            print(f"  {n} -> {os.path.join(outdir, f'pfam_{acc}{suffix}')}",
                  file=sys.stderr)
        return sum(counts.values())
    with _open_output(output_path, compress) as handle:
        return dump_jsonl(records, handle)


def main(argv=None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    split = bool(
        args.parser == "pfam"
        and getattr(args, "pfam_id", None)
        and not getattr(args, "join", False)
    )
    try:
        count = run(args.parser, args.input, args.output, args.progress,
                    args.gzip, _parser_kwargs(args), split=split)
    except ParseError as exc:
        print(f"bioparsers: {exc}", file=sys.stderr)
        return 1
    except BrokenPipeError:
        # A downstream consumer (e.g. `head`) closed the pipe early.
        # Redirect stdout to devnull so the interpreter's final flush
        # doesn't re-raise on exit, then report success.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0
    print(f"{count} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
