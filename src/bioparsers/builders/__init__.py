"""A small framework for composing curated datasets from parsed UniProt JSONL.

This subpackage is a sibling to ``bioparsers.parsers`` and consumes its
output. Where the parser layer turns raw ``.dat`` bytes into faithful
``Record`` dicts, the builder layer turns a stream of those record dicts
into a curated dataset — selecting, cleaning, and projecting fields.

It provides only the *framework*:

- :class:`Builder` — the abstract base every composition subclasses. Each
  concrete builder declares a stable ``name`` (e.g. ``swissprot_legacy``)
  and a long-form ``description`` documenting its output record form, both
  enforced at definition time.
- :func:`load_jsonl` / :func:`write_jsonl` / :func:`jsonl_writer` /
  :func:`materialize` — streaming, gzip-aware JSONL I/O.
This framework is database-agnostic. The *record-shaped* logic — helpers
that read specific fields and filters that select on them — is
database-specific and lives in a per-database subpackage instead:
``bioparsers.builders.uniprot`` (``helpers``, ``filters``). A new source
database gets its own sibling subpackage (e.g. ``builders.pfam``).

Concrete builders are **not** defined here either — they live in the
``recipes/`` scripts, each demonstrating how to define a custom
:class:`Builder` and use it (with a database's helpers/filters and this
framework's io) to compose a dataset::

    from bioparsers.builders import Builder, load_jsonl, write_jsonl
    from bioparsers.builders.uniprot import helpers

    class MyBuilder(Builder):
        '''{accession, sequence, function} records.'''
        name = "my_v1"
        def build(self, records):
            for rec in records:
                fn = helpers.joined_comment(rec, "FUNCTION")
                if fn:
                    yield {"accession": rec["primary_accession"],
                           "sequence": rec["sequence"], "function": fn}

    write_jsonl(MyBuilder().build(load_jsonl("data/uniprot_sprot.jsonl")),
                "outputs/my.jsonl")
"""

from bioparsers.builders.base import Builder
from bioparsers.builders.concatenate import ConcatenatedDataset, concatenate
from bioparsers.builders.io import (
    iter_text_lines,
    jsonl_writer,
    load_jsonl,
    materialize,
    write_jsonl,
)
from bioparsers.builders.manifest import Stage, generate_manifest, write_manifest

__all__ = [
    "Builder",
    "load_jsonl",
    "write_jsonl",
    "jsonl_writer",
    "iter_text_lines",
    "materialize",
    "concatenate",
    "ConcatenatedDataset",
    "Stage",
    "generate_manifest",
    "write_manifest",
]
