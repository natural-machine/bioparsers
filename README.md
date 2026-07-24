# bioparsers

Parsers for biological reference databases, reading raw flat-files into
faithful, typed Python records.

## Description

`bioparsers` has three layers:

- **`parsers`** — Processes raw database inputs and produces iterables
  of structured records.
- **`fetch`** — retrieves specific entries by identifier over the network,
  yielding the same records as the corresponding parser.
- **`builders`** — a composition layer that produces curated datasets 
  from the parsed records.

Implemented:

| Database | Module | Record |
|---|---|---|
| UniProtKB Swiss-Prot / TrEMBL `.dat` | `bioparsers.parsers.uniprot_dat` | `UniProtRecord` |
| UniProtKB FASTA (`sp\|...` / `tr\|...` headers) | `bioparsers.parsers.uniprot_fasta` | `UniProtFastaRecord` |
| Pfam-A Stockholm (`Pfam-A.full`) | `bioparsers.parsers.pfam_stockholm` | `PfamRecord` |
| Pfam-A member FASTA (`Pfam-A.fasta`) | `bioparsers.parsers.pfam_fasta` | `PfamFastaRecord` |
| Delimited table (CSV / TSV) | `bioparsers.parsers.csv_table` | `CsvRecord` |
| UniProtKB REST (by accession) | `bioparsers.fetch.uniprot_rest` | `UniProtRecord` |

## Setup

Requires Python 3.12+.

With conda (in-tree env at `./env`):

```bash
conda env create -p ./env -f environment.yml
conda activate ./env
```

Or with pip into an existing environment:

```bash
pip install -e '.[dev]'
```

## Usage

### Parsers

The parser layer reads database files into `Record`s. It can be used as a library
or through the `bioparsers` console script. Each parser exposes
`iter_records(path) -> Iterator[Record]`:

- **`uniprot_dat`** — consumes UniProtKB Swiss-Prot / TrEMBL `.dat` flat files;
  yields one **`UniProtRecord`** per entry: accessions, reviewed/unreviewed
  status, names, gene names, organism + lineage + taxon, references, comments,
  features, cross-references, keywords, and the amino-acid sequence (validated
  against the ID/SQ length and CRC64).
- **`pfam_stockholm`** — consumes a Pfam-A Stockholm release (`Pfam-A.full`),
  or `Pfam-A.hmm` for just the name table; yields one **`PfamRecord`** per
  family: accession, name, description, type, clan, references, GA/TC/NC
  thresholds, cross-references, and member count — with the member list and each
  member's ungapped sequence available via opt-in.
- **`pfam_fasta`** — consumes the Pfam-A member FASTA (`Pfam-A.fasta`, the
  redundancy-reduced member set); yields one **`PfamFastaRecord`** per member
  sequence: member accession + name, aligned region, its Pfam family, and the
  ungapped residues.
- **`uniprot_fasta`** — consumes UniProtKB FASTA, or anything derived from it
  that kept the headers (BLAST/HMMER hit sets, alignment exports); yields one
  **`UniProtFastaRecord`** per entry: `sp`/`tr` section, accession, entry name,
  protein name, organism + taxon, gene, PE/SV, and the sequence. Since the
  header carries the accession, such a file is enough to recover full
  annotation — see [Fetch](#fetch).

#### Command line

The `bioparsers` console script parses a flat-file to JSONL (one compact
object per line) on stdout, or to a file with `-o`. Input may be plain
or gzipped:

```bash
bioparsers uniprot uniprot_sprot.dat.gz > out.jsonl
bioparsers uniprot in.dat -o out.jsonl
bioparsers uniprot in.dat.gz --gzip -o out.jsonl.gz   # compress output
bioparsers uniprot in.dat.gz --progress > out.jsonl   # heartbeat to stderr
bioparsers uniprot-fasta hits.fasta > hits.jsonl      # UniProt FASTA headers
```

Pass `--gzip` (`-z`) to compress the output, and `--progress [N]` for a
record-count heartbeat on stderr (every N records, default 100000). The
record count is reported on stderr; corrupt or truncated input exits
non-zero with a message on stderr.

The `pfam` and `pfam-fasta` subcommands add Pfam options. `--pfam-id`
(repeatable) restricts to given families (scanning stops once found). For
`pfam`, `--with-member-accessions` / `--with-member-sequences` opt the
per-member list into the output, and multiple `--pfam-id` write one file per
family (`pfam_<accession>.jsonl`) under the `-o` directory unless `--join` is
given:

```bash
bioparsers pfam Pfam-A.full.gz > pfam.jsonl                       # family metadata
bioparsers pfam Pfam-A.full.gz --pfam-id PF00018 --pfam-id PF07714 \
    --with-member-sequences -o out_dir/                           # one file per family
bioparsers pfam Pfam-A.full.gz --pfam-id PF00018 --join > sh3.jsonl
bioparsers pfam-fasta Pfam-A.fasta.gz --pfam-id PF00018 > sh3_members.jsonl
```

The `csv` subcommand converts a delimited table (CSV/TSV) to JSONL, one object
per row keyed by the header. `--delimiter` overrides the extension-based default
(tab for `.tsv`/`.tab`, else comma):

```bash
bioparsers csv SH3_supplement_data.csv > supplement.jsonl
bioparsers csv table.tsv > table.jsonl                 # tab inferred from .tsv
```

#### Library

```python
from bioparsers.parsers.uniprot_dat import iter_records
from bioparsers.parsers import dump_jsonl

for record in iter_records("uniprot_sprot.dat.gz"):
    print(record.primary_accession, record.organism)
    print(record.description["rec_name"])

# Stream a whole file to JSONL:
with open("out.jsonl", "w") as f:
    n = dump_jsonl(iter_records("uniprot_sprot.dat.gz"), f)
```

A `Record` is a dict-backed field-bag: access fields by attribute
(`record.sequence`) or item (`record["sequence"]`), and serialize with
`record.as_dict()` or `record.to_json()`.

The Pfam parsers work the same way. `pfam_stockholm` yields one `PfamRecord`
per family from `Pfam-A.full`; pass `accessions=` to extract only certain
families (scanning stops once they are found) and `with_member_sequences=True`
to attach each member's sequence. `pfam_fasta` yields one `PfamFastaRecord` per
member sequence from the lighter `Pfam-A.fasta`.

```python
from bioparsers.parsers import pfam_stockholm, pfam_fasta

# Family metadata + member sequences for selected families:
for fam in pfam_stockholm.iter_records("Pfam-A.full.gz", accessions=["PF00018"],
                                       with_member_sequences=True):
    print(fam.accession, fam.name, len(fam.members))

# Redundancy-reduced member sequences (one record per sequence):
for member in pfam_fasta.iter_records("Pfam-A.fasta.gz", accessions=["PF00018"]):
    print(member.accession, member.region, member.sequence)
```

`csv_table` handles sources that already ship as a structured table (e.g. the SH3
Legacy supplemental data). It is **general**: each row becomes an open-bag
`CsvRecord` keyed by the header columns, values kept verbatim as strings. The
delimiter defaults to a tab for `.tsv`/`.tab` and a comma otherwise.

```python
from bioparsers.parsers import csv_table

for row in csv_table.iter_records("SH3_supplement_data.csv"):
    print(row["primary_Accession"], row["protein_name"], row["sh3_paralog_name"])
```

### Fetch

`bioparsers.fetch` is the complement to `parsers`: where a parser reads a local
distribution file top to bottom, a **fetcher** pulls the specific entries named
by a list of identifiers. Both yield the same `Record` types, so nothing
downstream can tell which was used.

Reach for a fetcher when the wanted subset is small relative to the
distribution file. Retrieving a few thousand TrEMBL entries over the REST API
takes seconds; finding them by scanning the 160 GB local `uniprot_trembl.dat.gz`
takes hours. The trade is currency and reachability — the API serves whatever
UniProt currently publishes, so a fetch is not pinned to a frozen snapshot.

`uniprot_rest` batches accessions (100 per request, the endpoint maximum) and
asks for `format=txt`, which **is** the `.dat` flat-file format — so the payload
goes straight to the `uniprot_dat` parser and the records are ordinary
`UniProtRecord`s. Swiss-Prot and TrEMBL accessions may be mixed.

```python
from bioparsers.fetch.uniprot_rest import iter_records
from bioparsers.parsers import dump_jsonl

with open("out.jsonl", "w") as f:
    dump_jsonl(iter_records(["P00441", "A0A072PZ83"]), f)
```

Because the payload is the `.dat` format, it is worth keeping. `raw_sink` tees
each batch to disk as it arrives; the concatenation is a valid flat file that
reparses offline to the same records, so a later parser fix costs a reparse
rather than a refetch:

```python
with open("entries.dat", "w") as raw, open("out.jsonl", "w") as out:
    dump_jsonl(iter_records(accessions, raw_sink=raw.write), out)
```

Accessions that yield no record are reported, never silently dropped, through
two callbacks that distinguish the causes: `on_invalid` (not a valid UniProt
accession — excluded client-side, since one malformed id makes the API reject
its whole batch of 100) and `on_missing` (well-formed but unknown — deleted, or
never assigned). `on_release` reports which UniProt release answered.

[`recipes/fetch_uniprot_annotations.py`](recipes/fetch_uniprot_annotations.py)
wires this together: FASTA in, annotated JSONL out, with `--save-flat` for the
archival flat file, `--summary` for per-field annotation coverage, and a
manifest sidecar recording the release and every unresolved accession.

```bash
python recipes/fetch_uniprot_annotations.py hits.fasta \
    -o outputs/annotations.jsonl --save-flat outputs/entries.dat.gz --summary
```

### Builders

`bioparsers.builders` is a framework for turning parsed JSONL into
curated datasets. The `Builder` base class, streaming I/O 
(`load_jsonl` / `write_jsonl` / `jsonl_writer` / `materialize`) is 
database-agnostic. The record-shaped logic lives in a per-database
subpackage e.g. `bioparsers.builders.uniprot` (helpers/filters + the
`run_by_pfam` single-pass runner) and `bioparsers.builders.pfam` (the
`run_pfam_legacy` runner, which joins Pfam members with family metadata and
UniProt annotation). A new source database gets its own sibling subpackage.

Concrete builders are directly provided in the package, but can be defined
by users in an on-demand fashion. Each is
a `Builder` subclass with a versioned `name` and a long-form `description`
documenting its output record form. Builders are streaming-first
(constant memory); `materialize()` collects streamed results into a list.

```python
from bioparsers.builders import Builder, load_jsonl, write_jsonl
from bioparsers.builders.uniprot import helpers

class SwissProtFunction(Builder):
    """Flat {accession, sequence, function} records."""
    name = "swissprot_function_v1"
    def build(self, records):
        for rec in records:
            fn = helpers.joined_comment(rec, "FUNCTION")   # cleaned, evidence-free
            if fn:
                yield {"accession": rec["primary_accession"],
                       "sequence": rec["sequence"], "function": fn}

records = load_jsonl("data/uniprot_sprot.jsonl")          # streaming, gz-aware
n = write_jsonl(SwissProtFunction().build(records), "outputs/sprot_function.jsonl")
```

For reproducibility, `write_manifest(builder, path, ...)` writes a JSON
sidecar recording the bioparsers version + git state, the builder's name and
description, the environment, and optional run details (output path, record
count, a custom `description`):

```python
from bioparsers.builders import write_manifest
write_manifest(SwissProtFunction(), "outputs/sprot_function.jsonl.manifest.json",
               output="outputs/sprot_function.jsonl", record_count=n,
               description="flat sequence/function pairs")
```

Only `name` and `description` are read, so any step declaring those two can be
manifested the same way. `Stage` is the minimal such producer — used by the
fetch recipe, whose provenance (which UniProt release answered, which
accessions came back empty) is not recoverable from its output afterwards:

```python
from bioparsers.builders import Stage, write_manifest
write_manifest(Stage("uniprot_rest_fetch", "Entries retrieved by accession."),
               "outputs/annotations.jsonl.manifest.json", record_count=n,
               extra={"unresolved": {"invalid_format": [], "not_returned": []}})
```

The [`recipes/`](recipes/) scripts are runnable, worked examples — each defines
a `Builder` and writes a `<output>.manifest.json` sidecar per output. The
`swissprot_*` recipes **consume** parsed Swiss-Prot JSONL filtered to one or
more Pfam IDs; the `pfam_*` recipes consume parsed Pfam member FASTA JSONL and
join family metadata + UniProt annotation; the `supplement_legacy` recipe is a
flat per-row transform of the parsed supplement table. Per kept entry they
**produce**:

| Builder | Output record |
|---|---|
| `swissprot_legacy` | the sequence + a legacy-style text `caption` (and the `fields` it is built from) |
| `swissprot_caption_fields` | the sequence + raw `fields` + a `caption_fields` dict (cleaned per-field text; no assembled caption) |
| `swissprot_demo_fields` | `{accession, sequence, fields:{name?, function?, domains?}}` — a minimal demo |
| `pfam_legacy` | the domain `sequence` + `region` + a legacy-style Pfam `caption` (FAMILY NAME/DESCRIPTION + UniProt annotation joined on the member accession) and the `fields` it is built from |
| `pfam_caption_fields` | the domain `sequence` + `region` + raw `fields` + a `caption_fields` dict (cleaned per-field text; no assembled caption) |
| `supplement_legacy` | the sequence + a legacy-style Supplemental `caption` (PROTEIN NAME / LINEAGE + optional paralog fields) and the `fields` it is built from |
| `supplement_caption_fields` | the sequence + raw `fields` + a `caption_fields` dict (no assembled caption) |
| `legacy_sh3_dataset` | assembles the three section outputs into one 4-column legacy **CSV** (`build_legacy_dataset.py`), in order Supplemental → Swiss-Prot → Pfam |
| `concatenated_dataset` | concatenates any number of labeled `NAME=PATH` JSONL sources into one **JSONL** (`concatenate_datasets.py`), each record tagged with a root `source` |

Optional fields are omitted when the source has no value.

## Tests

```bash
python -m pytest tests/
```
