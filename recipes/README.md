# recipes

Runnable recipe scripts that turn **parsed reference-database JSONL** into
curated datasets with a `Builder`. The Swiss-Prot recipes are a single streaming
pass over UniProt selected by Pfam domain; the Pfam recipes join Pfam members,
family metadata, and UniProt annotation; the supplement recipe is a flat per-row
transform of the parsed supplement table (no Pfam filter, no join). The shared
runners and field logic live in the package under `bioparsers.builders.uniprot`
(`run_by_pfam`) and `bioparsers.builders.pfam` (`run_pfam_join`); the
database-agnostic framework lives in `bioparsers.builders`.

(`swissprot` here means the Swiss-Prot — reviewed — section of UniProt, as
opposed to TrEMBL.)

| Recipe | Builder | Output record |
|---|---|---|
| `build_swissprot_legacy_by_pfam.py` | `swissprot_legacy` | `{accession, sequence, pfam_ids, caption, fields}` — a legacy-style `[final]text_caption` plus the fields it is built from |
| `build_swissprot_caption_fields_by_pfam.py` | `swissprot_caption_fields` | `{accession, sequence, pfam_ids, fields, caption_fields}` — raw fields and a parallel dict of cleaned, concatenated per-field text (no assembled caption) |
| `build_swissprot_demo_fields_by_pfam.py` | `swissprot_demo_fields` | `{accession, sequence, fields:{name?, function?, domains?}}` — a minimal demo |
| `build_pfam_legacy_by_pfam.py` | `pfam_legacy` | `{accession, sequence, region, pfam_ids, caption, fields}` — a legacy-style Pfam `[final]text_caption` (FAMILY NAME/DESCRIPTION + UniProt annotation joined on the member accession) plus the fields it is built from |
| `build_pfam_caption_fields_by_pfam.py` | `pfam_caption_fields` | `{accession, sequence, region, pfam_ids, fields, caption_fields}` — raw fields and a parallel dict of cleaned, concatenated per-field text (no assembled caption) |
| `build_supplement_legacy.py` | `supplement_legacy` | `{accession, sequence, pfam_ids, caption, fields}` — a legacy-style Supplemental `[final]text_caption` (PROTEIN NAME / LINEAGE + optional SH3 PARALOG NAME / PARALOG FUNCTION) plus the fields it is built from |
| `build_supplement_caption_fields.py` | `supplement_caption_fields` | `{accession, sequence, pfam_ids, fields, caption_fields}` — raw fields + bare caption-ready text (no assembled caption) |
| `build_legacy_dataset.py` | `legacy_sh3_dataset` | a 4-column **CSV** (`primary_Accession, protein_sequence, [final]text_caption, pfam_label`) concatenating the three section outputs in legacy order |
| `concatenate_datasets.py` | `concatenated_dataset` | concatenate any number of labeled `NAME=PATH` JSONL[.gz] sources into one JSONL, each record tagged with a root `source` |
| `fetch_uniprot_annotations.py` | *(fetch stage)* | one parsed `UniProtRecord` per line, retrieved by accession from the UniProt REST API — the accessions read from a UniProt-style FASTA |

Optional text fields are omitted when the source entry has no value.

## Getting annotation for a FASTA (the fetch recipe)

`fetch_uniprot_annotations.py` is the odd one out: it *produces* parsed
records rather than consuming them, and is the entry point when you have
sequences but no annotation — a jackhmmer/hmmsearch hit set, an alignment
export, any FASTA that kept its UniProt headers. It reads the accessions out of
the headers, retrieves each entry, and writes the same JSONL a local
`bioparsers uniprot` parse would.

```bash
python recipes/fetch_uniprot_annotations.py data/sod1s_final.fasta \
    -o outputs/sod1s_final_annotations.jsonl \
    --save-flat outputs/sod1s_final.dat.gz --summary
```

Its output feeds the builder recipes directly — it is parsed UniProt JSONL:

```bash
python recipes/build_swissprot_legacy_by_pfam.py \
    outputs/sod1s_final_annotations.jsonl \
    --pfam-ids PF00080 --pfam-names data/pfam_names.tsv \
    --include-unreviewed --join -o outputs/sod1s_final_legacy.jsonl
```

Prefer `--save-flat` for anything you keep: it stores the raw flat file the API
served, which reparses offline (`bioparsers uniprot <path>`) to identical
records — so a parser fix costs a reparse, not a refetch. The manifest records
which UniProt release answered and every accession that yielded no record,
split into `invalid_format` and `not_returned`; neither is recoverable from the
JSONL alone.

Note `--include-unreviewed` above. The Swiss-Prot recipes are reviewed-only by
default, matching the legacy dataset's Swiss-Prot section; a fetched hit set is
usually mostly TrEMBL, so it needs the opt-in. Field extraction is identical for
both sections (`PROTEIN NAME` falls back to `SubName`), but TrEMBL annotation is
automatic rather than curated — `--summary` reports per-field coverage split by
section so the difference is visible up front rather than discovered later.

## The supplement recipe

`build_supplement_legacy.py` reproduces the **Supplemental** section. Its input
is the parsed supplement table — first convert the CSV to JSONL with the `csv`
parser, then build (no Pfam filter, no UniProt join; `pfam_label` is empty for
every supplement entry):

```bash
bioparsers csv databases/misc/SH3_supplement_data.csv -o data/supplement.jsonl
python recipes/build_supplement_legacy.py data/supplement.jsonl \
    -o outputs/supplement_legacy_SH3.jsonl
```

Because the supplement table is its own source (no external release to drift
against), this reproduces the section's captions exactly.

## A deliberate deviation from legacy: COFACTOR

COFACTOR is the one caption-relevant CC topic whose text is a **structured
record** rather than prose:

```
Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence={ECO:...};
```

The four caption recipes take it via `helpers.cofactor_names()`, which keeps
only the `Name=` values and drops the ChEBI cross-reference, evidence codes and
`Note=` prose — matching legacy, which renders `COFACTOR: Mg(2+).` (verified
against all 23 COFACTOR-bearing rows of `FINAL_SH3_swissprot.csv`). Taking the
block verbatim, as the prose topics are taken, dragged the whole structured
record into the caption.

Where an entry declares **several** cofactors, all are kept as a list and
comma-joined in the caption (`COFACTOR: Cu cation, Zn(2+).`). Legacy has no
multi-cofactor row to match against — the three SH3 entries that now differ
(`P41239`, `G5EE56`, `O45539`) each gained a second cofactor in UniProt since,
and legacy kept only the first. This follows the same principle as DOMAIN
NAMES: capture everything in list form and filter downstream, rather than
bake in a first-only truncation.

## Assembling the complete legacy dataset

`build_legacy_dataset.py` stitches the three section outputs into the single
4-column CSV of `FINAL_SH3_all_dataset_with_prompts.csv`, concatenated in legacy
order (**Supplemental → Swiss-Prot → Pfam**). Build each section first, then
assemble:

```bash
python recipes/build_legacy_dataset.py \
    --supplement outputs/supplement_legacy_SH3.jsonl \
    --swissprot  outputs/swissprot_legacy_SH3.PF00018.jsonl \
    --pfam       outputs/pfam_legacy_SH3.PF00018.jsonl \
    -o outputs/FINAL_SH3_all_reproduced.csv
```

`pfam_label` is filled per the legacy convention (empty for supplement, the
list-repr of all Pfam IDs for Swiss-Prot, the single family accession for Pfam).
The Supplemental section reproduces exactly; the Swiss-Prot and Pfam sections are
approximate (Pfam-release drift), so the assembled file matches the legacy CSV's
form and ordering but not its exact row count.

## Concatenating datasets (the combined caption_fields dataset)

`concatenate_datasets.py` is a generic utility (core in
`bioparsers.builders.concatenate`): it concatenates any number of labeled
`NAME=PATH` JSONL[.gz] sources, in order, adding the source `NAME` at the root of
each record (key `source`, override with `--source-key`). It is the improved
counterpart to the legacy CSV — instead of an assembled caption, combine the
three sources' **caption_fields** outputs into one source-tagged JSONL so a
trainer can compose captions on the fly (annotation dropout, field
randomization). Build each section's `caption_fields` first, then concatenate:

```bash
python recipes/concatenate_datasets.py \
    supplemental=outputs/supplement_caption_fields_SH3.jsonl \
    swissprot=outputs/swissprot_caption_fields_SH3.PF00018.jsonl \
    pfam=outputs/pfam_caption_fields_SH3.PF00018.jsonl \
    -o outputs/SH3_caption_fields_all.jsonl
```

Each output record is the input record with `source` added at the root, e.g.
`{source, accession, sequence, region?, pfam_ids, fields, caption_fields}`. The
sources are arbitrary — any labeled JSONL files concatenate the same way.

## The Pfam recipes (a three-way join)

The two Pfam recipes share one join runner (`run_pfam_join`):
`build_pfam_legacy_by_pfam.py` assembles the legacy `[final]text_caption`, and
`build_pfam_caption_fields_by_pfam.py` instead keeps the structured `fields` and
cleaned `caption_fields` text side by side (no assembled caption) for on-the-fly
caption composition — the Pfam analogue of the `swissprot_caption_fields`
recipe. Both reproduce the **Pfam** section of the legacy BioM3 SH3 dataset.
Unlike the Swiss-Prot recipes, the input is the parsed Pfam **member FASTA**
JSONL (one redundancy-reduced domain sequence per member — the row's
`protein_sequence`), and each joins two more sources:

- `--pfam-families` — a parsed Pfam **family** JSONL (`bioparsers pfam`),
  supplying `FAMILY NAME` / `FAMILY DESCRIPTION` (the same for every member of a
  family).
- `--uniprot` — one or more parsed UniProt JSONL files, joined on the member
  accession for `PROTEIN NAME … GENE ONTOLOGY, LINEAGE`. Most Pfam members are
  TrEMBL, so the default is Swiss-Prot **then** TrEMBL
  (`data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz`); a member with
  no accession / no UniProt match keeps only the FAMILY section. The join uses
  an in-memory accession index (the ~24k accessions of a family fit easily),
  streaming the UniProt file(s) once and stopping early once all are found.

The caption is `FAMILY NAME: … . FAMILY DESCRIPTION: …` with the UniProt-derived
section (`GENE ONTOLOGY` ordered by aspect C→F→P) appended directly. As with the
Swiss-Prot recipes this is an *approximate* reproduction — the published dataset
was built against an older Pfam release, so the redundancy-reduced member set
has drifted.

### Runtime

The member FASTA and (especially) UniProt are large, so the join is optimized:

- Both scans **prefilter each line** (a cheap accession/family substring check)
  and only `json.loads` matches, and decompress `.gz` with **`pigz`** when the
  binary is on `PATH` (falling back to stdlib `gzip`) — so the cost is roughly
  decompression-bound rather than JSON-bound.
- The member scan **stops** once the requested families' (contiguous) blocks
  have passed.
- `--uniprot-cache PATH` writes the resolved UniProt subset (gzipped, with the
  requested accession set + sources it was built from). A later run whose
  accessions are a **subset** of the cached set, with the same `--uniprot`
  sources, reuses it and **skips the UniProt scan entirely** — e.g. running this
  recipe and then the `caption_fields` variant on the same family pays the
  TrEMBL scan only once. Point it at a **plain file** to cache the union of all
  `--pfam-ids`, or at a **directory** (trailing `/`, e.g. `data/.uniprot_cache/`)
  to keep one cache file per Pfam ID — then each family is independently
  reusable, and a run mixing cached and new families scans only the new ones.
  A per-family cache is small (only the *resolved* records — ~10–20 MB gzipped
  for SH3); a global UniProt index would instead be ~5–10 GB up to ~1 TB.

```bash
python recipes/build_pfam_legacy_by_pfam.py data/pfam_fasta.jsonl.gz \
    --pfam-ids PF00018 --pfam-families data/pfam.jsonl.gz \
    --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \
    --uniprot-cache data/.uniprot_cache/ \
    -o outputs/pfam_legacy_SH3.jsonl
```

The first run still pays one full TrEMBL pass (~tens of minutes — there is no
way to find scattered accessions in a 143 GB stream without reading it); every
later run on cached families is seconds.

Capture the Pfam field set plus caption-ready text (no assembled caption). With
the same `--uniprot-cache` directory it reuses the cache the legacy run built —
the TrEMBL scan is paid once across both recipes:

```bash
python recipes/build_pfam_caption_fields_by_pfam.py data/pfam_fasta.jsonl.gz \
    --pfam-ids PF00018 --pfam-families data/pfam.jsonl.gz \
    --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \
    --uniprot-cache data/.uniprot_cache/ \
    -o outputs/pfam_caption_fields_SH3.jsonl
```

## Pfam filtering

Pass one or more Pfam accessions with `--pfam-ids`. The output mode depends
on `--join`:

- **default (per ID):** one output file per Pfam ID. The ID is inserted
  before the extension, e.g. `-o outputs/sprot.jsonl` with
  `--pfam-ids PF00069 PF00027` writes `sprot.PF00069.jsonl` and
  `sprot.PF00027.jsonl`. An entry carrying several of the requested domains
  appears in each of their files. Made in a single streaming pass.
- **`--join`:** one output file containing the **union** of entries matching
  *any* of the IDs, each written once (no duplication).

```bash
# one file per domain
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00069 PF00027 -o outputs/sprot.jsonl

# union into a single file
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00069 PF00027 --join -o outputs/sprot_kinases.jsonl
```

For the Pfam recipe `--pfam-ids` selects which **families** to emit (members of
those families); the same per-ID vs `--join` output rules apply.

Inputs and outputs may be gzipped (`.gz` input is auto-detected; pass `--gzip`
to compress the output). `--min-length` applies on top of the Pfam filter for
all recipes (for `pfam_legacy` it measures the **domain region**); the demo
recipe also takes `--reviewed-only`. The `swissprot_legacy` and
`swissprot_caption_fields` recipes additionally **require** `--pfam-names`
(a `PF<TAB>name` TSV written by `scripts/parse_pfam_names.sh`) to fill the
FAMILY NAMES / `family_names` field; the `pfam_legacy` recipe instead
**requires** `--pfam-families` (a parsed Pfam family JSONL) and joins
`--uniprot` file(s) on the member accession.

## Build manifests

Every output file gets a `<output>.manifest.json` reproducibility sidecar
recording the bioparsers version + git state, the builder's name and
description, the environment, the Pfam IDs, and the record count. Pass
`--description "..."` to add a free-text note to each manifest.

## Examples

Reproduce the Swiss-Prot portion of the BioM3 legacy SH3 dataset:

```bash
python recipes/build_swissprot_legacy_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_legacy_SH3.jsonl
```

Capture the full field set plus caption-ready text (no assembled caption):

```bash
python recipes/build_swissprot_caption_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-names data/pfam_names.tsv \
    -o outputs/swissprot_caption_fields_SH3.jsonl
```

A minimal demo (no `--pfam-names`):

```bash
python recipes/build_swissprot_demo_fields_by_pfam.py data/uniprot_sprot.jsonl.gz \
    --pfam-ids PF00018 -o outputs/swissprot_demo.jsonl
```

Reproduce the Pfam portion of the BioM3 legacy SH3 dataset (member FASTA in,
family metadata + UniProt joined on the member accession). Point both recipes at
the same `--uniprot-cache` directory so the one-time TrEMBL scan is shared — the
first call builds `data/.uniprot_cache/PF00018.jsonl.gz`, the second reuses it:

```bash
# legacy [final]text_caption
python recipes/build_pfam_legacy_by_pfam.py data/pfam_fasta.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-families data/pfam.jsonl.gz \
    --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \
    --uniprot-cache data/.uniprot_cache/ \
    -o outputs/pfam_legacy_SH3.jsonl

# structured fields + caption-ready text, no assembled caption (reuses the cache)
python recipes/build_pfam_caption_fields_by_pfam.py data/pfam_fasta.jsonl.gz \
    --pfam-ids PF00018 \
    --pfam-families data/pfam.jsonl.gz \
    --uniprot data/uniprot_sprot.jsonl.gz data/uniprot_trembl.jsonl.gz \
    --uniprot-cache data/.uniprot_cache/ \
    -o outputs/pfam_caption_fields_SH3.jsonl
```

Reproduce the Supplemental section from the parsed supplement table
(`data/supplement.jsonl`):

```bash
# legacy [final]text_caption
python recipes/build_supplement_legacy.py data/supplement.jsonl \
    -o outputs/supplement_legacy_SH3.jsonl

# structured fields + caption-ready text, no assembled caption
python recipes/build_supplement_caption_fields.py data/supplement.jsonl \
    -o outputs/supplement_caption_fields_SH3.jsonl
```

Assemble the three legacy sections into the complete dataset CSV (in order
Supplemental → Swiss-Prot → Pfam):

```bash
python recipes/build_legacy_dataset.py \
    --supplement outputs/supplement_legacy_SH3.jsonl \
    --swissprot  outputs/swissprot_legacy_SH3.PF00018.jsonl \
    --pfam       outputs/pfam_legacy_SH3.PF00018.jsonl \
    -o outputs/FINAL_SH3_all_reproduced.csv
```

Concatenate the three sources' `caption_fields` into one JSONL, each record
tagged with a root `source`:

```bash
python recipes/concatenate_datasets.py \
    supplemental=outputs/supplement_caption_fields_SH3.jsonl \
    swissprot=outputs/swissprot_caption_fields_SH3.PF00018.jsonl \
    pfam=outputs/pfam_caption_fields_SH3.PF00018.jsonl \
    -o outputs/SH3_caption_fields_all.jsonl
```
