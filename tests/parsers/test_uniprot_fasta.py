"""Unit tests for the UniProtKB FASTA header parser.

The fixture ``uniprot_fasta_mini.fasta`` holds five real entries chosen to
cover the header variations that actually differ: Swiss-Prot and TrEMBL, with
and without ``GN=``, a protein name containing brackets, and an organism whose
parenthesised strain contains ``/`` and ``=``-adjacent text.
"""

import gzip
import os

import pytest

from bioparsers.parsers import ParseError
from bioparsers.parsers.uniprot_fasta import (
    UniProtFastaRecord,
    iter_accessions,
    iter_records,
    parse_header,
)

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
FASTA = os.path.join(DATDIR, "uniprot_fasta_mini.fasta")

PROMISED = (
    "db", "status", "accession", "entry_name", "protein_name", "organism",
    "taxon_id", "gene_name", "protein_existence", "sequence_version",
    "sequence", "sequence_length", "header",
)


@pytest.fixture(scope="module")
def entries():
    return list(iter_records(FASTA))


class TestStructuralInvariants:

    def test_counts_and_types(self, entries):
        assert len(entries) == 5
        assert all(isinstance(r, UniProtFastaRecord) for r in entries)
        assert all(r.record_type == "uniprot_fasta" for r in entries)

    def test_promised_schema_is_exact(self, entries):
        assert set(entries[0].as_dict()) == set(PROMISED)

    def test_sequence_length_matches_sequence(self, entries):
        assert all(r.sequence_length == len(r.sequence) for r in entries)

    def test_sequences_have_no_whitespace(self, entries):
        # The fixture wraps at 60 columns, as UniProt distributes.
        assert all(r.sequence == r.sequence.strip() for r in entries)
        assert all(" " not in r.sequence and "\n" not in r.sequence
                   for r in entries)


class TestParsing:

    def test_swissprot_entry_fields(self, entries):
        r = entries[0]
        assert r.db == "sp"
        assert r.status == "Reviewed"
        assert r.accession == "P00441"
        assert r.entry_name == "SODC_HUMAN"
        assert r.protein_name == "Superoxide dismutase [Cu-Zn]"
        assert r.organism == "Homo sapiens"
        assert r.taxon_id == 9606
        assert r.gene_name == "SOD1"
        assert r.protein_existence == 1
        assert r.sequence_version == 2
        assert r.sequence.startswith("MATKAVCVLKGDGPVQGIINFEQK")

    def test_trembl_entry_maps_to_unreviewed(self, entries):
        r = entries[3]
        assert r.db == "tr"
        assert r.status == "Unreviewed"
        assert r.accession == "A0A072PZ83"
        assert r.taxon_id == 1182545

    def test_missing_gene_name_is_none(self, entries):
        assert entries[1].accession == "P82205"
        assert entries[1].gene_name is None
        assert entries[4].gene_name is None

    def test_header_is_kept_verbatim(self, entries):
        assert entries[0].header.startswith("sp|P00441|SODC_HUMAN ")
        assert not entries[0].header.startswith(">")


class TestHeaderGrammar:
    """``parse_header`` on constructed headers — the awkward cases."""

    def test_protein_name_stops_at_first_key(self):
        f = parse_header("sp|P1AAA1|X_TEST Some protein OS=Escherichia coli OX=562")
        assert f["protein_name"] == "Some protein"
        assert f["organism"] == "Escherichia coli"

    def test_organism_may_contain_parentheses_and_slashes(self, entries):
        # Q6CPE2: "Kluyveromyces lactis (strain ATCC 8585 / CBS 2359 / ...)"
        r = entries[2]
        assert r.organism.startswith("Kluyveromyces lactis (strain ATCC 8585 /")
        assert r.organism.endswith(")")
        assert r.taxon_id == 284590

    def test_equals_inside_protein_name_is_not_a_field_break(self):
        # A bare "=" mid-name must not be mistaken for a KEY= token, and only
        # the documented keys anchor a field.
        f = parse_header("tr|A0A1AAA1|X_TEST Protein A=B thing OS=Homo sapiens OX=9606")
        assert f["protein_name"] == "Protein A=B thing"
        assert f["organism"] == "Homo sapiens"

    def test_absent_optional_fields_are_none(self):
        f = parse_header("tr|A0A1AAA1|X_TEST Uncharacterized protein")
        assert f["protein_name"] == "Uncharacterized protein"
        assert f["organism"] is None
        assert f["taxon_id"] is None
        assert f["gene_name"] is None
        assert f["protein_existence"] is None
        assert f["sequence_version"] is None

    def test_empty_protein_name(self):
        f = parse_header("sp|P1AAA1|X_TEST OS=Homo sapiens OX=9606")
        assert f["protein_name"] == ""
        assert f["organism"] == "Homo sapiens"


class TestAccessionsOnly:

    def test_iter_accessions_matches_full_parse(self, entries):
        assert list(iter_accessions(FASTA)) == [r.accession for r in entries]

    def test_iter_accessions_keeps_file_order_and_duplicates(self, tmp_path):
        p = tmp_path / "dupes.fasta"
        p.write_text(">sp|P1AAA1|A_TEST A OS=X OX=1\nMS\n"
                     ">tr|A0A1AAA1|B_TEST B OS=X OX=1\nMS\n"
                     ">sp|P1AAA1|A_TEST A OS=X OX=1\nMS\n")
        assert list(iter_accessions(str(p))) == ["P1AAA1", "A0A1AAA1", "P1AAA1"]


class TestGzipTransparency:

    def test_reads_gzipped_input(self, tmp_path, entries):
        p = tmp_path / "mini.fasta.gz"
        with open(FASTA, "rb") as src, gzip.open(p, "wb") as dst:
            dst.write(src.read())
        assert [r.as_dict() for r in iter_records(str(p))] == \
            [r.as_dict() for r in entries]


class TestFailLoud:

    def test_file_not_starting_at_a_header(self, tmp_path):
        p = tmp_path / "bad.fasta"
        p.write_text("MSEQUENCE\n>sp|P1AAA1|X_TEST A OS=X OX=1\nMS\n")
        with pytest.raises(ParseError, match="expected a '>' header"):
            list(iter_records(str(p)))

    def test_header_with_no_sequence(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text(">sp|P1AAA1|X_TEST A OS=X OX=1\n")
        with pytest.raises(ParseError, match="has no sequence"):
            list(iter_records(str(p)))

    def test_header_missing_pipe_fields(self):
        with pytest.raises(ParseError, match="malformed UniProt FASTA header"):
            parse_header("P00441 Superoxide dismutase")

    def test_unknown_db_token(self):
        with pytest.raises(ParseError, match="unknown UniProt db token"):
            parse_header("xx|P00441|SODC_HUMAN Superoxide dismutase")

    def test_empty_accession(self):
        with pytest.raises(ParseError, match="empty accession"):
            parse_header("sp||SODC_HUMAN Superoxide dismutase")

    @pytest.mark.parametrize("bad", ["OX=notanumber", "PE=high", "SV=x"])
    def test_non_integer_numeric_fields(self, bad):
        with pytest.raises(ParseError, match="non-integer"):
            parse_header(f"sp|P1AAA1|X_TEST A OS=Homo sapiens {bad}")

    def test_empty_file_yields_nothing(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text("")
        assert list(iter_records(str(p))) == []
