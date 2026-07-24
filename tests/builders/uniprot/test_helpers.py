"""Unit tests for UniProt builder field-level helpers."""

from bioparsers.builders.uniprot import helpers


class TestStripEvidence:

    def test_removes_tag_and_leading_space(self):
        assert helpers.strip_evidence("Heme {ECO:0000256|ARBA:ARBA00022617}") == "Heme"

    def test_removes_multiple_tags(self):
        s = "A {ECO:0000269|PubMed:1} and B {ECO:0000305}"
        assert helpers.strip_evidence(s) == "A and B"

    def test_no_tag_is_unchanged(self):
        assert helpers.strip_evidence("plain text.") == "plain text."

    def test_clean_text_normalizes_whitespace(self):
        assert helpers.clean_text("a   b\n c {ECO:0000305}") == "a b c"


class TestStripCitations:

    def test_single_pubmed_group(self):
        assert helpers.strip_citations("Binds zinc (PubMed:10707082).") == \
            "Binds zinc."

    def test_multiple_pubmed_group(self):
        s = "Active (PubMed:10707082, PubMed:11449278, PubMed:1936965)."
        assert helpers.strip_citations(s) == "Active."

    def test_pubmed_and_ref_mix(self):
        assert helpers.strip_citations("Cleaves (PubMed:30944155, Ref.2) it.") == \
            "Cleaves it."

    def test_and_joined_and_leading_comma_artifact(self):
        assert helpers.strip_citations("X (Ref.1 and PubMed:18301989) Y.") == "X Y."
        assert helpers.strip_citations("X (, PubMed:16895553) Y.") == "X Y."

    def test_descriptive_parenthetical_is_preserved(self):
        # An inline PubMed inside descriptive prose must NOT be excised.
        s = "Active (at pH 7.5, in PubMed:15629119) only."
        assert helpers.strip_citations(s) == s

    def test_clean_text_strips_citations_and_evidence(self):
        s = "Mixed {ECO:0000269|PubMed:1} and cite (PubMed:99)."
        assert helpers.clean_text(s) == "Mixed and cite."


class TestFullName:

    def _desc(self, **blocks):
        base = {"rec_name": None, "sub_name": None, "alt_names": [],
                "includes": [], "contains": [], "flags": []}
        base.update(blocks)
        return {"description": base}

    def test_prefers_rec_name(self):
        rec = self._desc(rec_name={"full": "Catalase", "short": [],
                                   "ec_numbers": [], "evidence": []})
        assert helpers.full_name(rec) == "Catalase"

    def test_falls_back_to_sub_name(self):
        rec = self._desc(sub_name={"full": "Uncharacterized protein",
                                   "short": [], "ec_numbers": [], "evidence": []})
        assert helpers.full_name(rec) == "Uncharacterized protein"

    def test_none_when_no_names(self):
        assert helpers.full_name(self._desc()) is None
        assert helpers.full_name({}) is None


class TestComments:

    rec = {"comments": [
        {"topic": "FUNCTION", "text": "Does a thing. {ECO:0000305}"},
        {"topic": "FUNCTION", "text": "Also another thing."},
        {"topic": "DOMAIN", "text": "Has two SH3 domains."},
    ]}

    def test_comment_texts_filters_by_topic(self):
        assert helpers.comment_texts(self.rec, "DOMAIN") == ["Has two SH3 domains."]

    def test_joined_comment_cleans_and_joins(self):
        assert helpers.joined_comment(self.rec, "FUNCTION") == \
            "Does a thing. Also another thing."

    def test_joined_comment_empty_when_absent(self):
        assert helpers.joined_comment(self.rec, "PATHWAY") == ""

    def test_joined_comment_concatenates_multiple_same_topic_blocks(self):
        # Several distinct DOMAIN comment blocks: each cleaned independently
        # (the ". {ECO}." double-period artifact repaired per block) and joined.
        rec = {"comments": [
            {"topic": "DOMAIN", "text": "First domain note. {ECO:0000305}."},
            {"topic": "DOMAIN", "text": "Second domain note."},
        ]}
        assert helpers.joined_comment(rec, "DOMAIN") == \
            "First domain note. Second domain note."

    def test_joined_comment_custom_separator(self):
        rec = {"comments": [
            {"topic": "FUNCTION", "text": "A."},
            {"topic": "FUNCTION", "text": "B."},
        ]}
        assert helpers.joined_comment(rec, "FUNCTION", sep=" | ") == "A. | B."


class TestKeywordsAndPfam:

    def test_keywords_strip_evidence(self):
        rec = {"keywords": ["Heme {ECO:0000256|ARBA:1}", "Transferase"]}
        assert helpers.keywords(rec) == ["Heme", "Transferase"]

    def test_pfam_ids_extracts_accession(self):
        rec = {"cross_references": {"Pfam": [
            "Pfam; PF00199; Catalase; 1.",
            "Pfam; PF06628; Catalase-rel; 1.",
        ]}}
        assert helpers.pfam_ids(rec) == ["PF00199", "PF06628"]

    def test_pfam_ids_empty_without_pfam(self):
        assert helpers.pfam_ids({"cross_references": {}}) == []
        assert helpers.pfam_ids({}) == []


class TestCofactorNames:
    """COFACTOR text is a structured record, not prose, so only the ``Name=``
    values belong in a caption. Legacy BioM3 renders ``COFACTOR: Mg(2+).``;
    taking the block verbatim would drag in ChEBI xrefs and evidence codes.
    """

    def test_extracts_name_dropping_xref_and_evidence(self):
        rec = {"comments": [{
            "topic": "COFACTOR",
            "text": "Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence=;",
        }]}
        assert helpers.cofactor_names(rec) == ["Mg(2+)"]

    def test_drops_note_prose(self):
        rec = {"comments": [{
            "topic": "COFACTOR",
            "text": ("Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Evidence=; "
                     "Note=Mg(2+) and Mn(2+) were both present in the kinase "
                     "buffer but Mg(2+) is likely the in vivo cofactor.;"),
        }]}
        assert helpers.cofactor_names(rec) == ["Mg(2+)"]

    def test_multiple_names_in_one_block(self):
        # Cu/Zn superoxide dismutase declares both metals in a single block.
        rec = {"comments": [{
            "topic": "COFACTOR",
            "text": ("Name=Cu cation; Xref=ChEBI:CHEBI:23378; Evidence=; "
                     "Note=Binds 1 copper ion per subunit.; "
                     "Name=Zn(2+); Xref=ChEBI:CHEBI:29105; Evidence=; "
                     "Note=Binds 1 zinc ion per subunit.;"),
        }]}
        assert helpers.cofactor_names(rec) == ["Cu cation", "Zn(2+)"]

    def test_multiple_blocks_are_collected_in_order(self):
        rec = {"comments": [
            {"topic": "COFACTOR", "text": "Name=Mn(2+); Xref=ChEBI:CHEBI:29035;"},
            {"topic": "FUNCTION", "text": "Unrelated."},
            {"topic": "COFACTOR", "text": "Name=Ca(2+); Xref=ChEBI:CHEBI:29108;"},
        ]}
        assert helpers.cofactor_names(rec) == ["Mn(2+)", "Ca(2+)"]

    def test_strips_evidence_tags_inside_the_name(self):
        rec = {"comments": [{
            "topic": "COFACTOR",
            "text": "Name=Zn(2+) {ECO:0000256|ARBA:ARBA00001947}; Xref=ChEBI:CHEBI:29105;",
        }]}
        assert helpers.cofactor_names(rec) == ["Zn(2+)"]

    def test_no_cofactor_comment_returns_empty(self):
        assert helpers.cofactor_names({"comments": [
            {"topic": "FUNCTION", "text": "Name=not a cofactor block;"},
        ]}) == []
        assert helpers.cofactor_names({}) == []

    def test_block_without_a_name_field_returns_empty(self):
        rec = {"comments": [{"topic": "COFACTOR", "text": "Note=Unspecified.;"}]}
        assert helpers.cofactor_names(rec) == []
