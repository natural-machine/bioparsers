"""Tests for build-manifest generation."""

import json

from bioparsers.builders import Builder, Stage, generate_manifest, write_manifest


class _Demo(Builder):
    """Demo output: {x}."""
    name = "demo_v0"

    def build(self, records):
        yield from records


class TestGenerateManifest:

    def test_core_fields(self):
        m = generate_manifest(_Demo(), description="a note",
                              output="out.jsonl", record_count=42)
        assert m["builder"] == {"name": "demo_v0", "description": "Demo output: {x}."}
        assert m["description"] == "a note"
        assert m["output"] == "out.jsonl"
        assert m["record_count"] == 42
        assert m["bioparsers_version"]                      # non-empty string
        assert m["bioparsers_build"].startswith(m["bioparsers_version"])
        assert set(m["git"]) == {"hash", "dirty", "branch", "remote"}
        assert "python" in m["environment"] and "platform" in m["environment"]
        assert m["timestamp"] and m["command"]

    def test_defaults_are_none(self):
        m = generate_manifest(_Demo())
        assert m["description"] is None
        assert m["output"] is None
        assert m["record_count"] is None

    def test_extra_merged_at_top_level(self):
        m = generate_manifest(_Demo(), extra={"pfam_ids": ["PF1", "PF2"], "join": True})
        assert m["pfam_ids"] == ["PF1", "PF2"]
        assert m["join"] is True

    def test_write_manifest_round_trips(self, tmp_path):
        p = tmp_path / "m.json"
        ret = write_manifest(_Demo(), str(p), description="x", record_count=3)
        assert ret == str(p)
        data = json.loads(p.read_text())
        assert data["builder"]["name"] == "demo_v0"
        assert data["record_count"] == 3
        assert data["description"] == "x"


class TestStageProducer:
    """``Stage`` lets a non-Builder step (a fetch, a download) get the same
    provenance sidecar — ``generate_manifest`` only ever reads name/description.
    """

    def test_stage_carries_name_and_description(self):
        stage = Stage("uniprot_rest_fetch", "Entries retrieved by accession.")
        m = generate_manifest(stage)
        assert m["builder"] == {
            "name": "uniprot_rest_fetch",
            "description": "Entries retrieved by accession.",
        }

    def test_stage_gets_the_standard_provenance_keys(self):
        m = generate_manifest(Stage("s", "d"))
        for key in ("bioparsers_version", "git", "timestamp", "command",
                    "environment"):
            assert key in m

    def test_extra_keys_merge_at_top_level(self):
        m = generate_manifest(
            Stage("s", "d"), record_count=2,
            extra={"unresolved": {"invalid_format": ["bad"],
                                  "not_returned": ["Q00000"]}},
        )
        assert m["record_count"] == 2
        assert m["unresolved"]["not_returned"] == ["Q00000"]

    def test_write_manifest_round_trips(self, tmp_path):
        path = tmp_path / "out.jsonl.manifest.json"
        write_manifest(Stage("s", "d"), str(path), extra={"counts": {"n": 1}})
        assert json.loads(path.read_text())["counts"] == {"n": 1}
