"""Unit tests for the UniProtKB REST fetcher.

No test here touches the network: :func:`_get_flat_text` is monkeypatched with
a stub that serves entries out of the real Swiss-Prot mini fixture, so the
batching, unresolved-accession, release-tracking and raw-tee logic is exercised
against genuine flat-file text without a live API.
"""

import os

import pytest

from bioparsers.fetch import uniprot_rest
from bioparsers.fetch.uniprot_rest import (
    FetchError,
    fetch_flat_text,
    is_valid_accession,
    iter_records,
)

DATDIR = os.path.join(os.path.dirname(__file__), "..", "_data")
SPROT = os.path.join(DATDIR, "uniprot_sprot_mini.dat")

RELEASE = {
    "release": "2026_02",
    "release_date": "10-June-2026",
    "api_deployment_date": "11-July-2026",
}


def _entry_blocks():
    """{primary accession: raw flat-file text} from the mini Swiss-Prot dump."""
    blocks, current = {}, []
    with open(SPROT) as handle:
        for line in handle:
            current.append(line)
            if line.rstrip("\n") == "//":
                accession = next(
                    ln[5:].split(";")[0].strip()
                    for ln in current if ln.startswith("AC   ")
                )
                blocks[accession] = "".join(current)
                current = []
    return blocks


@pytest.fixture(scope="module")
def blocks():
    return _entry_blocks()


@pytest.fixture
def fake_api(monkeypatch, blocks):
    """Stub the HTTP layer, recording each batch it was asked for.

    Returns the call log so tests can assert on batching. Unknown accessions
    are simply absent from the reply, mirroring the real endpoint.
    """
    calls = []

    def _stub(batch, context, retries, timeout):
        calls.append(list(batch))
        text = "".join(blocks[a] for a in batch if a in blocks)
        return text, dict(RELEASE)

    monkeypatch.setattr(uniprot_rest, "_get_flat_text", _stub)
    return calls


class TestAccessionSyntax:

    @pytest.mark.parametrize("accession", [
        "P00441", "Q6CPE2", "A0A072PZ83", "A0A385HVY4", "O95793",
    ])
    def test_valid(self, accession):
        assert is_valid_accession(accession)

    @pytest.mark.parametrize("accession", [
        "Q9NOTREAL", "not-an-id", "", "P0044", "p00441", "PF00080",
    ])
    def test_invalid(self, accession):
        assert not is_valid_accession(accession)


class TestRetrieval:

    def test_returns_records_for_known_accessions(self, fake_api, blocks):
        wanted = list(blocks)[:3]
        got = list(iter_records(wanted))
        assert [r["primary_accession"] for r in got] == wanted

    def test_records_are_fully_parsed(self, fake_api, blocks):
        accession = list(blocks)[0]
        (record,) = list(iter_records([accession]))
        assert record.record_type == "uniprot"
        assert record["sequence"]
        assert record["sequence_length"] == len(record["sequence"])

    def test_duplicates_collapsed_order_preserved(self, fake_api, blocks):
        a, b = list(blocks)[:2]
        got = list(iter_records([b, a, b]))
        assert [r["primary_accession"] for r in got] == [b, a]
        assert fake_api == [[b, a]]

    def test_empty_input_makes_no_request(self, fake_api):
        assert list(iter_records([])) == []
        assert fake_api == []


class TestBatching:

    def test_splits_into_batches(self, fake_api, blocks):
        wanted = list(blocks)[:10]
        list(iter_records(wanted, batch_size=4))
        assert fake_api == [wanted[0:4], wanted[4:8], wanted[8:10]]

    def test_batch_size_capped_at_endpoint_maximum(self, fake_api, blocks):
        wanted = list(blocks)[:5]
        list(iter_records(wanted, batch_size=10_000))
        assert len(fake_api) == 1

    def test_on_batch_reports_progress(self, fake_api, blocks):
        seen = []
        list(iter_records(list(blocks)[:5], batch_size=2,
                          on_batch=lambda i, n: seen.append((i, n))))
        assert seen == [(1, 3), (2, 3), (3, 3)]


class TestUnresolvedAccessions:
    """The two ways an accession yields no record mean different things and
    are reported separately — neither may be silently dropped."""

    def test_missing_reported_and_others_still_returned(self, fake_api, blocks):
        known = list(blocks)[:2]
        missing = []
        got = list(iter_records([known[0], "Q00000", known[1]],
                                on_missing=missing.extend))
        assert [r["primary_accession"] for r in got] == known
        assert missing == ["Q00000"]

    def test_invalid_excluded_from_request_not_sent_to_api(self, fake_api, blocks):
        known = list(blocks)[0]
        invalid = []
        got = list(iter_records([known, "Q9NOTREAL", "not-an-id"],
                                on_invalid=invalid.extend))
        assert [r["primary_accession"] for r in got] == [known]
        assert invalid == ["Q9NOTREAL", "not-an-id"]
        # The malformed ids never reach the endpoint: one bad accession would
        # otherwise make the API reject the whole batch.
        assert fake_api == [[known]]

    def test_callbacks_not_invoked_when_everything_resolves(self, fake_api, blocks):
        missing, invalid = [], []
        list(iter_records(list(blocks)[:3], on_missing=missing.append,
                          on_invalid=invalid.append))
        assert missing == [] and invalid == []

    def test_missing_reported_per_batch(self, fake_api, blocks):
        known = list(blocks)[:2]
        batches = []
        list(iter_records([known[0], "Q00000", known[1], "Q00001"],
                          batch_size=2, on_missing=batches.append))
        assert batches == [["Q00000"], ["Q00001"]]

    def test_secondary_accession_counts_as_found(self, fake_api, monkeypatch, blocks):
        """A record is matched on every AC it carries, not just the primary, so
        a secondary accession merged into another entry is not called missing.
        """
        primary = next(a for a in blocks if blocks[a].count("AC   ") and
                       len(_secondaries(blocks[a])) > 0)
        secondary = _secondaries(blocks[primary])[0]

        def _stub(batch, context, retries, timeout):
            return blocks[primary], dict(RELEASE)

        monkeypatch.setattr(uniprot_rest, "_get_flat_text", _stub)
        missing = []
        list(iter_records([secondary], on_missing=missing.extend))
        assert missing == []

    def test_all_invalid_makes_no_request(self, fake_api):
        invalid = []
        assert list(iter_records(["bad", "worse"], on_invalid=invalid.extend)) == []
        assert fake_api == []
        assert invalid == ["bad", "worse"]


def _secondaries(text):
    """Non-primary accessions declared on an entry's AC lines."""
    accessions = []
    for line in text.splitlines():
        if line.startswith("AC   "):
            accessions += [t.strip() for t in line[5:].split(";") if t.strip()]
    return accessions[1:]


class TestReleaseTracking:

    def test_release_reported_once(self, fake_api, blocks):
        releases = []
        list(iter_records(list(blocks)[:6], batch_size=2,
                          on_release=releases.append))
        assert releases == [RELEASE]

    def test_distinct_releases_each_reported(self, monkeypatch, blocks):
        """A long run can straddle a release switchover; both are recorded
        rather than the run being assumed homogeneous."""
        seen = {"n": 0}
        second = dict(RELEASE, release="2026_03")

        def _stub(batch, context, retries, timeout):
            seen["n"] += 1
            release = RELEASE if seen["n"] == 1 else second
            return "".join(blocks[a] for a in batch if a in blocks), dict(release)

        monkeypatch.setattr(uniprot_rest, "_get_flat_text", _stub)
        releases = []
        list(iter_records(list(blocks)[:4], batch_size=2,
                          on_release=releases.append))
        assert releases == [RELEASE, second]


class TestRawSink:

    def test_tees_flat_text_that_reparses_to_the_same_records(
            self, fake_api, blocks, tmp_path):
        from bioparsers.parsers.uniprot_dat import iter_records as parse_dat

        wanted = list(blocks)[:6]
        chunks = []
        fetched = list(iter_records(wanted, batch_size=2, raw_sink=chunks.append))

        path = tmp_path / "saved.dat"
        path.write_text("".join(chunks))
        reparsed = list(parse_dat(str(path)))

        assert [r.as_dict() for r in reparsed] == [r.as_dict() for r in fetched]

    def test_not_called_when_absent(self, fake_api, blocks):
        # Nothing to assert beyond "no exception"; the sink is optional.
        assert list(iter_records(list(blocks)[:2])) != []


class TestFetchFlatText:

    def test_returns_raw_text(self, fake_api, blocks):
        accession = list(blocks)[0]
        assert fetch_flat_text([accession]) == blocks[accession]

    def test_rejects_oversized_batch(self):
        with pytest.raises(ValueError, match="exceeds the endpoint maximum"):
            fetch_flat_text(["P00441"] * (uniprot_rest.MAX_BATCH + 1))


class TestErrorPropagation:

    def test_fetch_error_propagates(self, monkeypatch, blocks):
        def _boom(batch, context, retries, timeout):
            raise FetchError("HTTP 400 Bad Request")

        monkeypatch.setattr(uniprot_rest, "_get_flat_text", _boom)
        with pytest.raises(FetchError, match="400"):
            list(iter_records(["P00441"]))

    def test_malformed_payload_raises_parse_error(self, monkeypatch):
        def _garbage(batch, context, retries, timeout):
            return "this is not a flat file\n", dict(RELEASE)

        monkeypatch.setattr(uniprot_rest, "_get_flat_text", _garbage)
        with pytest.raises(uniprot_rest.ParseError, match="expected an ID line"):
            list(iter_records(["P00441"]))


class TestBackoff:

    def test_uses_retry_after_when_sent(self):
        assert uniprot_rest._backoff(1, "5") == 5.0

    def test_retry_after_capped(self):
        assert uniprot_rest._backoff(1, "9999") == 60.0

    def test_exponential_without_retry_after(self):
        assert [uniprot_rest._backoff(n, None) for n in (1, 2, 3, 4)] == \
            [1.0, 2.0, 4.0, 8.0]

    def test_unparseable_retry_after_falls_back(self):
        assert uniprot_rest._backoff(2, "in a while") == 2.0
