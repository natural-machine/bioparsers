"""Retrieve UniProtKB entries by accession from the UniProt REST API.

Serves the case a local scan handles badly: a known, modest list of accessions
scattered through a very large distribution file. The
``/uniprotkb/accessions`` endpoint takes up to 100 accessions per request and
can return ``format=txt`` — the **same flat-file format** as
``uniprot_sprot.dat`` / ``uniprot_trembl.dat``, line codes and ``//``
terminators intact. So nothing here parses anything: the payload goes straight
to :func:`bioparsers.parsers.uniprot_dat.iter_records_from_lines` and the
records are ordinary ``UniProtRecord``s, identical in shape to a local parse.

Contract: ``iter_records(accessions) -> Iterator[UniProtRecord]``. Swiss-Prot
and TrEMBL accessions may be mixed freely — the endpoint covers both sections,
and each record's ``status`` says which it came from.

    from bioparsers.fetch.uniprot_rest import iter_records
    from bioparsers.parsers import dump_jsonl

    with open("out.jsonl", "w") as f:
        dump_jsonl(iter_records(["P00441", "A0A072PZ83"]), f)

Because the payload *is* the ``.dat`` format, it is worth keeping. Pass
*raw_sink* to tee each batch's text to disk as it arrives; the concatenation
is a valid flat file that reparses offline to the same records, so a later
parser fix does not mean refetching, and the exact bytes the API served on a
given day stay pinned even as the release moves on::

    with open("entries.dat", "w") as raw, open("out.jsonl", "w") as out:
        dump_jsonl(iter_records(accessions, raw_sink=raw.write), out)

Two properties worth knowing before relying on the result:

**Unknown accessions are skipped, not raised.** The endpoint returns the
entries it recognizes and says nothing about the rest, so a batch of 100
containing one deleted accession yields 99 records. Nothing is silently lost
in a way the caller cannot see: pass *on_missing* to observe each such
accession as it is detected, and note that the count returned is authoritative
(records are never fabricated). Accessions are matched against every AC on the
entry, not just the primary, so a secondary accession that has since been
merged into another entry still counts as found.

**The API serves the current release.** A local ``.dat`` snapshot is frozen;
this is not. Accessions valid in an older snapshot may have been deleted or
merged, and annotation on the rest may have drifted. When exact agreement with
a specific release matters more than speed, scan the local file instead.

Fail-loud (raises ``FetchError``)
---------------------------------
- a non-retryable HTTP error (4xx other than 429)
- a persistent 429 / 5xx / connection failure, after *retries* attempts
- a payload that is not valid UniProt flat-file text (as ``ParseError``,
  which ``FetchError`` does not wrap — it propagates from the parser)
"""

from __future__ import annotations

import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Callable, Iterable, Iterator, Sequence

from bioparsers.parsers.base import ParseError
from bioparsers.parsers.uniprot_dat import UniProtRecord, iter_records_from_lines

#: The batch endpoint. Documented at https://www.uniprot.org/help/api_retrieve_entries
ENDPOINT = "https://rest.uniprot.org/uniprotkb/accessions"

#: Accessions per request. 100 is the server-side maximum for this endpoint;
#: a larger value is rejected with HTTP 400.
MAX_BATCH = 100

#: HTTP statuses worth another attempt: rate limiting and transient server or
#: gateway faults. Everything else (400, 404, ...) is a caller error.
_RETRY_STATUS = {429, 500, 502, 503, 504}

#: The official UniProtKB accession syntax, from
#: https://www.uniprot.org/help/accession_numbers. Checked client-side because
#: the endpoint rejects a whole request when any one accession is malformed —
#: a single bad identifier would otherwise cost the other 99 in its batch.
_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def is_valid_accession(accession: str) -> bool:
    """Return whether *accession* matches the UniProtKB accession syntax.

    A syntax check only — a well-formed accession may still be unknown to the
    API (deleted, or never assigned), which surfaces as a missing record
    rather than an error.
    """
    return bool(_ACCESSION_RE.match(accession))

_USER_AGENT = "bioparsers (https://github.com/; UniProtKB batch retrieval)"

#: Checked in order when the interpreter's default SSL context has no CA
#: certificates loaded — the usual state of a bare conda environment, which
#: ships no ``certifi``. Verification is never disabled; if none of these
#: exists the original verification error is allowed to surface.
_CA_BUNDLES = (
    "/etc/ssl/certs/ca-certificates.crt",   # Debian / Ubuntu / Alpine
    "/etc/pki/tls/certs/ca-bundle.crt",     # RHEL / Fedora
    "/etc/ssl/ca-bundle.pem",               # openSUSE
    "/etc/ssl/cert.pem",                    # macOS / BSD
)


class FetchError(Exception):
    """A request could not be completed — a non-retryable HTTP status, or
    retries exhausted against a rate-limited/failing server.

    Deliberately not a ``ParseError``: it means the bytes never arrived, not
    that they arrived malformed. A malformed payload still raises
    ``ParseError`` from the flat-file parser.
    """


# ===========================================================================
# Public API
# ===========================================================================

def iter_records(
    accessions: Iterable[str],
    *,
    batch_size: int = MAX_BATCH,
    on_missing: Callable[[list[str]], None] | None = None,
    on_invalid: Callable[[list[str]], None] | None = None,
    on_batch: Callable[[int, int], None] | None = None,
    on_release: Callable[[dict], None] | None = None,
    raw_sink: Callable[[str], None] | None = None,
    retries: int = 5,
    timeout: float = 120.0,
) -> Iterator[UniProtRecord]:
    """Yield a :class:`UniProtRecord` for each retrievable accession in
    *accessions*, in batches of *batch_size* (capped at :data:`MAX_BATCH`).

    Duplicates are collapsed and original order is preserved. Two kinds of
    accession yield no record, reported through separate callbacks because
    they mean different things — neither is silently dropped:

    - *on_invalid* is called once, before any request, with accessions that
      fail :func:`is_valid_accession`. These are excluded from the batches;
      sending one would make the API reject its whole batch of 100.
    - *on_missing* is called once per batch with that batch's well-formed
      accessions the API did not return — deleted, or never assigned.

    Neither callback is invoked with an empty list. *on_batch* receives
    ``(batch_index, batch_count)`` before each request, for progress
    reporting. *on_release* receives each distinct ``{release, release_date,
    api_deployment_date}`` the responses report — normally once, but more if a
    long run straddles a UniProt release switchover.

    *raw_sink* is handed each batch's flat-file text verbatim, before it is
    parsed. Concatenating those chunks gives a valid UniProt flat file — the
    archival artifact the API response actually was, reparseable offline with
    :func:`bioparsers.parsers.uniprot_dat.iter_records` and byte-comparable to
    the corresponding slice of a ``.dat`` distribution. Teeing rather than
    returning keeps the whole payload from being held in memory.

    Raises :class:`FetchError` on a non-retryable HTTP status or once
    *retries* attempts at a retryable one are exhausted.
    """
    unique = list(dict.fromkeys(a for a in accessions if a))
    wanted = [a for a in unique if is_valid_accession(a)]
    if on_invalid is not None:
        invalid = [a for a in unique if not is_valid_accession(a)]
        if invalid:
            on_invalid(invalid)
    if not wanted:
        return

    size = max(1, min(batch_size, MAX_BATCH))
    batches = [wanted[i : i + size] for i in range(0, len(wanted), size)]
    context = _ssl_context()

    seen_releases: list[dict] = []

    for index, batch in enumerate(batches, 1):
        if on_batch is not None:
            on_batch(index, len(batches))
        text, release = _get_flat_text(batch, context, retries, timeout)
        # A multi-batch fetch can straddle a UniProt release switchover; report
        # each distinct release so the caller records both rather than assuming
        # the run was homogeneous.
        if on_release is not None and release not in seen_releases:
            seen_releases.append(release)
            on_release(release)
        if raw_sink is not None:
            raw_sink(text)
        returned: set[str] = set()
        for record in iter_records_from_lines(
            text.splitlines(keepends=True),
            source=f"{ENDPOINT} (batch {index}/{len(batches)})",
        ):
            returned.update(record["accessions"])
            yield record
        if on_missing is not None:
            absent = [a for a in batch if a not in returned]
            if absent:
                on_missing(absent)


def fetch_flat_text(
    accessions: Sequence[str],
    *,
    retries: int = 5,
    timeout: float = 120.0,
) -> str:
    """Return the raw ``format=txt`` flat-file payload for up to
    :data:`MAX_BATCH` *accessions* — the escape hatch for callers that want
    the bytes rather than parsed records (caching them to disk, say).

    For a whole dataset prefer ``iter_records(..., raw_sink=...)``, which
    batches automatically and tees the same text to disk.

    Raises ``ValueError`` if given more than :data:`MAX_BATCH` accessions.
    """
    if len(accessions) > MAX_BATCH:
        raise ValueError(
            f"{len(accessions)} accessions exceeds the endpoint maximum of "
            f"{MAX_BATCH}; use iter_records() to batch automatically"
        )
    text, _release = _get_flat_text(
        list(accessions), _ssl_context(), retries, timeout
    )
    return text


# ===========================================================================
# Implementation details
# ===========================================================================

def _ssl_context() -> ssl.SSLContext:
    """Return a verifying SSL context, falling back to a system CA bundle when
    the default context loaded no certificates (a conda env without
    ``certifi``). Verification stays on in every branch.
    """
    context = ssl.create_default_context()
    if context.cert_store_stats().get("x509_ca", 0) > 0:
        return context
    for path in _CA_BUNDLES:
        try:
            context.load_verify_locations(cafile=path)
        except OSError:
            continue
        if context.cert_store_stats().get("x509_ca", 0) > 0:
            return context
    return context


def _get_flat_text(
    batch: list[str],
    context: ssl.SSLContext,
    retries: int,
    timeout: float,
) -> tuple[str, dict]:
    """GET one batch as flat-file text, retrying transient failures with
    exponential backoff (honouring ``Retry-After`` when the server sends it).

    Returns ``(text, release)`` where *release* is the served UniProt release
    read off the response headers — see :func:`_release_info`.
    """
    query = urllib.parse.urlencode(
        {"accessions": ",".join(batch), "format": "txt"}
    )
    request = urllib.request.Request(
        f"{ENDPOINT}?{query}", headers={"User-Agent": _USER_AGENT}
    )

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(
                request, timeout=timeout, context=context
            ) as response:
                text = response.read().decode("utf-8")
                return text, _release_info(response.headers)
        except urllib.error.HTTPError as exc:
            if exc.code not in _RETRY_STATUS:
                detail = _error_detail(exc)
                raise FetchError(
                    f"{ENDPOINT}: HTTP {exc.code} {exc.reason} for "
                    f"{len(batch)} accessions ({batch[0]}...){detail}"
                ) from exc
            if attempt == retries:
                raise FetchError(
                    f"{ENDPOINT}: HTTP {exc.code} {exc.reason} after "
                    f"{retries} attempts"
                ) from exc
            time.sleep(_backoff(attempt, exc.headers.get("Retry-After")))
        except (urllib.error.URLError, TimeoutError, ssl.SSLError) as exc:
            if attempt == retries:
                raise FetchError(
                    f"{ENDPOINT}: request failed after {retries} attempts: {exc}"
                ) from exc
            time.sleep(_backoff(attempt, None))

    raise FetchError(f"{ENDPOINT}: exhausted {retries} attempts")  # unreachable


def _release_info(headers) -> dict:
    """Extract the served UniProt release from a response's headers.

    The API tracks the current release rather than a frozen snapshot, so which
    release answered is the one provenance fact a fetch cannot reconstruct
    afterwards. UniProt reports it as ``X-UniProt-Release`` (e.g. ``2026_02``)
    with ``X-UniProt-Release-Date``. Values are ``None`` if absent — recorded
    as unknown rather than guessed.
    """
    return {
        "release": headers.get("X-UniProt-Release"),
        "release_date": headers.get("X-UniProt-Release-Date"),
        "api_deployment_date": headers.get("X-API-Deployment-Date"),
    }


def _backoff(attempt: int, retry_after: str | None) -> float:
    """Seconds to wait before *attempt*+1: the server's ``Retry-After`` when it
    sent a usable one, else 1, 2, 4, 8, ... capped at 60."""
    if retry_after:
        try:
            return min(float(retry_after), 60.0)
        except ValueError:
            pass
    return min(2.0 ** (attempt - 1), 60.0)


def _error_detail(exc: urllib.error.HTTPError) -> str:
    """Best-effort body snippet from an error response — UniProt explains 400s
    there (e.g. a malformed accession), which the status alone does not."""
    try:
        body = exc.read().decode("utf-8", "replace").strip()
    except Exception:
        return ""
    return f": {body[:300]}" if body else ""


#: ``ParseError`` is re-exported so a caller can catch both failure modes —
#: bytes that never arrived, and bytes that arrived malformed — from one import.
__all__ = [
    "ENDPOINT", "MAX_BATCH", "FetchError", "ParseError",
    "iter_records", "fetch_flat_text", "is_valid_accession",
]
