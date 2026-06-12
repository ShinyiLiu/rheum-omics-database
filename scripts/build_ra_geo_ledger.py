#!/usr/bin/env python3
"""Build a GEO ledger for human rheumatoid arthritis expression datasets."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import http.client
import html
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "ra_geo_transcriptome_datasets.csv"
DEFAULT_STATE = ROOT / "data" / "ra_geo_backfill_state.json"
DEFAULT_LATEST = ROOT / "data" / "ra_geo_latest_monthly_update.csv"
DEFAULT_MONTHLY_DIR = ROOT / "data" / "monthly"
NCBI_LAST_REQUEST_AT = 0.0
NCBI_MIN_INTERVAL_SECONDS = 0.45

FIELDS = [
    "accession",
    "title",
    "summary",
    "organism",
    "technology",
    "geo_study_type",
    "gpl_accessions",
    "gpl_titles",
    "sample_count",
    "sample_count_confidence",
    "score",
    "score_reason",
    "ra_relevance",
    "tissue",
    "case_control_hint",
    "treatment_hint",
    "clinical_info_hint",
    "processed_data_hint",
    "publication_date",
    "last_update_date",
    "url",
    "search_date",
    "notes",
]

RA_TERMS = [
    "rheumatoid arthritis",
    "rheumatoid-arthritis",
    "rheumatoid synovitis",
    "类风湿性关节炎",
    "类风湿关节炎",
]

EXCLUDE_TERMS = [
    "single-cell",
    "single cell",
    "single-nucleus",
    "single nucleus",
    "scrna",
    "snrna",
    "spatial transcriptomics",
    "spatial transcriptome",
    "visium",
    "merfish",
    "xenium",
    "cosmx",
    "methylation",
    "atac-seq",
    "atac seq",
    "chip-seq",
    "chip seq",
    "mirna",
    "microrna",
    "genotyping",
    "snp array",
]

NON_RA_AMBIGUOUS = [
    "retinoic acid",
    "all-trans retinoic acid",
    "atr a",
]

TISSUE_TERMS = [
    ("synovium", ["synovium", "synovial tissue", "synovial membrane"]),
    ("synovial fibroblast", ["fibroblast-like synoviocyte", "fibroblast like synoviocyte", "fls", "synovial fibroblast"]),
    ("blood", ["whole blood", "blood"]),
    ("PBMC", ["pbmc", "peripheral blood mononuclear"]),
    ("macrophage", ["macrophage"]),
    ("monocyte", ["monocyte"]),
    ("T cell", ["t cell", "t-cell", "cd4", "cd8"]),
    ("B cell", ["b cell", "b-cell"]),
]

CLINICAL_TERMS = ["das28", "crp", "esr", "acr", "remission", "activity", "response", "clinical"]
TREATMENT_TERMS = ["treatment", "therapy", "drug", "response", "tnf", "anti-tnf", "methotrexate", "mtx", "jak", "inhibitor", "tocilizumab", "abatacept", "rituximab"]
PROCESSED_TERMS = ["processed", "normalized", "counts", "count matrix", "matrix", "fpkm", "tpm", "raw counts", "supplementary", "series matrix"]
MIXED_DISEASE_TERMS = ["osteoarthritis", "psoriasis", "psoriatic", "lupus", "systemic sclerosis", "sjogren", "spondyloarthritis", "healthy and diseased"]

GSE_RE = re.compile(r"\bGSE\d+\b", re.IGNORECASE)
GSM_RE = re.compile(r"\bGSM\d+\b", re.IGNORECASE)
GPL_RE = re.compile(r"\bGPL\d+\b", re.IGNORECASE)


@dataclass
class GeoRecord:
    accession: str
    title: str
    summary: str
    organism: str
    geo_study_type: str
    publication_date: str
    last_update_date: str
    total_samples: int | None
    url: str
    soft_text: str = ""


@dataclass
class SampleCounts:
    sample_count: str
    confidence: str
    total: int | None
    ra: int | None
    control: int | None
    other: int | None
    notes: list[str]


def ncbi_wait() -> None:
    global NCBI_LAST_REQUEST_AT
    elapsed = time.monotonic() - NCBI_LAST_REQUEST_AT
    if elapsed < NCBI_MIN_INTERVAL_SECONDS:
        time.sleep(NCBI_MIN_INTERVAL_SECONDS - elapsed)
    NCBI_LAST_REQUEST_AT = time.monotonic()


def add_ncbi_identity(params: dict[str, str]) -> dict[str, str]:
    enriched = dict(params)
    enriched["tool"] = os.environ.get("NCBI_TOOL", "ra_geo_transcriptome_ledger")
    if os.environ.get("NCBI_EMAIL"):
        enriched["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        enriched["api_key"] = os.environ["NCBI_API_KEY"]
    return enriched


def http_text(url: str, timeout: int = 60, retries: int = 3) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "codex-ra-geo-ledger/1.0"})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            retryable = exc.code in {429, 500, 502, 503, 504}
            if not retryable or attempt >= retries:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = float(retry_after) if retry_after and retry_after.isdigit() else min(60.0, 2.0 * (2**attempt))
            time.sleep(delay)
        except (http.client.IncompleteRead, TimeoutError, OSError) as exc:
            if attempt >= retries:
                raise
            delay = min(60.0, 2.0 * (2**attempt))
            time.sleep(delay)
    raise RuntimeError("unreachable retry loop")


def http_json(url: str, timeout: int = 60, retries: int = 3) -> Any:
    return json.loads(http_text(url, timeout=timeout, retries=retries))


def parse_iso_date(value: str) -> dt.date:
    try:
        return dt.date.fromisoformat(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid date {value!r}; expected YYYY-MM-DD") from exc


def ncbi_date(value: dt.date) -> str:
    return value.strftime("%Y/%m/%d")


def window_key(start_date: dt.date | None, end_date: dt.date | None) -> str:
    if start_date and end_date:
        return f"{start_date.isoformat()}..{end_date.isoformat()}"
    return ""


def quarter_end(value: dt.date) -> dt.date:
    quarter_last_month = ((value.month - 1) // 3 + 1) * 3
    if quarter_last_month == 12:
        next_month = dt.date(value.year + 1, 1, 1)
    else:
        next_month = dt.date(value.year, quarter_last_month + 1, 1)
    return next_month - dt.timedelta(days=1)


def next_day(value: dt.date) -> dt.date:
    return value + dt.timedelta(days=1)


def quarterly_windows(start_date: dt.date, end_date: dt.date) -> list[tuple[dt.date, dt.date]]:
    windows: list[tuple[dt.date, dt.date]] = []
    cursor = start_date
    while cursor <= end_date:
        current_end = min(quarter_end(cursor), end_date)
        windows.append((cursor, current_end))
        cursor = next_day(current_end)
    return windows


def monthly_archive_path(monthly_dir: Path, today: dt.date | None = None) -> Path:
    stamp = (today or dt.date.today()).strftime("%Y-%m")
    return monthly_dir / f"ra_geo_monthly_update_{stamp}.csv"


def ncbi_search(
    term: str,
    retmax: int,
    since_days: int | None,
    start_date: dt.date | None = None,
    end_date: dt.date | None = None,
) -> list[str]:
    params = {
        "db": "gds",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "pub date",
    }
    if start_date and end_date:
        params["datetype"] = "pdat"
        params["mindate"] = ncbi_date(start_date)
        params["maxdate"] = ncbi_date(end_date)
    elif since_days and since_days > 0:
        params["datetype"] = "pdat"
        params["reldate"] = str(since_days)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(add_ncbi_identity(params))
    ncbi_wait()
    data = http_json(url)
    return data.get("esearchresult", {}).get("idlist", [])


def ncbi_summary(ids: list[str]) -> dict[str, Any]:
    if not ids:
        return {}
    params = add_ncbi_identity({"db": "gds", "id": ",".join(ids), "retmode": "json"})
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(params)
    ncbi_wait()
    return http_json(url).get("result", {})


def clean_text(value: Any) -> str:
    return " ".join(html.unescape(str(value or "")).split())


def text_has_any(text: str, terms: list[str]) -> bool:
    folded = text.casefold()
    return any(term.casefold() in folded for term in terms)


def strict_ra_match(text: str) -> tuple[bool, str]:
    folded = text.casefold()
    for term in RA_TERMS:
        if term.casefold() in folded:
            return True, f"matched {term}"
    return False, ""


def is_human_record(text: str) -> bool:
    folded = text.casefold()
    return "homo sapiens" in folded or re.search(r"\bhuman(s)?\b", folded) is not None


def infer_study_type(text: str) -> tuple[str, str]:
    folded = text.casefold()
    if "expression profiling by array" in folded:
        return "expression array", "Expression profiling by array"
    if "expression profiling by high throughput sequencing" in folded:
        return "bulk RNA-seq", "Expression profiling by high throughput sequencing"
    if "rna-seq" in folded or "rna seq" in folded or "transcriptome sequencing" in folded:
        return "bulk RNA-seq", "Expression profiling by high throughput sequencing"
    return "", ""


def positive_int(value: Any) -> int | None:
    if isinstance(value, int) and value >= 0:
        return value
    text = str(value or "").replace(",", "").strip()
    if re.fullmatch(r"\d+", text):
        return int(text)
    return None


def extract_total_samples(item: dict[str, Any], text: str) -> int | None:
    for key in ("n_samples", "nsamples", "sample_count", "samples", "gdsType"):
        count = positive_int(item.get(key))
        if count is not None:
            return count
    gsm = {match.upper() for match in GSM_RE.findall(text)}
    return len(gsm) if gsm else None


def fetch_soft(accession: str) -> str:
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urllib.parse.urlencode({"acc": accession, "form": "text"})
    return http_text(url, timeout=90, retries=3)


def extract_series_sample_ids(soft_text: str) -> list[str]:
    ids: list[str] = []
    for raw_line in soft_text.splitlines():
        if raw_line.startswith("!Series_sample_id") and "=" in raw_line:
            _, value = raw_line.split("=", 1)
            sample_id = clean_text(value).upper()
            if GSM_RE.fullmatch(sample_id):
                ids.append(sample_id)
    return list(dict.fromkeys(ids))


def fetch_sample_softs(sample_ids: list[str], max_samples: int, workers: int) -> tuple[str, list[str]]:
    warnings: list[str] = []
    if max_samples <= 0:
        if sample_ids:
            warnings.append(f"sample fetch capped at 0 of {len(sample_ids)} GSM records")
        return "", warnings

    capped_ids = sample_ids[:max_samples]
    chunks_by_id: dict[str, str] = {}

    def fetch_one(sample_id: str) -> tuple[str, str]:
        return sample_id, fetch_soft(sample_id)

    worker_count = max(1, workers)
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {executor.submit(fetch_one, sample_id): sample_id for sample_id in capped_ids}
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                fetched_id, text = future.result()
                chunks_by_id[fetched_id] = text
            except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
                warnings.append(f"{sample_id}: failed to fetch sample SOFT: {exc}")

    chunks: list[str] = []
    for sample_id in capped_ids:
        try:
            chunks.append(chunks_by_id[sample_id])
        except KeyError:
            pass
    if len(sample_ids) > max_samples:
        warnings.append(f"sample fetch capped at {max_samples} of {len(sample_ids)} GSM records")
    return "\n".join(chunks), warnings


def parse_soft_blocks(soft_text: str, prefix: str) -> list[dict[str, list[str]]]:
    blocks: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    for raw_line in soft_text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith(prefix):
            if current:
                blocks.append(current)
            _, value = line.split("=", 1)
            current = {"accession": [clean_text(value)]}
            continue
        if current is None or not line.startswith("!"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        current.setdefault(key.strip(), []).append(clean_text(value))
    if current:
        blocks.append(current)
    return blocks


def parse_soft_series_fields(soft_text: str) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    for raw_line in soft_text.splitlines():
        if raw_line.startswith("^SAMPLE"):
            break
        if raw_line.startswith("!Series_") and "=" in raw_line:
            key, value = raw_line.split("=", 1)
            result.setdefault(key.strip(), []).append(clean_text(value))
    return {key: "; ".join(values) for key, values in result.items()}


def sample_text(sample: dict[str, list[str]]) -> str:
    useful_keys = [
        "!Sample_title",
        "!Sample_source_name_ch1",
        "!Sample_characteristics_ch1",
        "!Sample_description",
        "!Sample_treatment_protocol_ch1",
        "!Sample_extract_protocol_ch1",
    ]
    pieces: list[str] = []
    for key in useful_keys:
        pieces.extend(sample.get(key, []))
    return "\n".join(pieces)


def classify_sample_group(text: str) -> str:
    folded = text.casefold()
    if any(term in folded for term in NON_RA_AMBIGUOUS):
        return "unclear"
    if "rheumatoid arthritis" in folded or "rheumatoid-arthritis" in folded or "类风湿" in folded:
        return "ra"
    if re.search(r"\bra\b", folded) and not text_has_any(folded, NON_RA_AMBIGUOUS):
        disease_context = any(token in folded for token in ["patient", "synovium", "synovial", "arthritis", "pbmc", "blood"])
        if disease_context:
            return "ra"
    control_patterns = [
        "healthy control",
        "healthy controls",
        "normal control",
        "normal controls",
        "control donor",
        "control donors",
        "healthy donor",
        "healthy donors",
        "normal donor",
        "normal donors",
        "non-ra",
        "non rheumatoid",
        "osteoarthritis control",
        "control subject",
        "control subjects",
    ]
    if any(pattern in folded for pattern in control_patterns):
        return "control"
    if re.search(r"(?<![a-z0-9])hd(?![a-z0-9])", folded):
        return "control"
    if re.search(r"\b(control|healthy|normal)\b", folded):
        return "control"
    return "unclear"


def build_sample_count(samples: list[dict[str, list[str]]], total_hint: int | None, gse_text: str) -> SampleCounts:
    notes: list[str] = []
    if samples:
        ra = control = unclear = 0
        for sample in samples:
            group = classify_sample_group(sample_text(sample))
            if group == "ra":
                ra += 1
            elif group == "control":
                control += 1
            else:
                unclear += 1
        total = len(samples)
        parts = [f"RA={ra}", f"control={control}"]
        if unclear:
            parts.append(f"other/unclear={unclear}")
        parts.append(f"total={total}")
        confidence = "high" if ra > 0 and control > 0 else ("medium" if ra > 0 else "low")
        if confidence != "high":
            notes.append("manual review needed for sample groups")
        return SampleCounts("; ".join(parts), confidence, total, ra, control, unclear, notes)

    total = total_hint
    if total is None:
        total = len({match.upper() for match in GSM_RE.findall(gse_text)}) or None
    if total is not None:
        return SampleCounts(f"RA=unknown; control=unknown; total={total}", "low", total, None, None, None, ["manual review needed for sample groups"])
    return SampleCounts("RA=unknown; control=unknown; total=unknown", "low", None, None, None, None, ["manual review needed for sample groups"])


def extract_gpl_info(soft_text: str, summary_text: str) -> tuple[str, str]:
    accessions: list[str] = []
    titles: dict[str, str] = {}
    for platform in parse_soft_blocks(soft_text, "^PLATFORM"):
        accession = platform.get("accession", [""])[0].upper()
        if accession:
            accessions.append(accession)
            title = "; ".join(platform.get("!Platform_title", []))
            if title:
                titles[accession] = title
    if not accessions:
        accessions = sorted({match.upper() for match in GPL_RE.findall(summary_text + "\n" + soft_text)})
    return "; ".join(dict.fromkeys(accessions)), "; ".join(titles.get(acc, "") for acc in dict.fromkeys(accessions) if titles.get(acc))


def infer_tissue(text: str) -> str:
    folded = text.casefold()
    found = [label for label, terms in TISSUE_TERMS if any(term in folded for term in terms)]
    return "; ".join(dict.fromkeys(found))


def yes_no_hint(text: str, terms: list[str]) -> str:
    return "yes" if text_has_any(text, terms) else "no"


def case_control_hint(counts: SampleCounts, text: str) -> str:
    if counts.control and counts.control > 0:
        return "yes"
    return "yes" if re.search(r"\b(control|healthy|normal)\b", text.casefold()) else "no"


def processed_hint(text: str) -> str:
    return "yes" if text_has_any(text, PROCESSED_TERMS) else "unknown"


def score_row(
    *,
    text: str,
    organism: str,
    technology: str,
    counts: SampleCounts,
    tissue: str,
    clinical_hint: str,
    processed: str,
) -> tuple[int, str]:
    score = 0
    reasons: list[str] = []
    if is_human_record(organism + "\n" + text):
        score += 20
        reasons.append("human")
    if counts.ra is not None and counts.ra > 0:
        score += 20
        reasons.append("RA samples clear")
    if counts.control is not None and counts.control > 0:
        score += 15
        reasons.append("controls clear")
    if counts.ra is not None and counts.control is not None and (counts.ra > 0 or counts.control > 0):
        score += 15
        reasons.append("RA/control counts shown")
    if counts.total is not None and counts.total >= 20:
        score += 10
        reasons.append("total samples >=20")
    if processed == "yes":
        score += 10
        reasons.append("processed data likely available")
    if clinical_hint == "yes":
        score += 5
        reasons.append("clinical metadata mentioned")
    if tissue:
        score += 5
        reasons.append("tissue/source clear")

    if counts.ra is None or counts.control is None or (counts.ra == 0 and counts.control == 0):
        score -= 20
        reasons.append("RA/control groups unclear")
    if counts.ra is not None and counts.ra == 0:
        score -= 20
        reasons.append("no RA samples parsed")
    if text_has_any(text, MIXED_DISEASE_TERMS) and not (counts.ra is not None and counts.ra > 0):
        score -= 15
        reasons.append("mixed disease cohort unclear")
    if not technology:
        score -= 15
        reasons.append("platform/data type unclear")
    return max(0, min(100, score)), "; ".join(reasons)


def row_for_record(record: GeoRecord) -> dict[str, str] | None:
    soft = record.soft_text
    series_fields = parse_soft_series_fields(soft)
    samples = parse_soft_blocks(soft, "^SAMPLE") if soft else []
    merged_text = "\n".join([
        record.title,
        record.summary,
        record.organism,
        record.geo_study_type,
        soft[:200000],
        "; ".join(series_fields.values()),
    ])

    ra_ok, ra_reason = strict_ra_match(merged_text)
    if not ra_ok:
        return None
    if not is_human_record(merged_text):
        return None
    if text_has_any(merged_text, EXCLUDE_TERMS):
        return None
    technology, geo_study_type = infer_study_type(merged_text)
    if not technology:
        return None

    counts = build_sample_count(samples, record.total_samples, merged_text)
    gpl_accessions, gpl_titles = extract_gpl_info(soft, merged_text)
    tissue = infer_tissue(merged_text)
    clinical = yes_no_hint(merged_text, CLINICAL_TERMS)
    treatment = yes_no_hint(merged_text, TREATMENT_TERMS)
    processed = processed_hint(merged_text)
    score, score_reason = score_row(
        text=merged_text,
        organism=record.organism,
        technology=technology,
        counts=counts,
        tissue=tissue,
        clinical_hint=clinical,
        processed=processed,
    )
    notes = ["auto-collected from GEO metadata"]
    notes.extend(counts.notes)
    if not gpl_accessions:
        notes.append("manual review needed for GPL platform")

    today = dt.date.today().isoformat()
    return {
        "accession": record.accession,
        "title": record.title,
        "summary": record.summary,
        "organism": record.organism,
        "technology": technology,
        "geo_study_type": geo_study_type or record.geo_study_type,
        "gpl_accessions": gpl_accessions,
        "gpl_titles": gpl_titles,
        "sample_count": counts.sample_count,
        "sample_count_confidence": counts.confidence,
        "score": str(score),
        "score_reason": score_reason,
        "ra_relevance": ra_reason,
        "tissue": tissue,
        "case_control_hint": case_control_hint(counts, merged_text),
        "treatment_hint": treatment,
        "clinical_info_hint": clinical,
        "processed_data_hint": processed,
        "publication_date": record.publication_date,
        "last_update_date": record.last_update_date,
        "url": record.url,
        "search_date": today,
        "notes": "; ".join(dict.fromkeys(note for note in notes if note)),
    }


def build_search_term() -> str:
    ra = '"rheumatoid arthritis"'
    human = '("Homo sapiens"[Organism] OR human)'
    expression = '("Expression profiling by array" OR "Expression profiling by high throughput sequencing" OR "RNA-seq" OR "transcriptome sequencing")'
    return f"{ra} AND {human} AND {expression}"


def records_from_summaries(summaries: dict[str, Any]) -> list[GeoRecord]:
    records: list[GeoRecord] = []
    for uid in summaries.get("uids", []):
        item = summaries.get(uid, {})
        title = clean_text(item.get("title") or item.get("summary") or "")
        summary = clean_text(item.get("summary") or "")
        item_payload = json.dumps(item, ensure_ascii=False)
        accession = clean_text(item.get("accession") or item.get("gse") or "")
        if not accession.upper().startswith("GSE"):
            match = GSE_RE.search(item_payload)
            accession = match.group(0).upper() if match else accession
        if not accession:
            accession = str(uid)
        if not accession.upper().startswith("GSE"):
            continue
        organism = clean_text(item.get("taxon") or item.get("organism") or "")
        type_text = clean_text(item.get("gdstype") or item.get("entryType") or item.get("gdsType") or "")
        combined = "\n".join([title, summary, organism, type_text, item_payload])
        if not strict_ra_match(combined)[0]:
            continue
        records.append(
            GeoRecord(
                accession=accession.upper(),
                title=title or accession.upper(),
                summary=summary,
                organism=organism,
                geo_study_type=type_text,
                publication_date=clean_text(item.get("PDAT") or item.get("pdat") or item.get("pubdate") or ""),
                last_update_date=clean_text(item.get("updatedate") or item.get("suppdate") or item.get("lastupdate") or ""),
                total_samples=extract_total_samples(item, combined),
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={urllib.parse.quote(accession.upper())}",
            )
        )
    dedup: dict[str, GeoRecord] = {}
    for record in records:
        dedup.setdefault(record.accession, record)
    return list(dedup.values())


def collect_rows(
    retmax: int,
    since_days: int | None,
    max_samples_per_series: int,
    sample_workers: int,
    start_date: dt.date | None = None,
    end_date: dt.date | None = None,
) -> tuple[list[dict[str, str]], list[str], int, int]:
    warnings: list[str] = []
    ids = ncbi_search(build_search_term(), retmax, since_days, start_date, end_date)
    summaries = ncbi_summary(ids)
    records = records_from_summaries(summaries)
    rows: list[dict[str, str]] = []
    for record in records:
        try:
            series_soft = fetch_soft(record.accession)
            sample_ids = extract_series_sample_ids(series_soft)
            sample_soft, sample_warnings = fetch_sample_softs(sample_ids, max_samples_per_series, sample_workers) if sample_ids else ("", [])
            record.soft_text = series_soft + "\n" + sample_soft
            warnings.extend(f"{record.accession}: {warning}" for warning in sample_warnings)
        except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
            warnings.append(f"{record.accession}: failed to fetch SOFT text: {exc}")
        row = row_for_record(record)
        if row:
            rows.append(row)
        time.sleep(0.2)
    rows.sort(key=lambda row: (-int(row.get("score") or 0), row.get("accession", "")))
    return rows, warnings, len(ids), len(records)


def load_state(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"version": 1, "windows": {}}
    with path.open("r", encoding="utf-8") as handle:
        state = json.load(handle)
    if not isinstance(state, dict):
        return {"version": 1, "windows": {}}
    state.setdefault("version", 1)
    state.setdefault("windows", {})
    return state


def write_state(path: Path, state: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(state, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def state_is_complete(state: dict[str, Any], key: str) -> bool:
    window = state.get("windows", {}).get(key, {})
    return window.get("status") in {"completed", "completed_with_warnings"}


def update_state(
    state: dict[str, Any],
    key: str,
    *,
    start_date: dt.date,
    end_date: dt.date,
    queried_ids: int,
    candidate_gse_count: int,
    rows_written_for_window: int,
    final_csv_rows: int,
    warnings: list[str],
    started_at: str,
    finished_at: str,
) -> None:
    state.setdefault("windows", {})[key] = {
        "status": "completed_with_warnings" if warnings else "completed",
        "start_date": start_date.isoformat(),
        "end_date": end_date.isoformat(),
        "started_at": started_at,
        "finished_at": finished_at,
        "queried_geo_ids": queried_ids,
        "candidate_gse_count": candidate_gse_count,
        "rows_written_for_window": rows_written_for_window,
        "final_csv_rows": final_csv_rows,
        "warnings": warnings,
    }


def update_recent_state(
    state: dict[str, Any],
    *,
    since_days: int,
    queried_ids: int,
    candidate_gse_count: int,
    rows_written_for_window: int,
    final_csv_rows: int,
    warnings: list[str],
    started_at: str,
    finished_at: str,
) -> None:
    state.setdefault("recent_updates", []).append(
        {
            "status": "completed_with_warnings" if warnings else "completed",
            "since_days": since_days,
            "started_at": started_at,
            "finished_at": finished_at,
            "queried_geo_ids": queried_ids,
            "candidate_gse_count": candidate_gse_count,
            "rows_written_for_window": rows_written_for_window,
            "final_csv_rows": final_csv_rows,
            "warnings": warnings,
        }
    )


def read_existing(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            accession = (row.get("accession") or "").upper()
            if accession:
                rows[accession] = {field: row.get(field, "") for field in FIELDS}
        return rows


def parsed_sample_total(sample_count: str) -> int | None:
    match = re.search(r"(?:^|;\s*)total=(\d+)", sample_count or "")
    return int(match.group(1)) if match else None


def preserve_more_complete_sample_count(old: dict[str, str], new: dict[str, str]) -> None:
    old_total = parsed_sample_total(old.get("sample_count", ""))
    new_total = parsed_sample_total(new.get("sample_count", ""))
    if old_total is None or new_total is None or new_total >= old_total:
        return
    for field in ("sample_count", "sample_count_confidence", "score", "score_reason"):
        if old.get(field):
            new[field] = old[field]
    note = f"preserved previous sample_count with total={old_total}; latest run parsed total={new_total}"
    if note not in new.get("notes", ""):
        new["notes"] = "; ".join(part for part in [new.get("notes", ""), note] if part)


def merge_rows(existing: dict[str, dict[str, str]], new_rows: list[dict[str, str]], include_existing: bool) -> list[dict[str, str]]:
    merged = dict(existing) if include_existing else {}
    for row in new_rows:
        accession = row["accession"].upper()
        old = merged.get(accession, {})
        combined = {field: row.get(field, "") for field in FIELDS}
        if old:
            preserve_more_complete_sample_count(old, combined)
        if old.get("notes") and old["notes"] not in combined.get("notes", ""):
            combined["notes"] = "; ".join(part for part in [combined.get("notes", ""), old["notes"]] if part)
        merged[accession] = combined
    rows = list(merged.values())
    rows.sort(key=lambda row: (-int(row.get("score") or 0), row.get("accession", "")))
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_monthly_update(latest_file: Path, monthly_dir: Path, rows: list[dict[str, str]]) -> Path:
    monthly_rows = sorted(rows, key=lambda row: (-int(row.get("score") or 0), row.get("accession", "")))
    archive_file = monthly_archive_path(monthly_dir)
    write_csv(latest_file, monthly_rows)
    write_csv(archive_file, monthly_rows)
    return archive_file


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=500, help="Maximum GEO records to request from ESearch.")
    parser.add_argument("--since-days", type=int, default=None, help="Restrict GEO search to records published within this many days.")
    parser.add_argument("--start-date", type=parse_iso_date, default=None, help="Start publication date for a fixed GEO search window, YYYY-MM-DD.")
    parser.add_argument("--end-date", type=parse_iso_date, default=None, help="End publication date for a fixed GEO search window, YYYY-MM-DD.")
    parser.add_argument("--backfill-quarterly-local", action="store_true", help="Run local quarterly backfill from start-date to end-date.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output CSV path.")
    parser.add_argument("--include-existing", action="store_true", help="Keep existing rows not found in the current search.")
    parser.add_argument("--state-file", type=Path, default=None, help="Optional backfill state JSON path for fixed date windows.")
    parser.add_argument("--force-window", action="store_true", help="Re-run a fixed date window even if state says it is complete.")
    parser.add_argument("--latest-file", type=Path, default=DEFAULT_LATEST, help="Latest monthly update CSV path.")
    parser.add_argument("--monthly-dir", type=Path, default=DEFAULT_MONTHLY_DIR, help="Directory for archived monthly update CSVs.")
    parser.add_argument("--write-monthly-update", action="store_true", help="Write latest and archived monthly update CSVs from current search rows.")
    parser.add_argument("--max-samples-per-series", type=int, default=250, help="Maximum GSM records to fetch per GSE for RA/control counts.")
    parser.add_argument("--sample-workers", type=int, default=4, help="Concurrent GSM fetch workers per GSE.")
    parser.add_argument("--ncbi-email", default="", help="Optional NCBI email identity.")
    parser.add_argument("--ncbi-api-key", default="", help="Optional NCBI API key.")
    parser.add_argument("--validate-csv", action="store_true", help="Validate output CSV after writing.")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if bool(args.start_date) != bool(args.end_date):
        raise SystemExit("--start-date and --end-date must be provided together")
    if args.start_date and args.end_date and args.start_date > args.end_date:
        raise SystemExit("--start-date must be on or before --end-date")
    if args.since_days and args.start_date:
        raise SystemExit("--since-days cannot be combined with --start-date/--end-date")
    if args.force_window and not (args.start_date and args.end_date and args.state_file):
        raise SystemExit("--force-window requires --start-date, --end-date, and --state-file")
    if args.backfill_quarterly_local and not (args.start_date and args.end_date and args.state_file):
        raise SystemExit("--backfill-quarterly-local requires --start-date, --end-date, and --state-file")
    if args.backfill_quarterly_local and args.write_monthly_update:
        raise SystemExit("--backfill-quarterly-local cannot be combined with --write-monthly-update")


def validate_csv(path: Path) -> None:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    accessions = [row.get("accession", "") for row in rows]
    duplicates = sorted({accession for accession in accessions if accessions.count(accession) > 1})
    if duplicates:
        raise RuntimeError(f"duplicate accessions in CSV: {', '.join(duplicates[:10])}")
    scores = [int(row.get("score") or 0) for row in rows]
    if scores != sorted(scores, reverse=True):
        raise RuntimeError("CSV scores are not sorted in descending order")
    missing_counts = [row.get("accession", "") for row in rows if "RA=" not in row.get("sample_count", "") or "control=" not in row.get("sample_count", "")]
    if missing_counts:
        raise RuntimeError(f"sample_count missing RA/control format: {', '.join(missing_counts[:10])}")


def run_one_window(
    *,
    args: argparse.Namespace,
    state: dict[str, Any] | None,
    start_date: dt.date | None,
    end_date: dt.date | None,
) -> tuple[list[dict[str, str]], list[str], int, int, list[dict[str, str]]]:
    started_at = dt.datetime.now(dt.timezone.utc).isoformat()
    rows, warnings, queried_ids, candidate_count = collect_rows(
        args.retmax,
        args.since_days if start_date is None else None,
        args.max_samples_per_series,
        args.sample_workers,
        start_date,
        end_date,
    )
    existing = read_existing(args.out)
    final_rows = merge_rows(existing, rows, True if args.backfill_quarterly_local else args.include_existing)
    write_csv(args.out, final_rows)
    if args.write_monthly_update:
        archive_file = write_monthly_update(args.latest_file, args.monthly_dir, rows)
        print(f"Wrote monthly update rows to {args.latest_file} and {archive_file}")
    if args.validate_csv:
        validate_csv(args.out)

    key = window_key(start_date, end_date)
    if state is not None and key and start_date and end_date:
        update_state(
            state,
            key,
            start_date=start_date,
            end_date=end_date,
            queried_ids=queried_ids,
            candidate_gse_count=candidate_count,
            rows_written_for_window=len(rows),
            final_csv_rows=len(final_rows),
            warnings=warnings,
            started_at=started_at,
            finished_at=dt.datetime.now(dt.timezone.utc).isoformat(),
        )
        write_state(args.state_file, state)
    elif state is not None and args.since_days:
        update_recent_state(
            state,
            since_days=args.since_days,
            queried_ids=queried_ids,
            candidate_gse_count=candidate_count,
            rows_written_for_window=len(rows),
            final_csv_rows=len(final_rows),
            warnings=warnings,
            started_at=started_at,
            finished_at=dt.datetime.now(dt.timezone.utc).isoformat(),
        )
        write_state(args.state_file, state)

    return rows, warnings, queried_ids, candidate_count, final_rows


def run_quarterly_backfill(args: argparse.Namespace, state: dict[str, Any]) -> int:
    total_new_rows = 0
    for start_date, end_date in quarterly_windows(args.start_date, args.end_date):
        key = window_key(start_date, end_date)
        if state_is_complete(state, key) and not args.force_window:
            print(f"Skipping completed window {key}; use --force-window to re-run it.")
            continue
        print(f"Running window {key}")
        rows, warnings, queried_ids, candidate_count, final_rows = run_one_window(
            args=args,
            state=state,
            start_date=start_date,
            end_date=end_date,
        )
        total_new_rows += len(rows)
        print(f"Window {key}: {len(rows)} rows from {candidate_count} candidate GSE records ({queried_ids} GEO ids)")
        if warnings:
            print("Warnings:")
            for warning in warnings:
                print(f"- {warning}")
        print(f"Total CSV rows now: {len(final_rows)}")
    print(f"Quarterly backfill completed; wrote {total_new_rows} window rows across requested range.")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    validate_args(args)
    if args.ncbi_email:
        os.environ["NCBI_EMAIL"] = args.ncbi_email
    if args.ncbi_api_key:
        os.environ["NCBI_API_KEY"] = args.ncbi_api_key

    state: dict[str, Any] | None = None
    key = window_key(args.start_date, args.end_date)
    if args.state_file:
        state = load_state(args.state_file)
    if args.backfill_quarterly_local:
        return run_quarterly_backfill(args, state or {"version": 1, "windows": {}})
    if state is not None and key:
        if state_is_complete(state, key) and not args.force_window:
            print(f"Skipping completed window {key}; use --force-window to re-run it.")
            return 0

    rows, warnings, queried_ids, candidate_count, final_rows = run_one_window(
        args=args,
        state=state,
        start_date=args.start_date,
        end_date=args.end_date,
    )
    print(f"Wrote {len(final_rows)} rows to {args.out}")
    if key:
        print(f"Window {key}: {len(rows)} rows from {candidate_count} candidate GSE records ({queried_ids} GEO ids)")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
