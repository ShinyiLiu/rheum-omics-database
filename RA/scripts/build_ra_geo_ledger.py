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
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "ra_geo_transcriptome_datasets.csv"
DEFAULT_STATE = ROOT / "data" / "ra_geo_backfill_state.json"
DEFAULT_LATEST = ROOT / "data" / "ra_geo_latest_monthly_update.csv"
DEFAULT_MONTHLY_DIR = ROOT / "data" / "monthly"
DEFAULT_IMPACT_FACTOR_CACHE = ROOT / "data" / "online_impact_factor_cache.json"
DEFAULT_LOCAL_JCR_DIR = ROOT.parent / "IF reference"
NCBI_LAST_REQUEST_AT = 0.0
NCBI_MIN_INTERVAL_SECONDS = 0.45

FIELDS = [
    "accession",
    "title",
    "organism",
    "technology",
    "geo_study_type",
    "gpl_accessions",
    "sample_count",
    "sample_count_confidence",
    "score",
    "score_reason",
    "ra_relevance",
    "tissue",
    "case_control_hint",
    "treatment_hint",
    "clinical_info_hint",
    "source_journals",
    "source_journal_impact_factor",
    "source_journal_impact_factor_year",
    "gse_reuse_total_count",
    "publication_date",
    "last_update_date",
    "url",
    "search_date",
    "manual_review",
    "notes",
]

MANUAL_REVIEW_VALUES = {"NULL", "YES", "NO"}


def normalize_manual_review(value: str | None) -> str:
    normalized = (value or "NULL").strip().upper()
    if normalized not in MANUAL_REVIEW_VALUES:
        raise ValueError(f"manual_review must be one of {sorted(MANUAL_REVIEW_VALUES)}, got {value!r}")
    return normalized

MESH_HEADING = "Arthritis, Rheumatoid"
MESH_ENTRY_TERMS = [
    "Rheumatoid Arthritis",
]
CLINICAL_SEARCH_TERMS = [
    "rheumatoid arthritis",
    "rheumatoid-arthritis",
    "rheumatoid synovitis",
    "RA",
    "???????",
    "??????",
]
RA_TERMS = [MESH_HEADING, *MESH_ENTRY_TERMS, *CLINICAL_SEARCH_TERMS]

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


@dataclass
class GseReuseCounts:
    total: int


@dataclass(frozen=True)
class JournalImpactFactor:
    impact_factor: str
    year: str
    source: str
    source_url: str = ""


@dataclass
class ImpactFactorAudit:
    static_matches: dict[str, JournalImpactFactor]
    local_jcr_matches: dict[str, JournalImpactFactor]
    cache_matches: dict[str, JournalImpactFactor]
    online_matches: dict[str, JournalImpactFactor]
    missing: dict[str, list[str]]
    skipped: dict[str, list[str]]
    warnings: list[str]


PREPRINT_JOURNAL_PATTERNS = [
    "biorxiv",
    "medrxiv",
    "preprint server",
    "research square",
]


# Publicly verifiable JIF values collected from journal/publisher pages or
# pages that explicitly cite Journal Citation Reports/Web of Science.
ONLINE_IMPACT_FACTORS = {
    "arthritisrheumatologyhobokennj": JournalImpactFactor("10.9", "2024", "Wiley Online Library journal metrics"),
    "arthritisrheumatology": JournalImpactFactor("10.9", "2024", "Wiley Online Library journal metrics"),
    "haematologica": JournalImpactFactor("7.9", "2024", "Haematologica journal metrics"),
    "thejournalofbiologicalchemistry": JournalImpactFactor("3.9", "2024", "ASBMB/JBC journal metrics"),
    "journalofbiologicalchemistry": JournalImpactFactor("3.9", "2024", "ASBMB/JBC journal metrics"),
    "immunologicresearch": JournalImpactFactor("3.1", "2024", "Springer Nature journal metrics"),
    "thejournalofphysiology": JournalImpactFactor("4.4", "2024", "Wiley Online Library journal metrics"),
    "journalofphysiology": JournalImpactFactor("4.4", "2024", "Wiley Online Library journal metrics"),
    "thejournalofinvestigativedermatology": JournalImpactFactor("5.7", "2024", "Journal of Investigative Dermatology journal metrics"),
    "journalofinvestigativedermatology": JournalImpactFactor("5.7", "2024", "Journal of Investigative Dermatology journal metrics"),
    "proceedingsofthenationalacademyofsciencesoftheunitedstatesofamerica": JournalImpactFactor("9.1", "2024", "PNAS journal metrics"),
    "proceedingsofthenationalacademyofsciences": JournalImpactFactor("9.1", "2024", "PNAS journal metrics"),
    "naturecommunications": JournalImpactFactor("15.7", "2024", "Nature Communications journal metrics"),
    "scientificreports": JournalImpactFactor("3.9", "2024", "Scientific Reports journal metrics"),
    "frontiersinimmunology": JournalImpactFactor("5.9", "", "Frontiers in Immunology journal page"),
    "arthritisresearchtherapy": JournalImpactFactor("4.6", "2025", "LetPub/科研通 journal metrics"),
    "plosone": JournalImpactFactor("2.6", "2024-2025", "LetPub journal metrics"),
    "sciencetranslationalmedicine": JournalImpactFactor("14.7", "2024-2025", "LetPub journal metrics"),
    "annalsoftherheumaticdiseases": JournalImpactFactor("20.6", "2024", "RheumNow/JCR-labelled report"),
    "immunity": JournalImpactFactor("26.3", "2024-2025", "LetPub journal metrics"),
    "journalofautoimmunity": JournalImpactFactor("7", "2025", "科研通 journal metrics"),
    "cellularmolecularimmunology": JournalImpactFactor("19.8", "2024-2025", "LetPub journal metrics"),
    "jciinsight": JournalImpactFactor("6.1", "2024", "JCI Insight journal page"),
    "communicationsbiology": JournalImpactFactor("5.1", "2025", "科研通 journal metrics"),
    "scienceadvances": JournalImpactFactor("12.5", "2024-2025", "LetPub journal metrics"),
    "internationaljournalofrheumaticdiseases": JournalImpactFactor("2", "2024-2025", "LetPub journal metrics"),
    "rheumatologyoxfordengland": JournalImpactFactor("4.4", "2024-2025", "LetPub journal metrics"),
    "clinicalandexperimentalimmunology": JournalImpactFactor("3.8", "2024", "Oxford Academic/Ovid JCR-labelled metrics"),
    "cellularandmolecularlifesciencescmls": JournalImpactFactor("6.2", "2024-2025", "LetPub journal metrics"),
    "cellularandmolecularlifesciences": JournalImpactFactor("6.2", "2024-2025", "LetPub journal metrics"),
    "iscience": JournalImpactFactor("4.1", "2024-2025", "LetPub journal metrics"),
    "cells": JournalImpactFactor("5.2", "2024-2025", "LetPub journal metrics"),
    "experimentalandtherapeuticmedicine": JournalImpactFactor("2.3", "2024-2025", "LetPub journal metrics"),
    "natureimmunology": JournalImpactFactor("27.6", "2024", "Nature Immunology journal metrics"),
    "journalofleukocytebiology": JournalImpactFactor("3.1", "2025", "Journal Metrics JCR-labelled listing"),
    "advancedscienceweinheimbadenwurttemberggermany": JournalImpactFactor("14.3", "2024-2025", "LetPub journal metrics"),
    "advancedscience": JournalImpactFactor("14.3", "2024-2025", "LetPub journal metrics"),
    "internationaljournalofmolecularsciences": JournalImpactFactor("4.9", "2024", "MDPI journal metrics"),
    "elife": JournalImpactFactor("0", "2024-2025", "LetPub journal metrics"),
    "scientificdata": JournalImpactFactor("6.9", "2024", "Nature Portfolio journal metrics"),
    "europeanjournalofimmunology": JournalImpactFactor("4.4", "2024-2025", "LetPub/online journal metrics"),
    "thejournalofmoleculardiagnosticsjmd": JournalImpactFactor("3.4", "2024", "Association for Molecular Pathology journal metrics"),
    "journalofmoleculardiagnosticsjmd": JournalImpactFactor("3.4", "2024", "Association for Molecular Pathology journal metrics"),
    "scienceimmunology": JournalImpactFactor("16.4", "2024-2025", "LetPub journal metrics"),
    "moleculartherapythejournaloftheamericansocietyofgenetherapy": JournalImpactFactor("12", "2024-2025", "LetPub journal metrics"),
    "moleculartherapy": JournalImpactFactor("12", "2024-2025", "LetPub journal metrics"),
    "actabiomaterialia": JournalImpactFactor("9.6", "2024-2025", "LetPub journal metrics"),
    "cellreportsmedicine": JournalImpactFactor("10.6", "2024-2025", "LetPub/科研通 journal metrics"),
    "peerj": JournalImpactFactor("2.4", "2024-2025", "LetPub/PeerJ journal metrics"),
    "journalofnanobiotechnology": JournalImpactFactor("12.6", "2025", "科研通 journal metrics"),
    "molecularpain": JournalImpactFactor("3.3", "2024", "PMC/Web of Science-based review"),
    "clinicalandtranslationalmedicine": JournalImpactFactor("6.8", "2024-2025", "LetPub/科研通 journal metrics"),
    "journaloftranslationalautoimmunity": JournalImpactFactor("3.6", "2025", "Journal Metrics JCR-labelled listing"),
    "clinicalandexperimentalrheumatology": JournalImpactFactor("3.3", "2024-2025", "LetPub journal metrics"),
    "biomedicines": JournalImpactFactor("3.9", "2024-2025", "LetPub journal metrics"),
    "chinesemedicine": JournalImpactFactor("5.7", "2024-2025", "LetPub journal metrics"),
    "humangenetherapy": JournalImpactFactor("4.5", "2024-2025", "LetPub journal metrics"),
    "journaloftranslationalmedicine": JournalImpactFactor("7.5", "2025", "科研通/Springer Nature journal metrics"),
    "journalofimmunologybaltimoremd1950": JournalImpactFactor("3.4", "2024", "Oxford Academic journal metrics"),
    "ebiomedicine": JournalImpactFactor("10.8", "2024-2025", "LetPub/科研通 journal metrics"),
    "emboreports": JournalImpactFactor("6.2", "2024-2025", "LetPub journal metrics"),
    "jhepreportsinnovationinhepatology": JournalImpactFactor("7.5", "2025", "EASL/JCR-labelled report"),
    "frontiersinpharmacology": JournalImpactFactor("4.8", "2024-2025", "LetPub/Frontiers journal metrics"),
    "biomedicalreports": JournalImpactFactor("1.9", "2025", "Journal Metrics JCR-labelled listing"),
    "genesandimmunity": JournalImpactFactor("4.5", "2024-2025", "LetPub journal metrics"),
    "molecularbiologyreports": JournalImpactFactor("2.8", "2024", "Springer Nature journal metrics"),
}


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
    merged: dict[str, Any] = {}
    for index in range(0, len(ids), 100):
        batch = ids[index : index + 100]
        params = add_ncbi_identity({"db": "gds", "id": ",".join(batch), "retmode": "json"})
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(params)
        ncbi_wait()
        merged.update(http_json(url).get("result", {}))
    return merged


def ncbi_esearch(db: str, term: str, retmax: int = 0) -> dict[str, Any]:
    params = {
        "db": db,
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
    }
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(add_ncbi_identity(params))
    ncbi_wait()
    return http_json(url).get("esearchresult", {})


def ncbi_pubmed_summary(ids: list[str]) -> dict[str, Any]:
    if not ids:
        return {}
    params = add_ncbi_identity({"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
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
    if any(term in folded for term in NON_RA_AMBIGUOUS):
        return False, ""
    for term in RA_TERMS:
        term_folded = term.casefold()
        if term_folded == "ra":
            if re.search(r"ra", folded):
                disease_context = any(token in folded for token in ["patient", "synovium", "synovial", "arthritis", "pbmc", "blood"])
                if disease_context:
                    return True, "matched RA in disease context"
        elif term_folded in folded:
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


def extract_source_pubmed_ids(soft_text: str) -> list[str]:
    ids: list[str] = []
    for raw_line in soft_text.splitlines():
        if raw_line.startswith("^SAMPLE"):
            break
        if raw_line.startswith("!Series_pubmed_id") and "=" in raw_line:
            _, value = raw_line.split("=", 1)
            ids.extend(re.findall(r"\d+", value))
    return list(dict.fromkeys(ids))


def normalize_journal_name(value: str) -> str:
    folded = clean_text(value).casefold()
    return re.sub(r"[^a-z0-9]+", "", folded)


def xlsx_column_index(cell_reference: str) -> int:
    letters = "".join(char for char in cell_reference if char.isalpha())
    index = 0
    for char in letters.upper():
        index = index * 26 + (ord(char) - ord("A") + 1)
    return max(0, index - 1)


def xlsx_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for item in root:
        parts = [node.text or "" for node in item.iter() if node.tag.endswith("}t") or node.tag == "t"]
        strings.append("".join(parts))
    return strings


def xlsx_cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return clean_text("".join(node.text or "" for node in cell.iter() if node.tag.endswith("}t") or node.tag == "t"))
    value_node = next((child for child in cell if child.tag.endswith("}v") or child.tag == "v"), None)
    if value_node is None or value_node.text is None:
        return ""
    value = value_node.text
    if cell_type == "s":
        try:
            return clean_text(shared_strings[int(value)])
        except (ValueError, IndexError):
            return ""
    return clean_text(value)


def iter_xlsx_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with zipfile.ZipFile(path) as archive:
        sheet_names = sorted(name for name in archive.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", name))
        if not sheet_names:
            return rows
        shared_strings = xlsx_shared_strings(archive)
        root = ET.fromstring(archive.read(sheet_names[0]))
        for row_node in root.iter():
            if not (row_node.tag.endswith("}row") or row_node.tag == "row"):
                continue
            row_values: list[str] = []
            for cell in row_node:
                if not (cell.tag.endswith("}c") or cell.tag == "c"):
                    continue
                column = xlsx_column_index(cell.attrib.get("r", ""))
                while len(row_values) <= column:
                    row_values.append("")
                row_values[column] = xlsx_cell_text(cell, shared_strings)
            rows.append(row_values)
    return rows


def load_local_jcr_impact_factors(reference_dir: Path = DEFAULT_LOCAL_JCR_DIR) -> dict[str, JournalImpactFactor]:
    if not reference_dir.exists():
        return {}
    for path in sorted(reference_dir.glob("*.xlsx")):
        if path.name.startswith("~$"):
            continue
        rows = iter_xlsx_rows(path)
        header_index = -1
        journal_index = -1
        impact_index = -1
        year = ""
        for index, row in enumerate(rows):
            folded = [clean_text(value).casefold() for value in row]
            for column, value in enumerate(folded):
                if "journal name" in value:
                    journal_index = column
                if "impact factor" in value or "影响因子" in value:
                    impact_index = column
                    year_match = re.search(r"(20\d{2})", value)
                    if year_match:
                        year = year_match.group(1)
            if journal_index >= 0 and impact_index >= 0:
                header_index = index
                break
        if header_index < 0:
            continue
        mapping: dict[str, JournalImpactFactor] = {}
        source = f"local JCR XLSX: {path.name}"
        for row in rows[header_index + 1 :]:
            if max(journal_index, impact_index) >= len(row):
                continue
            journal = clean_text(row[journal_index])
            impact_factor = clean_text(row[impact_index])
            if not journal or not impact_factor or impact_factor.upper() in {"N/A", "NA"}:
                continue
            key = normalize_journal_name(journal)
            if key:
                mapping[key] = JournalImpactFactor(impact_factor=impact_factor, year=year, source=source)
        if mapping:
            return mapping
    return {}


def is_preprint_or_non_jcr_journal(value: str) -> bool:
    folded = clean_text(value).casefold()
    return any(pattern in folded for pattern in PREPRINT_JOURNAL_PATTERNS)


def source_impact_factor(
    journals: list[str],
    cache: dict[str, JournalImpactFactor] | None = None,
    local_jcr: dict[str, JournalImpactFactor] | None = None,
) -> JournalImpactFactor | None:
    matches: list[JournalImpactFactor] = []
    for journal in journals:
        key = normalize_journal_name(journal)
        if local_jcr and key in local_jcr:
            matches.append(local_jcr[key])
        elif key in ONLINE_IMPACT_FACTORS:
            matches.append(ONLINE_IMPACT_FACTORS[key])
        elif cache and key in cache:
            matches.append(cache[key])
    if not matches:
        return None
    return max(matches, key=lambda impact: float(impact.impact_factor or 0))


def load_impact_factor_cache(path: Path) -> dict[str, JournalImpactFactor]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    raw_journals = payload.get("journals", payload) if isinstance(payload, dict) else {}
    cache: dict[str, JournalImpactFactor] = {}
    if not isinstance(raw_journals, dict):
        return cache
    for key, item in raw_journals.items():
        if not isinstance(item, dict):
            continue
        impact_factor = clean_text(item.get("impact_factor") or "")
        year = clean_text(item.get("year") or "")
        if impact_factor and year:
            cache[str(key)] = JournalImpactFactor(
                impact_factor=impact_factor,
                year=year,
                source=clean_text(item.get("source") or "online cache"),
                source_url=clean_text(item.get("source_url") or ""),
            )
    return cache


def write_impact_factor_cache(path: Path, cache: dict[str, JournalImpactFactor], journal_names: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": 1,
        "updated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "journals": {
            key: {
                "journal": journal_names.get(key, key),
                "impact_factor": value.impact_factor,
                "year": value.year,
                "source": value.source,
                "source_url": value.source_url,
            }
            for key, value in sorted(cache.items())
        },
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def impact_factor_candidate_urls(journal: str) -> list[tuple[str, str]]:
    quoted = urllib.parse.quote(journal)
    plus = urllib.parse.quote_plus(journal)
    return [
        ("LetPub", f"https://www.letpub.com.cn/index.php?page=journalapp&view=search&searchname={plus}"),
        ("ScholarScope", f"https://www.scholarscope.cn/journal/search?search={quoted}"),
        ("Google-style publisher search", f"https://www.google.com/search?q={plus}+journal+impact+factor"),
    ]


def parse_impact_factor_from_text(text: str) -> tuple[str, str] | None:
    cleaned = clean_text(re.sub(r"<[^>]+>", " ", html.unescape(text)))
    if not re.search(r"(impact factor|影响因子|JIF|Journal Impact Factor)", cleaned, flags=re.IGNORECASE):
        return None
    year_match = re.search(r"\b(20[12]\d)(?:[-/](20[12]\d))?\b", cleaned)
    if not year_match:
        return None
    year = year_match.group(0)
    patterns = [
        r"(?:impact factor|影响因子|JIF|Journal Impact Factor)[^\d]{0,80}(\d+(?:\.\d+)?)",
        r"(\d+(?:\.\d+)?)[^\d]{0,80}(?:impact factor|影响因子|JIF|Journal Impact Factor)",
    ]
    for pattern in patterns:
        match = re.search(pattern, cleaned, flags=re.IGNORECASE)
        if not match:
            continue
        value = match.group(1)
        try:
            number = float(value)
        except ValueError:
            continue
        if 0 <= number <= 100:
            return value.rstrip("0").rstrip(".") if "." in value else value, year
    return None


def lookup_online_impact_factor(journal: str) -> JournalImpactFactor | None:
    for source, url in impact_factor_candidate_urls(journal):
        try:
            parsed = parse_impact_factor_from_text(http_text(url, timeout=30, retries=1))
        except (urllib.error.URLError, urllib.error.HTTPError, http.client.IncompleteRead, TimeoutError, OSError):
            continue
        if parsed:
            impact_factor, year = parsed
            return JournalImpactFactor(impact_factor=impact_factor, year=year, source=source, source_url=url)
    return None


def journal_accessions(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for row in rows:
        accession = clean_text(row.get("accession") or "")
        for journal in (row.get("source_journals") or "").split(";"):
            journal = clean_text(journal)
            if journal:
                result.setdefault(journal, [])
                if accession and accession not in result[journal]:
                    result[journal].append(accession)
    return result


def apply_impact_factor_lookup(
    rows: list[dict[str, str]],
    *,
    cache_path: Path,
    online_lookup: bool,
    warnings: list[str],
    local_jcr_dir: Path = DEFAULT_LOCAL_JCR_DIR,
) -> ImpactFactorAudit:
    cache = load_impact_factor_cache(cache_path)
    local_jcr = load_local_jcr_impact_factors(local_jcr_dir)
    journal_names = {key: key for key in cache}
    audit = ImpactFactorAudit({}, {}, {}, {}, {}, {}, [])
    changed_cache = False
    accessions_by_journal = journal_accessions(rows)

    for journal, accessions in sorted(accessions_by_journal.items(), key=lambda item: normalize_journal_name(item[0])):
        key = normalize_journal_name(journal)
        journal_names[key] = journal
        if is_preprint_or_non_jcr_journal(journal):
            audit.skipped[journal] = accessions
            continue
        impact = local_jcr.get(key)
        if impact:
            audit.local_jcr_matches[journal] = impact
            continue
        impact = ONLINE_IMPACT_FACTORS.get(key)
        if impact:
            audit.static_matches[journal] = impact
            continue
        impact = cache.get(key)
        if impact:
            audit.cache_matches[journal] = impact
            continue
        if online_lookup:
            impact = lookup_online_impact_factor(journal)
            if impact:
                cache[key] = impact
                audit.online_matches[journal] = impact
                changed_cache = True
                continue
        audit.missing[journal] = accessions

    if changed_cache or (online_lookup and not cache_path.exists()):
        write_impact_factor_cache(cache_path, cache, journal_names)

    for row in rows:
        journals = [clean_text(journal) for journal in (row.get("source_journals") or "").split(";") if clean_text(journal)]
        impact = source_impact_factor(journals, cache=cache, local_jcr=local_jcr)
        if impact:
            row["source_journal_impact_factor"] = impact.impact_factor
            row["source_journal_impact_factor_year"] = impact.year

    for journal, accessions in audit.missing.items():
        message = f"missing impact factor for {journal} ({', '.join(accessions[:5])})"
        warnings.append(message)
        audit.warnings.append(message)
    return audit


def impact_factor_audit_path(monthly_dir: Path, today: dt.date | None = None) -> Path:
    stamp = (today or dt.date.today()).strftime("%Y-%m")
    return monthly_dir / f"impact_factor_audit_{stamp}.md"


def write_impact_factor_audit(path: Path, audit: ImpactFactorAudit) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Impact factor audit {dt.date.today().isoformat()}",
        "",
        "Impact factors in CSV are informational only. Source URLs are kept here or in data/online_impact_factor_cache.json, not in CSV columns.",
        "",
    ]
    sections = [
        ("Local JCR XLSX matches", audit.local_jcr_matches),
        ("Static mapping matches", audit.static_matches),
        ("Cache matches", audit.cache_matches),
        ("New online matches", audit.online_matches),
    ]
    for title, matches in sections:
        lines.extend([f"## {title}", ""])
        if matches:
            for journal, impact in sorted(matches.items()):
                source = f"; {impact.source_url}" if impact.source_url else ""
                lines.append(f"- {journal}: {impact.impact_factor} ({impact.year}); {impact.source}{source}")
        else:
            lines.append("- None")
        lines.append("")
    lines.extend(["## Missing impact factors", ""])
    if audit.missing:
        for journal, accessions in sorted(audit.missing.items()):
            query = urllib.parse.quote_plus(journal)
            lines.append(f"- {journal}: {', '.join(accessions)}")
            lines.append(f"  - LetPub: https://www.letpub.com.cn/index.php?page=journalapp&view=search&searchname={query}")
            lines.append(f"  - ScholarScope: https://www.scholarscope.cn/journal/search?search={urllib.parse.quote(journal)}")
    else:
        lines.append("- None")
    lines.extend(["", "## Skipped non-JCR/preprint sources", ""])
    if audit.skipped:
        for journal, accessions in sorted(audit.skipped.items()):
            lines.append(f"- {journal}: {', '.join(accessions)}")
    else:
        lines.append("- None")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fetch_source_journals(pubmed_ids: list[str]) -> list[str]:
    summaries = ncbi_pubmed_summary(pubmed_ids)
    journals: list[str] = []
    for pubmed_id in pubmed_ids:
        item = summaries.get(pubmed_id, {})
        journal = clean_text(item.get("fulljournalname") or item.get("source") or "")
        if journal:
            journals.append(journal)
    return list(dict.fromkeys(journals))


def europe_pmc_search(accession: str) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {
            "query": accession.upper(),
            "format": "json",
            "pageSize": "1000",
            "synonym": "false",
        }
    )
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{params}"
    return http_json(url, timeout=60, retries=3)


def europe_pmc_reuse_ids(accession: str) -> set[str]:
    data = europe_pmc_search(accession)
    results = data.get("resultList", {}).get("result", [])
    ids: set[str] = set()
    for item in results:
        source = clean_text(item.get("source") or "EPMC")
        identifier = clean_text(item.get("pmid") or item.get("pmcid") or item.get("id") or "")
        if identifier:
            ids.add(f"{source}:{identifier}")
    if not ids:
        hit_count = positive_int(data.get("hitCount")) or 0
        ids.update(f"EPMC:{accession.upper()}:{index}" for index in range(hit_count))
    return ids


def ncbi_reuse_ids(accession: str) -> set[str]:
    term = f'"{accession.upper()}"[All Fields]'
    pubmed_result = ncbi_esearch("pubmed", term, retmax=1000)
    pubmed_ids = {f"MED:{value}" for value in pubmed_result.get("idlist", [])}

    pmc_result = ncbi_esearch("pmc", term, retmax=1000)
    pmc_ids = {f"PMC:{value}" for value in pmc_result.get("idlist", [])}
    return pubmed_ids | pmc_ids


def count_gse_reuse(accession: str, source_pubmed_ids: list[str]) -> GseReuseCounts:
    try:
        ids = europe_pmc_reuse_ids(accession)
    except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError):
        ids = ncbi_reuse_ids(accession)
    if not ids:
        ids = ncbi_reuse_ids(accession)
    source_ids = {f"MED:{pubmed_id}" for pubmed_id in source_pubmed_ids}
    return GseReuseCounts(total=max(0, len(ids - source_ids)))


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
    if any(term.casefold() in folded for term in RA_TERMS if term.casefold() != "ra"):
        return "ra"
    if re.search(r"(?<![a-z0-9])ra[- _]?fls(?![a-z0-9])", folded) or "mh7a" in folded:
        return "ra"
    if re.search(r"\bra\b", folded):
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


def extract_gpl_accessions(soft_text: str, summary_text: str) -> str:
    accessions: list[str] = []
    for platform in parse_soft_blocks(soft_text, "^PLATFORM"):
        accession = platform.get("accession", [""])[0].upper()
        if accession:
            accessions.append(accession)
    if not accessions:
        accessions = sorted({match.upper() for match in GPL_RE.findall(summary_text + "\n" + soft_text)})
    return "; ".join(dict.fromkeys(accessions))


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


def score_row(
    *,
    text: str,
    organism: str,
    technology: str,
    counts: SampleCounts,
    tissue: str,
    clinical_hint: str,
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


def row_for_record(
    record: GeoRecord,
    *,
    source_journals: list[str] | None = None,
    reuse_counts: GseReuseCounts | None = None,
) -> dict[str, str] | None:
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
    gpl_accessions = extract_gpl_accessions(soft, merged_text)
    tissue = infer_tissue(merged_text)
    clinical = yes_no_hint(merged_text, CLINICAL_TERMS)
    treatment = yes_no_hint(merged_text, TREATMENT_TERMS)
    score, score_reason = score_row(
        text=merged_text,
        organism=record.organism,
        technology=technology,
        counts=counts,
        tissue=tissue,
        clinical_hint=clinical,
    )
    notes = ["auto-collected from GEO metadata"]
    notes.extend(counts.notes)
    if not gpl_accessions:
        notes.append("manual review needed for GPL platform")

    journals = source_journals or []
    impact_factor = source_impact_factor(journals)
    reuse = reuse_counts or GseReuseCounts(total=0)
    today = dt.date.today().isoformat()
    return {
        "accession": record.accession,
        "title": record.title,
        "organism": record.organism,
        "technology": technology,
        "geo_study_type": geo_study_type or record.geo_study_type,
        "gpl_accessions": gpl_accessions,
        "sample_count": counts.sample_count,
        "sample_count_confidence": counts.confidence,
        "score": str(score),
        "score_reason": score_reason,
        "ra_relevance": ra_reason,
        "tissue": tissue,
        "case_control_hint": case_control_hint(counts, merged_text),
        "treatment_hint": treatment,
        "clinical_info_hint": clinical,
        "source_journals": "; ".join(journals),
        "source_journal_impact_factor": impact_factor.impact_factor if impact_factor else "",
        "source_journal_impact_factor_year": impact_factor.year if impact_factor else "",
        "gse_reuse_total_count": str(reuse.total),
        "publication_date": record.publication_date,
        "last_update_date": record.last_update_date,
        "url": record.url,
        "search_date": today,
        "manual_review": "NULL",
        "notes": "; ".join(dict.fromkeys(note for note in notes if note)),
    }


def quoted_search_terms(terms: list[str]) -> str:
    return "(" + " OR ".join(f'"{term}"' for term in terms) + ")"


def build_search_term() -> str:
    disease = quoted_search_terms([MESH_HEADING, *MESH_ENTRY_TERMS, *CLINICAL_SEARCH_TERMS[:-3]])
    human = '("Homo sapiens"[Organism] OR human)'
    expression = '("Expression profiling by array" OR "Expression profiling by high throughput sequencing" OR "RNA-seq" OR "transcriptome sequencing")'
    return f"{disease} AND {human} AND {expression}"


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
        source_pubmed_ids = extract_source_pubmed_ids(record.soft_text)
        source_journals: list[str] = []
        try:
            source_journals = fetch_source_journals(source_pubmed_ids)
        except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
            warnings.append(f"{record.accession}: failed to fetch source PubMed metadata: {exc}")
        reuse_counts = GseReuseCounts(total=0)
        try:
            reuse_counts = count_gse_reuse(record.accession, source_pubmed_ids)
        except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
            warnings.append(f"{record.accession}: failed to count GSE reuse mentions: {exc}")
        row = row_for_record(record, source_journals=source_journals, reuse_counts=reuse_counts)
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
                existing_row = {field: row.get(field, "") for field in FIELDS}
                existing_row["manual_review"] = normalize_manual_review(row.get("manual_review"))
                rows[accession] = existing_row
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
        old = existing.get(accession, {})
        combined = {field: row.get(field, "") for field in FIELDS}
        combined["manual_review"] = normalize_manual_review(combined.get("manual_review"))
        if old:
            preserve_more_complete_sample_count(old, combined)
            combined["manual_review"] = normalize_manual_review(old.get("manual_review"))
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
        for row in rows:
            output_row = {field: row.get(field, "") for field in FIELDS}
            output_row["manual_review"] = normalize_manual_review(row.get("manual_review"))
            writer.writerow(output_row)


def write_monthly_update(latest_file: Path, monthly_dir: Path, rows: list[dict[str, str]]) -> Path:
    monthly_rows = sorted(rows, key=lambda row: (-int(row.get("score") or 0), row.get("accession", "")))
    archive_file = monthly_archive_path(monthly_dir)
    write_csv(latest_file, monthly_rows)
    write_csv(archive_file, monthly_rows)
    return archive_file


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = [{field: row.get(field, "") for field in FIELDS} for row in reader]
        for row in rows:
            row["manual_review"] = normalize_manual_review(row.get("manual_review"))
        return rows


def citation_fields_for_accession(
    accession: str,
    *,
    warnings: list[str],
) -> dict[str, str]:
    source_pubmed_ids: list[str] = []
    source_journals: list[str] = []
    try:
        source_pubmed_ids = extract_source_pubmed_ids(fetch_soft(accession))
    except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
        warnings.append(f"{accession}: failed to fetch SOFT text for citation fields: {exc}")

    if source_pubmed_ids:
        try:
            source_journals = fetch_source_journals(source_pubmed_ids)
        except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
            warnings.append(f"{accession}: failed to fetch source PubMed metadata: {exc}")

    reuse_counts: GseReuseCounts | None = None
    try:
        reuse_counts = count_gse_reuse(accession, source_pubmed_ids)
    except (urllib.error.URLError, http.client.IncompleteRead, TimeoutError, OSError) as exc:
        warnings.append(f"{accession}: failed to count GSE reuse mentions: {exc}")

    impact_factor = source_impact_factor(source_journals)
    return {
        "source_journals": "; ".join(source_journals),
        "source_journal_impact_factor": impact_factor.impact_factor if impact_factor else "",
        "source_journal_impact_factor_year": impact_factor.year if impact_factor else "",
        "gse_reuse_total_count": str(reuse_counts.total) if reuse_counts else "",
    }


def enrich_existing_csv_files(
    paths: list[Path],
    *,
    force: bool = False,
    accessions: set[str] | None = None,
    verbose: bool = True,
) -> tuple[int, list[str]]:
    warnings: list[str] = []
    cache: dict[str, dict[str, str]] = {}
    total_rows = 0

    for path in paths:
        rows = read_csv_rows(path)
        changed = False
        for index, row in enumerate(rows, start=1):
            accession = (row.get("accession") or "").upper()
            if not accession:
                continue
            if accessions is not None and accession not in accessions:
                continue
            if row.get("gse_reuse_total_count", "") and not force:
                cache.setdefault(
                    accession,
                    {
                        "source_journals": row.get("source_journals", ""),
                        "source_journal_impact_factor": row.get("source_journal_impact_factor", ""),
                        "source_journal_impact_factor_year": row.get("source_journal_impact_factor_year", ""),
                        "gse_reuse_total_count": row.get("gse_reuse_total_count", ""),
                    },
                )
                continue
            if accession not in cache:
                cache[accession] = citation_fields_for_accession(accession, warnings=warnings)
            row.update(cache[accession])
            changed = True
            write_csv(path, rows)
            if verbose:
                print(f"Updated {path} row {index}/{len(rows)} ({accession})", flush=True)
        if not changed:
            write_csv(path, rows)
        total_rows += len(rows)
    return total_rows, warnings


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
    parser.add_argument("--online-impact-factor-lookup", action="store_true", help="Look up missing source journal impact factors online and cache successful matches.")
    parser.add_argument("--impact-factor-cache", type=Path, default=DEFAULT_IMPACT_FACTOR_CACHE, help="JSON cache for online source journal impact factors.")
    parser.add_argument("--enrich-existing-csv", type=Path, nargs="+", default=[], help="Update existing CSV file(s) with source journal, impact factor, and GSE reuse fields, without running a GEO search.")
    parser.add_argument("--force-enrich-existing-csv", action="store_true", help="Refresh citation fields even when existing reuse counts are already populated.")
    parser.add_argument("--enrich-accession", nargs="+", default=[], help="Restrict --enrich-existing-csv updates to specific GSE accessions.")
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
    if args.enrich_existing_csv and (args.backfill_quarterly_local or args.write_monthly_update or args.start_date or args.end_date or args.since_days):
        raise SystemExit("--enrich-existing-csv cannot be combined with search, backfill, or monthly-update options")
    if (args.force_enrich_existing_csv or args.enrich_accession) and not args.enrich_existing_csv:
        raise SystemExit("--force-enrich-existing-csv and --enrich-accession require --enrich-existing-csv")


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
    if "manual_review" not in (reader.fieldnames or []):
        raise RuntimeError("CSV is missing manual_review column")
    invalid_reviews = [row.get("accession", "") for row in rows if row.get("manual_review") not in MANUAL_REVIEW_VALUES]
    if invalid_reviews:
        raise RuntimeError(f"manual_review must be NULL, YES, or NO: {', '.join(invalid_reviews[:10])}")
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
    audit = apply_impact_factor_lookup(
        final_rows,
        cache_path=args.impact_factor_cache,
        online_lookup=args.online_impact_factor_lookup,
        warnings=warnings,
    )
    if rows:
        monthly_lookup_warnings: list[str] = []
        apply_impact_factor_lookup(
            rows,
            cache_path=args.impact_factor_cache,
            online_lookup=False,
            warnings=monthly_lookup_warnings,
        )
    write_csv(args.out, final_rows)
    if args.write_monthly_update:
        archive_file = write_monthly_update(args.latest_file, args.monthly_dir, rows)
        audit_file = impact_factor_audit_path(args.monthly_dir)
        write_impact_factor_audit(audit_file, audit)
        print(f"Wrote monthly update rows to {args.latest_file} and {archive_file}")
        print(f"Wrote impact factor audit to {audit_file}")
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
    if args.enrich_existing_csv:
        accessions = {accession.upper() for accession in args.enrich_accession} if args.enrich_accession else None
        total_rows, warnings = enrich_existing_csv_files(
            args.enrich_existing_csv,
            force=args.force_enrich_existing_csv,
            accessions=accessions,
        )
        print(f"Updated {total_rows} rows across {len(args.enrich_existing_csv)} CSV file(s)")
        if warnings:
            print("Warnings:")
            for warning in warnings:
                print(f"- {warning}")
        return 0
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
