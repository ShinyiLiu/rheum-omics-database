import datetime as dt
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

DISEASE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DISEASE_ROOT))

from scripts.build_as_geo_ledger import (
    FIELDS,
    GeoRecord,
    GseReuseCounts,
    JournalImpactFactor,
    SampleCounts,
    apply_impact_factor_lookup,
    build_sample_count,
    build_search_term,
    classify_sample_group,
    enrich_existing_csv_files,
    infer_study_type,
    merge_rows,
    monthly_archive_path,
    ncbi_date,
    normalize_manual_review,
    row_for_record,
    quarterly_windows,
    source_impact_factor,
    state_is_complete,
    score_row,
    strict_ra_match,
    window_key,
    write_impact_factor_audit,
    write_monthly_update,
)


def write_sample_jcr_xlsx(path: Path) -> None:
    sheet = """<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetData>
    <row r="3">
      <c r="B3" t="inlineStr"><is><t>期刊名称 Journal Name</t></is></c>
      <c r="G3" t="inlineStr"><is><t>2025影响因子 Impact Factor</t></is></c>
    </row>
    <row r="5">
      <c r="B5" t="inlineStr"><is><t>JOURNAL OF TRANSLATIONAL MEDICINE</t></is></c>
      <c r="G5" t="inlineStr"><is><t>7.5</t></is></c>
    </row>
  </sheetData>
</worksheet>"""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


class LedgerRulesTest(unittest.TestCase):
    def test_strict_disease_match(self):
        self.assertTrue(strict_ra_match("ankylosing spondylitis patient blood RNA-seq")[0])
        self.assertFalse(strict_ra_match("unrelated healthy donor blood RNA-seq")[0])

    def test_sample_group_classification(self):
        self.assertEqual(classify_sample_group("ankylosing spondylitis patient sample"), "ra")
        self.assertEqual(classify_sample_group("healthy control PBMC"), "control")

    def test_sample_count_format(self):
        samples = [
            {"!Sample_title": ["ankylosing spondylitis patient 1"]},
            {"!Sample_title": ["ankylosing spondylitis patient 2"]},
            {"!Sample_title": ["healthy control 1"]},
        ]
        counts = build_sample_count(samples, None, "")
        self.assertEqual(counts.sample_count, "AS=2; control=1; total=3")
        self.assertEqual(counts.confidence, "high")


    def test_mesh_search_term_includes_heading_and_entry_terms(self):
        term = build_search_term()
        self.assertIn('"Spondylitis, Ankylosing"', term)
        self.assertIn('"Ankylosing Spondylitis"', term)
        self.assertIn('"Bechterew Disease"', term)
        self.assertNotIn('"Axial Spondyloarthritis"', term)

    def test_mesh_heading_and_entry_terms_match_disease(self):
        self.assertTrue(strict_ra_match("Spondylitis, Ankylosing patient")[0])
        self.assertTrue(strict_ra_match("Bechterew Disease patient")[0])

    def test_study_type_includes_array(self):
        self.assertEqual(infer_study_type("Expression profiling by array")[0], "expression array")
        self.assertEqual(infer_study_type("Expression profiling by high throughput sequencing")[0], "bulk RNA-seq")

    def test_score_rewards_clear_counts(self):
        score, reason = score_row(
            text="ankylosing spondylitis healthy control clinical blood",
            organism="Homo sapiens",
            technology="bulk RNA-seq",
            counts=SampleCounts("AS=12; control=8; total=20", "high", 20, 12, 8, 0, []),
            tissue="blood",
            clinical_hint="yes",
        )
        self.assertEqual(score, 90)
        self.assertIn("AS/control counts shown", reason)

    def test_score_penalizes_no_parsed_cases(self):
        score, reason = score_row(
            text="ankylosing spondylitis blood",
            organism="Homo sapiens",
            technology="bulk RNA-seq",
            counts=SampleCounts("AS=0; control=0; other/unclear=8; total=8", "low", 8, 0, 0, 8, []),
            tissue="blood",
            clinical_hint="no",
        )
        self.assertLess(score, 50)
        self.assertIn("no AS samples parsed", reason)

    def test_fixed_window_helpers(self):
        start = dt.date(2016, 1, 1)
        end = dt.date(2016, 3, 31)
        self.assertEqual(ncbi_date(start), "2016/01/01")
        self.assertEqual(window_key(start, end), "2016-01-01..2016-03-31")

    def test_state_complete_detection(self):
        state = {"windows": {"2016-01-01..2016-03-31": {"status": "completed_with_warnings"}}}
        self.assertTrue(state_is_complete(state, "2016-01-01..2016-03-31"))
        self.assertFalse(state_is_complete(state, "2016-04-01..2016-06-30"))

    def test_merge_rows_deduplicates_and_sorts_by_score(self):
        rows = merge_rows(
            {"GSE1": {"accession": "GSE1", "score": "10", "notes": "manual note"}},
            [{"accession": "GSE2", "score": "80", "notes": "new"}, {"accession": "GSE1", "score": "90", "notes": "updated"}],
            include_existing=True,
        )
        self.assertEqual([row["accession"] for row in rows], ["GSE1", "GSE2"])
        self.assertIn("manual note", rows[0]["notes"])

    def test_manual_review_values_and_preservation(self):
        self.assertEqual([normalize_manual_review(value) for value in ("NULL", "YES", "NO")], ["NULL", "YES", "NO"])
        with self.assertRaises(ValueError):
            normalize_manual_review("MAYBE")
        rows = merge_rows(
            {"GSE1": {"accession": "GSE1", "score": "10", "manual_review": "YES"}},
            [{"accession": "GSE1", "score": "90", "manual_review": "NULL"}],
            include_existing=False,
        )
        self.assertEqual(rows[0]["manual_review"], "YES")

    def test_quarterly_windows(self):
        windows = quarterly_windows(dt.date(2016, 1, 1), dt.date(2016, 7, 15))
        self.assertEqual(windows[-1], (dt.date(2016, 7, 1), dt.date(2016, 7, 15)))

    def test_monthly_archive_path(self):
        path = monthly_archive_path(Path("data/monthly"), dt.date(2026, 6, 12))
        self.assertEqual(str(path).replace("\\", "/"), "data/monthly/as_geo_monthly_update_2026-06.csv")

    def test_write_monthly_update_allows_empty_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            latest = root / "latest.csv"
            archive = write_monthly_update(latest, root / "monthly", [])
            self.assertTrue(latest.exists())
            self.assertTrue(archive.exists())
            self.assertEqual(latest.read_text(encoding="utf-8").splitlines()[0], ",".join(FIELDS))

    def test_relevance_column_name(self):
        self.assertIn("manual_review", FIELDS)
        self.assertIn("as_relevance", FIELDS)
        self.assertNotIn("ra_relevance", FIELDS)

    def test_row_without_source_pubmed_ids_leaves_source_fields_blank(self):
        record = GeoRecord(
            accession="GSE1",
            title="ankylosing spondylitis blood RNA-seq",
            summary="Homo sapiens Expression profiling by high throughput sequencing processed counts blood",
            organism="Homo sapiens",
            geo_study_type="Expression profiling by high throughput sequencing",
            publication_date="2026/01/01",
            last_update_date="",
            total_samples=3,
            url="https://example.test/GSE1",
            soft_text="\n".join([
                "^SERIES = GSE1",
                "!Series_title = ankylosing spondylitis blood RNA-seq",
                "^SAMPLE = GSM1",
                "!Sample_title = ankylosing spondylitis patient 1",
                "^SAMPLE = GSM2",
                "!Sample_title = ankylosing spondylitis patient 2",
                "^SAMPLE = GSM3",
                "!Sample_title = healthy control 1",
            ]),
        )
        row = row_for_record(record)
        self.assertIsNotNone(row)
        assert row is not None
        self.assertEqual(row["source_journals"], "")
        self.assertEqual(row["gse_reuse_total_count"], "0")
        self.assertIn("as_relevance", row)

    def test_row_writes_reuse_counts_and_impact_factor(self):
        record = GeoRecord(
            accession="GSE1",
            title="ankylosing spondylitis blood RNA-seq",
            summary="Homo sapiens Expression profiling by high throughput sequencing processed counts blood",
            organism="Homo sapiens",
            geo_study_type="Expression profiling by high throughput sequencing",
            publication_date="2026/01/01",
            last_update_date="",
            total_samples=3,
            url="https://example.test/GSE1",
            soft_text="\n".join([
                "^SERIES = GSE1",
                "!Series_pubmed_id = 12345",
                "^SAMPLE = GSM1",
                "!Sample_title = ankylosing spondylitis patient 1",
                "^SAMPLE = GSM2",
                "!Sample_title = ankylosing spondylitis patient 2",
                "^SAMPLE = GSM3",
                "!Sample_title = healthy control 1",
            ]),
        )
        row = row_for_record(record, source_journals=["Arthritis & Rheumatology (Hoboken, N.J.)"], reuse_counts=GseReuseCounts(total=5))
        self.assertIsNotNone(row)
        assert row is not None
        self.assertEqual(row["source_journal_impact_factor"], "10.9")
        self.assertEqual(row["gse_reuse_total_count"], "5")

    def test_enrich_existing_csv_files_updates_new_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            csv_path = root / "ledger.csv"
            csv_path.write_text("accession,title,score,sample_count\nGSE1,Title,100,AS=1; control=1; total=2\n", encoding="utf-8")
            with patch(
                "scripts.build_as_geo_ledger.citation_fields_for_accession",
                return_value={
                    "source_journals": "Arthritis & Rheumatology (Hoboken, N.J.)",
                    "source_journal_impact_factor": "10.9",
                    "source_journal_impact_factor_year": "2024",
                    "gse_reuse_total_count": "5",
                },
            ):
                total_rows, warnings = enrich_existing_csv_files([csv_path], verbose=False)
            self.assertEqual(total_rows, 1)
            self.assertEqual(warnings, [])
            self.assertIn("10.9", csv_path.read_text(encoding="utf-8"))

    def test_online_impact_factor_lookup_updates_cache_and_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"
            rows = [{"accession": "GSE2", "source_journals": "New Journal", "source_journal_impact_factor": "", "source_journal_impact_factor_year": ""}]
            with patch(
                "scripts.build_as_geo_ledger.lookup_online_impact_factor",
                return_value=JournalImpactFactor("4.2", "2025", "LetPub", "https://example.test/letpub"),
            ):
                audit = apply_impact_factor_lookup(rows, cache_path=cache_path, online_lookup=True, warnings=[])
            self.assertEqual(rows[0]["source_journal_impact_factor"], "4.2")
            self.assertTrue(cache_path.exists())
            self.assertIn("New Journal", audit.online_matches)

    def test_local_jcr_impact_factor_skips_online_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_dir = Path(tmp) / "IF reference"
            reference_dir.mkdir()
            write_sample_jcr_xlsx(reference_dir / "sample_jcr.xlsx")
            rows = [{"accession": "GSE6", "source_journals": "Journal of translational medicine", "source_journal_impact_factor": "", "source_journal_impact_factor_year": ""}]
            with patch("scripts.build_as_geo_ledger.lookup_online_impact_factor") as mocked:
                audit = apply_impact_factor_lookup(rows, cache_path=Path(tmp) / "cache.json", online_lookup=True, warnings=[], local_jcr_dir=reference_dir)
            mocked.assert_not_called()
            self.assertEqual(rows[0]["source_journal_impact_factor"], "7.5")
            self.assertEqual(rows[0]["source_journal_impact_factor_year"], "2025")
            self.assertIn("Journal of translational medicine", audit.local_jcr_matches)

    def test_missing_impact_factor_audit_skips_preprint_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = [
                {"accession": "GSE4", "source_journals": "Unknown Journal", "source_journal_impact_factor": "", "source_journal_impact_factor_year": ""},
                {"accession": "GSE5", "source_journals": "bioRxiv : the preprint server for biology", "source_journal_impact_factor": "", "source_journal_impact_factor_year": ""},
            ]
            with patch("scripts.build_as_geo_ledger.lookup_online_impact_factor", return_value=None):
                audit = apply_impact_factor_lookup(rows, cache_path=Path(tmp) / "cache.json", online_lookup=True, warnings=[])
            report = Path(tmp) / "audit.md"
            write_impact_factor_audit(report, audit)
            text = report.read_text(encoding="utf-8")
            self.assertIn("Unknown Journal", audit.missing)
            self.assertIn("Skipped non-JCR/preprint sources", text)

    def test_as_strict_scope(self):
        self.assertTrue(strict_ra_match("ankylosing spondylitis patient")[0])
        self.assertFalse(strict_ra_match("AS patient blood sample")[0])
        self.assertFalse(strict_ra_match("axial spondyloarthritis patient")[0])
        self.assertFalse(strict_ra_match("axSpA patient")[0])


if __name__ == "__main__":
    unittest.main()
