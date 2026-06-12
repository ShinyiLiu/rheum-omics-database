import unittest
import datetime as dt
import tempfile
from pathlib import Path

from scripts.build_ra_geo_ledger import (
    FIELDS,
    SampleCounts,
    build_sample_count,
    classify_sample_group,
    infer_study_type,
    merge_rows,
    monthly_archive_path,
    ncbi_date,
    quarterly_windows,
    state_is_complete,
    score_row,
    strict_ra_match,
    window_key,
    write_monthly_update,
)


class LedgerRulesTest(unittest.TestCase):
    def test_strict_ra_does_not_match_retinoic_acid(self):
        self.assertFalse(strict_ra_match("retinoic acid response in human cells")[0])
        self.assertTrue(strict_ra_match("rheumatoid arthritis synovium")[0])

    def test_sample_group_classification(self):
        self.assertEqual(classify_sample_group("rheumatoid arthritis patient synovium"), "ra")
        self.assertEqual(classify_sample_group("healthy control PBMC"), "control")
        self.assertEqual(classify_sample_group("retinoic acid treated sample"), "unclear")

    def test_sample_count_format(self):
        samples = [
            {"!Sample_title": ["rheumatoid arthritis patient 1"]},
            {"!Sample_title": ["rheumatoid arthritis patient 2"]},
            {"!Sample_title": ["healthy control 1"]},
        ]
        counts = build_sample_count(samples, None, "")
        self.assertEqual(counts.sample_count, "RA=2; control=1; total=3")
        self.assertEqual(counts.confidence, "high")

    def test_study_type_includes_array(self):
        self.assertEqual(infer_study_type("Expression profiling by array")[0], "expression array")
        self.assertEqual(infer_study_type("Expression profiling by high throughput sequencing")[0], "bulk RNA-seq")

    def test_score_rewards_clear_counts(self):
        score, reason = score_row(
            text="rheumatoid arthritis healthy control processed counts DAS28 synovium",
            organism="Homo sapiens",
            technology="bulk RNA-seq",
            counts=SampleCounts("RA=12; control=8; total=20", "high", 20, 12, 8, 0, []),
            tissue="synovium",
            clinical_hint="yes",
            processed="yes",
        )
        self.assertEqual(score, 100)
        self.assertIn("RA/control counts shown", reason)

    def test_score_penalizes_no_parsed_ra_samples(self):
        score, reason = score_row(
            text="rheumatoid arthritis processed counts synovium",
            organism="Homo sapiens",
            technology="bulk RNA-seq",
            counts=SampleCounts("RA=0; control=0; other/unclear=8; total=8", "low", 8, 0, 0, 8, []),
            tissue="synovium",
            clinical_hint="no",
            processed="yes",
        )
        self.assertLess(score, 50)
        self.assertIn("no RA samples parsed", reason)

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
            {
                "GSE1": {"accession": "GSE1", "score": "10", "notes": "manual note"},
            },
            [
                {"accession": "GSE2", "score": "80", "notes": "new"},
                {"accession": "GSE1", "score": "90", "notes": "updated"},
            ],
            include_existing=True,
        )
        self.assertEqual([row["accession"] for row in rows], ["GSE1", "GSE2"])
        self.assertIn("manual note", rows[0]["notes"])

    def test_merge_preserves_more_complete_sample_count(self):
        rows = merge_rows(
            {
                "GSE1": {
                    "accession": "GSE1",
                    "sample_count": "RA=61; control=104; total=165",
                    "sample_count_confidence": "high",
                    "score": "100",
                    "score_reason": "old complete",
                    "notes": "",
                },
            },
            [
                {
                    "accession": "GSE1",
                    "sample_count": "RA=61; control=103; total=164",
                    "sample_count_confidence": "high",
                    "score": "100",
                    "score_reason": "new partial",
                    "notes": "",
                },
            ],
            include_existing=True,
        )
        self.assertEqual(rows[0]["sample_count"], "RA=61; control=104; total=165")
        self.assertIn("preserved previous sample_count", rows[0]["notes"])

    def test_quarterly_windows(self):
        windows = quarterly_windows(dt.date(2016, 1, 1), dt.date(2016, 7, 15))
        self.assertEqual(
            windows,
            [
                (dt.date(2016, 1, 1), dt.date(2016, 3, 31)),
                (dt.date(2016, 4, 1), dt.date(2016, 6, 30)),
                (dt.date(2016, 7, 1), dt.date(2016, 7, 15)),
            ],
        )

    def test_monthly_archive_path(self):
        path = monthly_archive_path(Path("data/monthly"), dt.date(2026, 6, 12))
        self.assertEqual(str(path).replace("\\", "/"), "data/monthly/ra_geo_monthly_update_2026-06.csv")

    def test_write_monthly_update_allows_empty_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            latest = root / "latest.csv"
            archive = write_monthly_update(latest, root / "monthly", [])
            self.assertTrue(latest.exists())
            self.assertTrue(archive.exists())
            self.assertEqual(latest.read_text(encoding="utf-8").splitlines()[0], ",".join(FIELDS))
            self.assertEqual(archive.read_text(encoding="utf-8").splitlines()[0], ",".join(FIELDS))


if __name__ == "__main__":
    unittest.main()
