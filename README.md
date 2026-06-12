# Human RA GEO transcriptome ledger

This project builds a CSV ledger of human rheumatoid arthritis expression transcriptome datasets from NCBI GEO.

It includes GEO Series records for:

- `Expression profiling by high throughput sequencing` as `bulk RNA-seq`
- `Expression profiling by array` as `expression array`

It excludes obvious single-cell, single-nucleus, spatial, methylation, ATAC-seq, ChIP-seq, miRNA, and genotyping records. The disease match is intentionally strict: standalone `RA` is not enough because it often means retinoic acid.

## Run

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200
```

For larger runs, GSM sample metadata is fetched concurrently so `sample_count` can show RA/control counts. You can tune that behavior:

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200 --max-samples-per-series 250 --sample-workers 4
```

Write to a custom file:

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200 --out .\data\ra_geo_transcriptome_datasets.csv
```

Search recent records only:

```powershell
python .\scripts\build_ra_geo_ledger.py --since-days 30 --retmax 200
```

Run a fixed publication-date window:

```powershell
python .\scripts\build_ra_geo_ledger.py --start-date 2016-01-01 --end-date 2016-03-31 --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --max-samples-per-series 250 --sample-workers 4 --validate-csv
```

The fixed-window state file is:

```text
data/ra_geo_backfill_state.json
```

If a completed window should be re-run, add `--force-window`.

Run the local quarterly backfill from 2016 to the current date:

```powershell
python .\scripts\build_ra_geo_ledger.py --backfill-quarterly-local --start-date 2016-01-01 --end-date 2026-06-12 --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --max-samples-per-series 250 --sample-workers 4 --validate-csv
```

Use NCBI identity parameters to reduce rate-limit risk:

```powershell
python .\scripts\build_ra_geo_ledger.py --ncbi-email you@example.com --ncbi-api-key YOUR_KEY
```

## Output

The main CSV is:

```text
data/ra_geo_transcriptome_datasets.csv
```

Each row is one GEO Series accession. The `sample_count` field is human-readable and directly shows RA/control counts when they can be parsed, for example:

```text
RA=12; control=8; total=20
RA=10; control=6; other/unclear=4; total=20
RA=unknown; control=unknown; total=24
```

The CSV keeps the GEO `title` but omits the long GEO `summary` field to keep the ledger readable. Use the `url` column to open the source GEO page when the full abstract/summary is needed.

`sample_count_confidence` is:

- `high`: GSM metadata supports per-sample RA/control parsing.
- `medium`: group counts are inferred from GSE-level text.
- `low`: only total samples are available or grouping is unclear.

## Score

`score` is a 0-100 priority score for deciding which datasets to inspect first:

- +20 human study is clear
- +20 RA samples are clear
- +15 control samples are clear
- +15 `sample_count` distinguishes RA and control
- +10 total samples >= 20
- +10 processed/normalized/counts/matrix/supplementary data is suggested
- +5 clinical metadata is suggested
- +5 tissue or cell source is clear
- -20 RA/control groups are not distinguishable
- -15 mixed disease cohort with unclear RA subset
- -15 platform or data type is unclear

The CSV is sorted by `score` descending, then `accession`.

## Notes

GEO metadata is not fully standardized. The script is conservative: when group labels are unclear, it writes `unknown` rather than inventing RA/control counts, lowers the score, and adds a manual-review note.

## Historical Backfill

Run the full `2016-01-01` to present backfill locally with `--backfill-quarterly-local`. The script runs one quarter at a time and writes state after each completed window, so interrupted runs can continue without repeating completed quarters.

The existing state file already records completed windows. By default, completed windows are skipped. Use `--force-window` only when a fixed window must be re-run.

The quarterly sequence is:

- `2016-01-01` to `2016-03-31`
- `2016-04-01` to `2016-06-30`
- `2016-07-01` to `2016-09-30`
- `2016-10-01` to `2016-12-31`

GitHub Actions is not used for historical backfill.

## Monthly Update

The scheduled workflow is:

```text
.github/workflows/ra-geo-ledger-monthly.yml
```

It runs on the first day of each month at 01:00 UTC, which is 09:00 Beijing time, and searches the recent 30-day window.

The monthly update writes:

- `data/ra_geo_transcriptome_datasets.csv`: long-term total ledger
- `data/ra_geo_latest_monthly_update.csv`: latest monthly update table
- `data/monthly/ra_geo_monthly_update_YYYY-MM.csv`: archived monthly update table

The monthly update tables contain all records that pass filtering in that run's recent 30-day search window, not only newly discovered accessions.

Optional GitHub Actions secrets:

- `NCBI_EMAIL`
- `NCBI_API_KEY`
