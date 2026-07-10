# 人类风湿免疫性疾病转录组数据库

## **数据库介绍**

本数据库整理了 NCBI Gene Expression Omnibus（GEO）中与人类风湿免疫性疾病相关的表达转录组数据集。数据库以 GEO Series（GSE）为基本单位，当前包含以下 4 个疾病子库：

| 缩写 | 疾病 | 目录 | 主数据表 | 当前记录数 |
| --- | --- | --- | --- | ---: |
| RA | 类风湿关节炎（rheumatoid arthritis） | [`RA/`](RA/) | [`ra_geo_transcriptome_datasets.csv`](RA/data/ra_geo_transcriptome_datasets.csv) | 201 |
| SLE | 系统性红斑狼疮（systemic lupus erythematosus，包含狼疮性肾炎） | [`SLE/`](SLE/) | [`sle_geo_transcriptome_datasets.csv`](SLE/data/sle_geo_transcriptome_datasets.csv) | 102 |
| AS | 强直性脊柱炎（ankylosing spondylitis） | [`AS/`](AS/) | [`as_geo_transcriptome_datasets.csv`](AS/data/as_geo_transcriptome_datasets.csv) | 16 |
| pSS | 原发性干燥综合征（primary Sjögren's syndrome） | [`pss/`](pss/) | [`pss_geo_transcriptome_datasets.csv`](pss/data/pss_geo_transcriptome_datasets.csv) | 35 |

> 上表是仓库当前数据快照。各子库最近一次写入的 `search_date` 分别为 RA：2026-06-12，SLE/AS/pSS：2026-07-11；记录数会随数据库更新而变化。

纳入的数据类型为：

- `Expression profiling by high throughput sequencing`，统一标记为 `bulk RNA-seq`；
- `Expression profiling by array`，统一标记为 `expression array`。

脚本采用疾病名称、MeSH 主题词及其同义词进行检索，并在本地执行较严格的疾病相关性过滤。明显属于单细胞/单核转录组、空间转录组、甲基化、ATAC-seq、ChIP-seq、miRNA、基因分型等类型的数据集会被排除。疾病缩写可能存在歧义，因此仅出现孤立的 `RA`、`AS` 或 `pSS` 并不一定会被判定为对应疾病。

当前目录结构如下：

```text
.
├── RA/                  # 类风湿关节炎子库
├── SLE/                 # 系统性红斑狼疮子库
├── AS/                  # 强直性脊柱炎子库
├── pss/                 # 原发性干燥综合征子库
├── IF reference/        # 本地 JCR 期刊影响因子参考表
└── README.md
```

每个疾病目录均包含：

```text
<disease>/
├── README.md            # 该疾病子库的详细说明
├── data/                # 汇总数据、月度数据、回填状态及影响因子审计
├── scripts/             # GEO 数据采集与整理脚本
└── tests/               # 脚本规则和输出测试
```

## **数据库内数据及列名解释**

### 数据文件

每个疾病子库的 `data/` 目录主要包含以下文件：

| 文件 | 说明 |
| --- | --- |
| `*_geo_transcriptome_datasets.csv` | 长期累计的主数据表，每行对应一个 GSE |
| `*_geo_latest_monthly_update.csv` | 最近一次更新任务检索并通过过滤的记录；不等同于“仅新增记录” |
| `monthly/*_geo_monthly_update_YYYY-MM.csv` | 按年月归档的更新结果 |
| `*_geo_backfill_state.json` | 固定时间窗或历史回填任务的完成状态，用于断点续跑和避免重复执行 |
| `online_impact_factor_cache.json` | 已成功查询的在线期刊影响因子缓存（存在时） |
| `monthly/impact_factor_audit_YYYY-MM.md` | 影响因子匹配来源、缺失项及警告的审计记录 |
| `*_geo_accession_import_audit_YYYY-MM-DD.csv` | 外部候选 GSE 清单与总表对账时生成的逐条纳入、排除和重复决策记录（存在时） |

仓库当前主数据表概况：

| 疾病 | GSE 数量 | bulk RNA-seq | expression array | GEO 发布日期范围 |
| --- | ---: | ---: | ---: | --- |
| RA | 201 | 140 | 61 | 2016-01-27 至 2026-06-10 |
| SLE | 102 | 45 | 57 | 2006-12-31 至 2026-06-04 |
| AS | 16 | 10 | 6 | 2008-06-27 至 2025-07-09 |
| pSS | 35 | 16 | 19 | 2007-12-31 至 2024-07-17 |

### CSV 列名

4 个主数据表具有相同的核心结构，仅疾病相关性字段不同：RA、SLE、AS 和 pSS 分别使用 `ra_relevance`、`sle_relevance`、`as_relevance` 和 `pss_relevance`。

| 列名 | 含义 |
| --- | --- |
| `accession` | GEO Series 登录号，例如 `GSE12345` |
| `title` | GEO 数据集标题 |
| `organism` | GEO 元数据中的物种 |
| `technology` | 统一后的技术类型：`bulk RNA-seq` 或 `expression array` |
| `geo_study_type` | GEO 原始研究类型描述 |
| `gpl_accessions` | 数据集关联的 GEO Platform（GPL）登录号 |
| `sample_count` | 可读的病例/对照/其他/总样本数，例如 `RA=12; control=8; total=20` |
| `sample_count_confidence` | 样本分组解析置信度：`high`、`medium` 或 `low` |
| `score` | 0–100 的人工复核优先级评分，越高表示人类病例、对照、样本量和临床信息越明确 |
| `score_reason` | 评分加分项、减分项及原因 |
| `*_relevance` | 该条记录被判定为相应疾病相关的文本匹配依据 |
| `tissue` | 根据 GEO 及样本元数据推断的组织或细胞来源 |
| `case_control_hint` | 是否能识别病例-对照设计，取值为 `yes`/`no` |
| `treatment_hint` | 是否出现治疗、药物或疗效反应相关信息，取值为 `yes`/`no` |
| `clinical_info_hint` | 是否出现疾病活动度、实验室指标或疗效等临床信息，取值为 `yes`/`no` |
| `source_journals` | 与 GEO 源论文 PMID 对应的期刊名称 |
| `source_journal_impact_factor` | 由本地 JCR 表、内置映射、缓存或可核查在线来源匹配的期刊影响因子；无法可靠匹配时留空 |
| `source_journal_impact_factor_year` | 上述影响因子对应的年份或年份范围 |
| `gse_reuse_total_count` | 公开文献中直接提及该 GSE 登录号的去重记录数，是数据集被复用或讨论次数的近似指标 |
| `publication_date` | GEO 数据集发布日期 |
| `last_update_date` | GEO 数据集最近更新日期（若可获得） |
| `url` | GEO 数据集详情页链接 |
| `search_date` | 脚本检索或刷新该记录的日期 |
| `manual_review` | 人工审核状态，只允许 `NULL`、`YES`、`NO`：`NULL` 表示尚未人工审核，`YES` 表示审核通过，`NO` 表示审核未通过 |
| `notes` | 自动采集警告及需要人工复核的说明 |

`sample_count_confidence` 的解释如下：

- `high`：GSM 样本级元数据支持病例/对照逐样本解析；
- `medium`：分组数主要由 GSE 层级文本推断；
- `low`：仅能得到总样本数，或病例/对照分组不清楚。

`score` 用于安排人工复核顺序，不代表数据集的科研质量、证据等级或统计学质量。评分主要奖励明确的人类研究、疾病样本、对照样本、可解析的分组数量、较大样本量、临床信息和组织来源；混合疾病队列、分组不清或技术平台不明会被扣分。

## **脚本信息**

### 运行环境

- Python 3.10 或更高版本（当前使用 Python 3.12 验证）；
- 脚本仅使用 Python 标准库，无需安装额外 Python 包；
- 在线检索需要访问 NCBI GEO、NCBI E-utilities、Europe PMC 及可能的期刊/影响因子公开页面；
- 为减少 NCBI 限流风险，建议通过 `--ncbi-email` 和 `--ncbi-api-key` 提供 NCBI 身份信息。

### 快速运行

以下命令均在仓库根目录执行：

```powershell
# RA
python .\RA\scripts\build_ra_geo_ledger.py --retmax 200 --validate-csv

# SLE
python .\SLE\scripts\build_sle_geo_ledger.py --retmax 200 --validate-csv

# AS
python .\AS\scripts\build_as_geo_ledger.py --retmax 200 --validate-csv

# pSS
python .\pss\scripts\build_pss_geo_ledger.py --retmax 200 --validate-csv
```

脚本默认将结果写入对应疾病目录下的主数据表。以 RA 为例，常用运行方式如下（其他疾病替换目录名、脚本名和文件名前缀即可）：

```powershell
# 检索最近 30 天，并保留主表中的既有记录
python .\RA\scripts\build_ra_geo_ledger.py `
  --since-days 30 `
  --retmax 500 `
  --include-existing `
  --write-monthly-update `
  --validate-csv

# 执行固定发布日期窗口
python .\RA\scripts\build_ra_geo_ledger.py `
  --start-date 2016-01-01 `
  --end-date 2016-03-31 `
  --retmax 500 `
  --include-existing `
  --state-file .\RA\data\ra_geo_backfill_state.json `
  --validate-csv

# 从指定日期开始按季度回填历史记录
python .\RA\scripts\build_ra_geo_ledger.py `
  --backfill-quarterly-local `
  --start-date 2016-01-01 `
  --end-date 2026-06-30 `
  --retmax 500 `
  --include-existing `
  --state-file .\RA\data\ra_geo_backfill_state.json `
  --validate-csv
```

常用参数：

| 参数 | 说明 |
| --- | --- |
| `--retmax N` | GEO ESearch 最多返回的记录数 |
| `--since-days N` | 仅检索最近 N 天发布的记录 |
| `--start-date` / `--end-date` | 设置固定发布日期窗口，格式为 `YYYY-MM-DD` |
| `--include-existing` | 合并并保留当前主表中未出现在本次检索结果内的记录 |
| `--out PATH` | 指定主 CSV 输出路径 |
| `--max-samples-per-series N` | 每个 GSE 最多抓取的 GSM 样本数，默认 250 |
| `--sample-workers N` | GSM 元数据并发抓取线程数，默认 4 |
| `--write-monthly-update` | 同时写入最新月度表和年月归档表 |
| `--online-impact-factor-lookup` | 在线补充缺失的期刊影响因子，并缓存可靠匹配结果 |
| `--enrich-existing-csv PATH...` | 不重新检索 GEO，仅补充既有 CSV 的期刊、影响因子和 GSE 复用字段 |
| `--state-file PATH` | 记录固定时间窗/历史回填状态，以便断点续跑 |
| `--force-window` | 即使状态文件显示已完成，仍重新运行指定时间窗 |
| `--ncbi-email` / `--ncbi-api-key` | NCBI 请求身份信息 |
| `--validate-csv` | 写入后检查 CSV 字段和记录有效性 |

查看某个脚本的完整参数：

```powershell
python .\RA\scripts\build_ra_geo_ledger.py --help
```

运行全部本地测试：

```powershell
python -m unittest discover -s .\RA\tests -p "test_*.py"
python -m unittest discover -s .\SLE\tests -p "test_*.py"
python -m unittest discover -s .\AS\tests -p "test_*.py"
python -m unittest discover -s .\pss\tests -p "test_*.py"
```

### GitHub Actions 月度更新

本仓库使用一个统一的单仓库工作流：

```text
.github/workflows/geo-ledgers-monthly.yml
```

工作流默认在每月 1 日 01:00 UTC（北京时间 09:00）运行。它通过动态矩阵并行更新 RA、SLE、AS、pSS，各任务将更新后的疾病数据作为短期构件上传，最后由统一的汇总任务集中提交并推送一次。这样既缩短总运行时间，也避免多个任务同时写入同一分支产生冲突。

也可通过 `workflow_dispatch` 手动执行，并设置：

- `disease`：选择 `all`、`RA`、`SLE`、`AS` 或 `pSS`；
- `since_days`：检索最近多少天的数据，默认 30；
- `retmax`：每个疾病最多请求的 GEO 记录数，默认 500。

仓库可配置以下可选密钥：`NCBI_EMAIL` 和 `NCBI_API_KEY`。工作流具有 `contents: write` 权限，用于提交更新后的数据表、月度归档、回填状态、影响因子缓存和审计报告。

更详细的疾病匹配规则、评分逻辑和示例请参阅各子库的 [`RA/README.md`](RA/README.md)、[`SLE/README.md`](SLE/README.md)、[`AS/README.md`](AS/README.md) 和 [`pss/README.md`](pss/README.md)。

## **注意事项**

1. GEO 元数据并非完全标准化。脚本无法可靠判断分组时会写入 `unknown`、降低评分并在 `notes` 中提示，而不会猜测病例/对照数。
2. 数据库是自动检索和规则筛选的结果，可能存在漏检、误纳入、重复队列或样本分组解析偏差。用于正式研究前，应通过 `url` 回到 GEO 页面，并结合论文及原始样本信息人工核验。
3. `sample_count` 是根据公开元数据生成的便读摘要；当 GSE 的 GSM 数量超过 `--max-samples-per-series`、网络请求失败或分组标签不规范时，计数可能不完整。
4. `gse_reuse_total_count` 统计公开文献记录中对 GSE 登录号的直接提及，不是 Google Scholar 引用数，也不是源论文的被引次数。
5. 影响因子具有年份属性，且不同来源可能存在差异。无法同时确认数值和年份时字段会留空；预印本平台不按 JCR 期刊处理。具体来源和未匹配项请查看相应的影响因子审计文件。
6. `manual_review` 是人工维护字段。新记录默认为 `NULL`，人工审核后只能改为 `YES` 或 `NO`；不要留空或填写其他值。脚本重新抓取同一 GSE 时会保留主表中已有的人工审核状态，不会用默认 `NULL` 覆盖它。
7. 不带 `--include-existing` 运行脚本时，输出主表只包含本次检索并通过过滤的记录；虽然重合 GSE 的 `manual_review` 会被保留，但未进入本次结果的既有记录仍会从输出中消失。更新长期主表时应谨慎检查命令参数并保留备份。
8. 月度更新表包含该次检索时间窗内所有通过过滤的记录，并不只包含相对主表新出现的 GSE。
9. CSV 和 Markdown 文件采用 UTF-8 编码。Excel 直接打开 CSV 时如出现中文显示异常，请使用“从文本/CSV 导入”并显式选择 UTF-8。
