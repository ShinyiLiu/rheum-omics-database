# 人类类风湿关节炎 GEO 转录组数据台账

本项目用于构建 NCBI GEO 中人类类风湿关节炎（rheumatoid arthritis，RA）表达转录组数据集的 CSV 台账。

纳入以下 GEO Series 数据：

- `Expression profiling by high throughput sequencing`，统一标记为 `bulk RNA-seq`
- `Expression profiling by array`，统一标记为 `expression array`

GEO 检索包含 MeSH 主题词 `Arthritis, Rheumatoid`、入口词 `Rheumatoid Arthritis`，以及供本地严格过滤使用的临床同义词。疾病匹配较为保守，因为独立出现的 `RA` 在非疾病语境中也可能指视黄酸（retinoic acid）。

本项目排除明显属于单细胞、单核、空间转录组、甲基化、ATAC-seq、ChIP-seq、miRNA 和基因分型的数据。疾病匹配有意采用严格规则：仅出现独立的 `RA` 不足以判定为类风湿关节炎。

## 运行方式

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200
```

运行较大规模任务时，脚本会并发获取 GSM 样本元数据，使 `sample_count` 能够显示 RA/对照样本数。可以通过以下参数调整该行为：

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200 --max-samples-per-series 250 --sample-workers 4
```

写入自定义文件：

```powershell
python .\scripts\build_ra_geo_ledger.py --retmax 200 --out .\data\ra_geo_transcriptome_datasets.csv
```

仅检索最近发布的记录：

```powershell
python .\scripts\build_ra_geo_ledger.py --since-days 30 --retmax 200
```

检索截至指定日期的最近 30 天（例如截至 2026-06-30）：

```powershell
python .\scripts\build_ra_geo_ledger.py --start-date 2026-06-01 --end-date 2026-06-30 --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --write-monthly-update --max-samples-per-series 250 --sample-workers 4 --validate-csv
```

运行固定发布日期窗口：

```powershell
python .\scripts\build_ra_geo_ledger.py --start-date 2016-01-01 --end-date 2016-03-31 --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --max-samples-per-series 250 --sample-workers 4 --validate-csv
```

固定窗口状态文件为：

```text
data/ra_geo_backfill_state.json
```

如需重新运行已完成的窗口，请添加 `--force-window`。

从 2016 年开始按季度在本地回填至当前日期：

```powershell
python .\scripts\build_ra_geo_ledger.py --backfill-quarterly-local --start-date 2016-01-01 --end-date YYYY-MM-DD --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --max-samples-per-series 250 --sample-workers 4 --validate-csv
```

建议提供 NCBI 身份参数以降低触发限流的风险：

```powershell
python .\scripts\build_ra_geo_ledger.py --ncbi-email you@example.com --ncbi-api-key YOUR_KEY
```

## 输出

主 CSV 文件为：

```text
data/ra_geo_transcriptome_datasets.csv
```

每一行对应一个 GEO Series 登录号。`sample_count` 字段采用便于阅读的格式；能够解析分组时会直接显示 RA/对照样本数，例如：

```text
RA=12; control=8; total=20
RA=10; control=6; other/unclear=4; total=20
RA=unknown; control=unknown; total=24
```

CSV 保留 GEO `title`，但省略较长的 GEO `summary` 字段，以保持台账易读。如需查看完整摘要，请通过 `url` 列打开 GEO 来源页面。

CSV 还包含源论文和数据集复用相关字段：

- `source_journals`：根据相关 PMID 从 PubMed 获取的期刊名称。
- `source_journal_impact_factor`：从可公开核查的在线来源获得的期刊影响因子；无法获取时留空。
- `source_journal_impact_factor_year`：影响因子对应年份。
- `gse_reuse_total_count`：公开文献记录中直接提及该 GSE 登录号（如 `GSE12345`）的去重数量。该值可作为 GEO 数据集被后续研究复用或讨论频率的实用近似指标。脚本优先使用 Europe PMC，因为其索引 PubMed、PMC、预印本和全文；失败时使用 NCBI E-utilities。`0` 表示未找到直接提及该登录号的公开文献记录。

影响因子不是从本地 `journal_impact_factors.csv` 文件读取的。无法找到可靠的公开在线来源时，影响因子字段留空。bioRxiv 等预印本平台不作为 JCR 期刊处理，因此其影响因子字段留空。`gse_reuse_total_count` 既不是 Google Scholar 引用次数，也不是源论文的被引次数。

月度更新使用 `--online-impact-factor-lookup` 时，脚本先使用内置期刊映射，再使用持久化在线缓存：

```text
data/online_impact_factor_cache.json
```

如果新期刊不在上述来源中，脚本会尝试公开在线页面，并且仅在能够同时解析数值和年份时写入影响因子。来源名称和 URL 保存在缓存及月度审计报告中，而不写入 CSV 列。审计报告写入：

```text
data/monthly/impact_factor_audit_YYYY-MM.md
```

脚本使用的影响因子来源页面包括：

- Haematologica 期刊介绍页：https://haematologica.org/about
- ASBMB Journal of Biological Chemistry 页面：https://www.asbmb.org/journals/journal-of-biological-chemistry
- Springer Nature Immunologic Research 页面：https://link.springer.com/journal/12026
- Wiley/ACR Arthritis & Rheumatology 页面：https://acrjournals.onlinelibrary.wiley.com/journal/23265205
- Wiley The Journal of Physiology 页面：https://physoc.onlinelibrary.wiley.com/journal/14697793
- PNAS 期刊信息页：https://www.pnas.org/about
- Journal of Investigative Dermatology 期刊信息页：https://www.jidonline.org/
- Nature Communications 页面：https://www.nature.com/ncomms/
- Scientific Reports 页面：https://www.nature.com/srep/
- Frontiers in Immunology 介绍页：https://www.frontiersin.org/journals/immunology/about
- Nature Portfolio 期刊指标页面，包括 Nature Immunology 和 Scientific Data。
- Oxford Academic 期刊指标页面，包括 Clinical and Experimental Immunology 和 The Journal of Immunology。
- Association for Molecular Pathology 的 Journal of Molecular Diagnostics 指标页面。
- Springer Nature 期刊指标页面，包括 Molecular Biology Reports 和 Journal of Translational Medicine。
- 当期刊或出版商页面未提供可直接解析的影响因子时使用 LetPub 页面或搜索结果，涉及 PLoS One、Science Translational Medicine、Immunity、Cellular & Molecular Immunology、Science Advances、International Journal of Rheumatic Diseases、Clinical and Experimental Rheumatology、Clinical and Translational Medicine、Biomedicines、Chinese Medicine、Human Gene Therapy、EBioMedicine、EMBO Reports 和 Genes and Immunity 等期刊。
- 对部分具有可检索影响因子历史的期刊使用科研通/ScholarScope 风格页面，包括 Arthritis Research & Therapy、Journal of Autoimmunity、Communications Biology、Journal of Nanobiotechnology、Journal of Translational Medicine 和 Cell Reports Medicine。

## CSV 列名

- `accession`：GEO Series 登录号，例如 `GSE12345`。
- `title`：GEO Series 标题。
- `organism`：GEO 元数据中的物种。
- `technology`：脚本统一后的技术类型，目前为 `bulk RNA-seq` 或 `expression array`。
- `geo_study_type`：GEO 原始研究类型文本，例如 `Expression profiling by high throughput sequencing`。
- `gpl_accessions`：该 Series 中检测到的 GEO Platform 登录号，例如 `GPL24676`。
- `sample_count`：便于阅读的 RA/对照/样本总数摘要，例如 `RA=12; control=8; total=20`。
- `sample_count_confidence`：样本分组解析置信度：`high`、`medium` 或 `low`。
- `score`：0–100 的人工复核优先级评分，用于决定优先检查哪些数据集。
- `score_reason`：以分号分隔的评分组成和扣分原因。
- `ra_relevance`：说明记录为何被判定为类风湿关节炎相关的文本匹配依据。
- `tissue`：根据 GEO 元数据和样本文本推断的组织或细胞来源。
- `case_control_hint`：是否可识别病例-对照设计的 `yes`/`no` 提示。
- `treatment_hint`：是否出现治疗、药物或疗效反应语境的 `yes`/`no` 提示。
- `clinical_info_hint`：是否出现 DAS28、CRP、ESR、缓解或疗效反应等临床词语的 `yes`/`no` 提示。
- `source_journals`：根据 PubMed 元数据获取的 GEO 源论文 PMID 对应期刊名称。
- `source_journal_impact_factor`：从公开在线来源或在线缓存获得的源期刊影响因子；无法获取时留空。
- `source_journal_impact_factor_year`：影响因子对应年份或年份范围。
- `gse_reuse_total_count`：公开文献记录中直接提及该 GSE 登录号的去重数量。
- `publication_date`：GEO 发布日期。
- `last_update_date`：GEO 最近更新日期（如可获得）。
- `url`：GEO 登录号详情页链接。
- `search_date`：脚本生成或刷新该行的日期。
- `manual_review`：人工审核状态。只允许 `NULL`（尚未审核）、`YES`（人工审核通过）和 `NO`（人工审核未通过）。新记录默认为 `NULL`；后续脚本更新同一 GSE 时会保留现有状态。
- `notes`：自动采集说明和人工复核警告。

`sample_count_confidence` 的含义：

- `high`：GSM 元数据支持逐样本解析 RA/对照分组。
- `medium`：分组数量由 GSE 层级文本推断。
- `low`：只能获得总样本数，或分组信息不明确。

## 评分规则

`score` 是用于确定数据集人工检查优先级的 0–100 分评分：

- +20：明确为人类研究
- +20：RA 样本明确
- +15：对照样本明确
- +15：`sample_count` 能区分 RA 与对照
- +10：总样本数 >= 20
- +5：可能包含临床元数据
- +5：组织或细胞来源明确
- -20：无法区分 RA/对照组
- -15：混合疾病队列中的 RA 子集不明确
- -15：平台或数据类型不明确

CSV 先按 `score` 降序排列，再按 `accession` 排列。

## 说明

GEO 元数据并非完全标准化。脚本采用保守策略：分组标签不明确时写入 `unknown`，而不猜测 RA/对照样本数，同时降低评分并添加人工复核提示。

仅应在人工整理时编辑 `manual_review`，并严格使用 `NULL`、`YES` 或 `NO`，不得留空。重新抓取已有 GSE 时，脚本会保留当前人工审核状态。

## 历史回填

使用 `--backfill-quarterly-local` 在本地运行从 `2016-01-01` 至今的完整回填。脚本每次处理一个季度，并在每个窗口完成后写入状态，因此中断后可以继续运行，而无需重复已完成的季度。

现有状态文件已经记录已完成的窗口（含 `2026-06-01` 至 `2026-06-30`），默认会跳过这些窗口。只有确实需要重新运行固定窗口时才使用 `--force-window`。

季度序列从以下窗口开始：

- `2016-01-01` 至 `2016-03-31`
- `2016-04-01` 至 `2016-06-30`
- `2016-07-01` 至 `2016-09-30`
- `2016-10-01` 至 `2016-12-31`

历史回填不使用 GitHub Actions。

## 月度更新

单仓库定时工作流为：

```text
.github/workflows/geo-ledgers-monthly.yml
```

该工作流位于仓库根目录，每月 1 日 01:00 UTC（北京时间 09:00）运行。各疾病通过并行矩阵任务更新并上传短期数据构件，最后合并为一次提交。若要手动仅运行 RA，请通过 `workflow_dispatch` 设置 `disease=RA`；也可以自定义 `since_days` 和 `retmax`。

对应的本地命令为：

```powershell
python .\scripts\build_ra_geo_ledger.py --start-date 2026-06-01 --end-date 2026-06-30 --retmax 500 --include-existing --state-file .\data\ra_geo_backfill_state.json --write-monthly-update --online-impact-factor-lookup --validate-csv --ncbi-email you@example.com --ncbi-api-key YOUR_KEY
```

月度更新写入：

- `data/ra_geo_transcriptome_datasets.csv`：长期累计总表
- `data/ra_geo_latest_monthly_update.csv`：最近一次月度更新表
- `data/monthly/ra_geo_monthly_update_YYYY-MM.csv`：按年月归档的更新表
- `data/online_impact_factor_cache.json`：可复用的在线影响因子缓存
- `data/monthly/impact_factor_audit_YYYY-MM.md`：影响因子查询审计报告

仓库当前最近一次月度归档为 `data/monthly/ra_geo_monthly_update_2026-06.csv`（检索窗口 `2026-06-01` 至 `2026-06-30`，命中 5 条），对应审计报告为 `data/monthly/impact_factor_audit_2026-06.md`。

月度更新表包含该次检索时间窗内所有通过过滤的记录，而不只是新发现的登录号。

可选的 GitHub Actions 密钥：

- `NCBI_EMAIL`
- `NCBI_API_KEY`
