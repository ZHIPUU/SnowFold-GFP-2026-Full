# Protein Design Atlas 数据字典

## 数据库

位置：`data/design.db`

当前导入规模：

- sequences: 238,691
- metrics: 402,650
- artifacts: 5,120
- rounds: 25
- documents: 157
- submissions: 198
- lineage_edges: 4,616

## 表

### rounds
轮次元数据。`round_key` 如 R22，`best_score` 为该轮最高有效 0-1 分数。

### sequences
按 `sha256(sequence)[:16]` 去重的序列主表。保留长度、M 开头、合法 AA、最佳分数、最佳 pTM/pLDDT/chromo、来源轮次。

### metrics
每个 artifact 中解析出的评分记录。同一序列可能出现多条 metrics，因为它可能同时出现在 progress、final、submission 或 consensus 文件中。

### lineage_edges
由 JSON/FASTA 中的 `parent` 字段推断的父代标签边。部分 `parent_sequence_id` 无法精确匹配，保留在 `parent_label` 中。

### artifacts
所有被扫描到的 JSON/CSV/FA/PDB/JSONL/MD/PY/LOG 文件证据。

### submissions
所有 submission CSV 中的提交序列记录。

### documents
所有 Markdown 文档正文，可在网站 Docs 页面浏览。

## 评分尺度注意

早期轮次存在分数 >1 的历史估计值；API 当前最佳与 Top metrics 默认过滤 `score BETWEEN 0 AND 1`，避免历史非同尺度估计污染榜单。

## 当前最佳

- R26 Top1 = 0.9449
- R24 当前 = 0.9447
- R22 = 0.9430

## API

- `/api/stats`
- `/api/rounds`
- `/api/sequences`
- `/api/metrics/top`
- `/api/graph/lineage`
- `/api/documents`
- `/api/plots/round-trend`
- `/api/plots/score-scatter`
