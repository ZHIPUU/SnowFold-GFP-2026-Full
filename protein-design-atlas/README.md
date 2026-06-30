# Protein Design Atlas

本项目把 R4-R26 蛋白设计轮次中的 JSON/CSV/FA/PDB/Markdown/Pipeline 脚本汇总为一个本地 SQLite 数据库，并提供三栏式 API 文档风格可视化网站。

## 当前导入规模

本地已成功导入：

- 238,691 条去重序列
- 402,650 条指标记录
- 5,120 个文件证据 artifact
- 25 个轮次

数据库位置：`data/design.db`

## 技术栈

- FastAPI + SQLite + pandas/networkx
- React + Vite + Cytoscape.js + Plotly.js
- R/ggplot2/ggraph/patchwork 图表 worker
- Docker Compose 部署

## 端口

默认使用高位端口，避免冲突：

- 前端：`http://localhost:18082`
- 后端 API：`http://localhost:18000`
- Swagger：`http://localhost:18000/docs`

当前本地开发服务已经通过端到端测试：

- `GET /api/stats` 返回 238,691 条序列与当前最佳 R26=0.9449
- `GET /api/metrics/top?limit=2` 返回 Top 指标
- 前端 `http://127.0.0.1:18082/` 返回 HTTP 200
- Swagger `http://127.0.0.1:18000/docs` 返回 HTTP 200

如果冲突，复制 `.env.example` 为 `.env` 后修改端口。

## 本地非 Docker 快速启动

```powershell
cd D:\workspace\protein-design-atlas\backend
$env:PYTHONPATH='D:\workspace\protein-design-atlas\backend'
$env:ATLAS_DB='D:\workspace\protein-design-atlas\data\design.db'
$env:ATLAS_SOURCE_ROOTS='D:/workspace;D:/生信/2026Protein Design'
python -m app.services.importer
uvicorn app.main:app --host 0.0.0.0 --port 18000
```

前端：

```powershell
cd D:\workspace\protein-design-atlas\frontend
npm install
$env:VITE_API_URL='http://localhost:18000'
npm run dev -- --port 18082
```

## Docker 启动

```powershell
cd D:\workspace\protein-design-atlas
copy .env.example .env
docker compose up --build
```

生成 R 图表：

```powershell
docker compose --profile plots run --rm r-worker
```

## 主要页面

- 总览 Dashboard：轮次趋势、pTM × chromophore 图、核心计数
- 序列库：按分数/轮次筛选所有序列
- 拓扑网络：Cytoscape.js 谱系图
- 轮次文档：直接浏览 docs/round* markdown
- 图表：Plotly 实时图 + R worker 产物索引

## API 示例

```text
GET /api/stats
GET /api/rounds
GET /api/sequences?min_score=0.94&limit=100
GET /api/metrics/top?limit=50
GET /api/graph/lineage?min_score=0.94
GET /api/documents?round_key=R22
```

## 数据模型

核心表：

- `rounds`
- `sequences`
- `metrics`
- `lineage_edges`
- `artifacts`
- `submissions`
- `documents`

序列使用 `sha256(sequence)[:16]` 去重。

## 设计风格

界面采用“黑曜石实验日志 + 金色航迹”的三栏布局：

- 左栏：项目导航 / 轮次
- 中栏：数据、图表、文档、网络
- 右栏：Inspector / API 示例 / 当前最佳

## 后续可扩展

- 增加 R 图表 artifact API 索引
- 增加多样性约束 Top6 选择器
- 增加 mutation hotspot network
- 增加高精度 r=20 对比视图
- 接入远程服务器任务日志
