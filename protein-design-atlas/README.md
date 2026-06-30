# Protein Design Atlas

> SnowFold GFP Design 2026 的实验数据全景数据库与交互式可视化平台。

## 在线访问

部署后的交互式网站可通过以下地址访问：

- **交互式网站**: http://120.48.98.164:18082
- **API 文档 (Swagger)**: http://120.48.98.164:18000/docs

## 项目简介

将 27 轮蛋白设计迭代产生的全部序列、指标、文档整合为 SQLite 数据库，并提供三栏式 API 文档风格可视化网站，支持序列检索、谱系追踪、网络拓扑、图表联动。

## 数据库规模

| 表 | 记录数 | 内容 |
|:---|------:|:-----|
| sequences | 243,386 | 去重氨基酸序列（SHA-256 主键） |
| metrics | 1,474,887 | 全部预测指标（多轮次、多 recycle） |
| rounds | 27 | 各轮次元数据 |
| artifacts | 5,539 | 文件证据索引 |
| lineage_edges | 16,000+ | 序列谱系关系 |
| documents | 157 | 实验文档全文 |

## 技术栈

- **后端**: FastAPI + SQLite
- **前端**: React + Vite + Cytoscape.js + Plotly.js
- **图表**: R/ggplot2 + patchwork
- **部署**: Docker Compose

## 功能模块

| 模块 | 说明 |
|:-----|:-----|
| Dashboard | 分数趋势图、pTM × chromo 散点图、指标卡片 |
| Sequence Vault | 24 万条序列检索，按轮次/分数/长度筛选 |
| Topology Network | Cytoscape 谱系网络图，节点颜色映射分数 |
| Docs Browser | 157 篇实验文档在线浏览 |
| R 图表报告 | 自动生成 6 张高质量分析图表 |

## 关联仓库

本项目是 SnowFold GFP Design 2026 的数据可视化子系统：

| 仓库 | 说明 | 链接 |
|:-----|:-----|:-----|
| **竞赛核心仓库** | 标准化 pipeline + 最终结果 | [SnowFold-GFP-2026](https://github.com/ZHIPUU/SnowFold-GFP-2026) |
| **完整实验记录** | R2-R27 全部脚本/数据/文档 | [SnowFold-GFP-2026-Full](https://github.com/ZHIPUU/SnowFold-GFP-2026-Full) |
| **gssh CLI** | 远程 GPU 服务器管理工具 | [gssh](https://github.com/ZHIPUU/gssh) |
| **Protein Design Atlas** (本仓库) | 数据库 + 可视化平台 | [protein-design-atlas](https://github.com/ZHIPUU/protein-design-atlas) |

## Docker 部署

```bash
# 构建并启动
docker compose up --build

# 访问
# 前端: http://localhost:18082
# API:  http://localhost:18000
# Swagger: http://localhost:18000/docs
```

## License

MIT
