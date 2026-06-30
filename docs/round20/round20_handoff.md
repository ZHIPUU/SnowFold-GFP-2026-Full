# Round 20 接手指南

## 当前状态

- ✅ **Phase 1 完成**: R20 Top 1 = 0.9396 (项目新高!)
- 🔄 **Phase 2 进行中**: R21 大规模 MPNN 探索 (P3/6 200/600)
- ⏳ **Phase 3 待启动**: Top 50 + r=20 高精度重算

## R20 Top 6 (✅ ALL PASS)

| Seq | Score | pTM | pLDDT | Chromo |
|:---:|:-----:|:---:|:-----:|:------:|
| 1 | **0.9396** | 0.9250 | 0.934 | 0.965 |
| 2 | 0.9391 | 0.9279 | 0.932 | 0.961 |
| 3 | 0.9391 | 0.9277 | 0.933 | 0.961 |
| 4 | 0.9389 | 0.9260 | 0.934 | 0.961 |
| 5 | 0.9387 | 0.9277 | 0.935 | 0.957 |
| 6 | 0.9382 | 0.9251 | 0.934 | 0.960 |

## 提交 CSV

`D:\workspace\round20\submission_r20.csv` (6 候选, 全部合规)

## 关键修复 (R20 vs R18/R19)

```python
# R18/R19 ❌ Buggy
FIXED = [65, 66, 67, 96, 222]

# R20 ✅ Fixed
FIXED = [1, 65, 66, 67, 96, 222]  # ← 关键: 加 1 (M) 到 fixed
```

通过率从 ~8% 提升到 62.8%。

## 下一步

R22 任务还在跑 (`a0885c7b13f8`)，Phase 2/3 完成后出最终 R22 Top 6。

## 复现命令

```bash
# 上传
gssh cp "D:\workspace\r22_long.py" 9ca7acb1b94c:/root/autodl-tmp/r22_long.py

# 启动
gssh run 9ca7acb1b94c "cd /root/autodl-tmp && stdbuf -oL -eL python3 -u r22_long.py > /tmp/r22.log 2>&1"

# 监控
gssh --json exec 9ca7acb1b94c "tail -30 /tmp/r22.log"
```

## 关键文件路径

- 远程: `/root/autodl-tmp/r22/results/r20_top6.json` — Phase 1 输出
- 远程: `/root/autodl-tmp/r22/results/phase1_progress.json` — 完整 1227 通过候选
- 本地: `D:\workspace\round20\` — 下载的所有 R20 数据
