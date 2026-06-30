# Round 22 接手指南

## 当前状态

- ✅ **R22 完成**: sort_score **0.9430** 🏆 项目新纪录
- ✅ **6/6 ALL PASS** 合规
- 🔄 **R23 在服务器上跑** (~4.5 小时, ~18:30 完成)

## R22 Top 6 (推荐提交)

| Seq | Score | pTM | pLDDT | Chromo |
|:---:|:-----:|:---:|:-----:|:------:|
| **1** | **0.9430** | 0.9276 | 0.943 | 0.964 |
| 2 | 0.9429 | 0.9300 | 0.938 | 0.965 |
| 3 | 0.9416 | 0.9295 | 0.939 | 0.961 |
| 4 | 0.9416 | 0.9282 | 0.943 | 0.958 |
| 5 | 0.9414 | 0.9284 | 0.939 | 0.961 |
| 6 | 0.9413 | 0.9281 | 0.938 | 0.962 |

## 提交 CSV

📄 `D:\workspace\round22\submission_r22.csv` (6 候选, 全部合规)

## 关键修复（同 R20）

```python
# 关键: fixed position 1 (M)
FIXED = [1, 65, 66, 67, 96, 222]  # ← 关键
```

R18/R19 之前缺 position 1 导致 97% 候选不以 M 开头。

## 与历届最佳对比

| 轮次 | Top 1 | Δ vs R19 |
|:----:|:-----:|:-------:|
| R19 | 0.9321 | 基准 |
| R20 | 0.9396 | +0.80% |
| **R22** | **0.9430** | **+1.17%** 🏆 |

## 下一步

R23 在服务器上跑（ID `0e6f316839ca`），用 R20 Top 3 父代 + 高温度找更多样性。
- 预计完成 ~18:30
- 完成时下载 + 对比 R22

## 复现命令

```bash
gssh cp "D:\workspace\r22_long.py" 9ca7acb1b94c:/root/autodl-tmp/r22_long.py
gssh cp "C:\proteinmpnn_r10\*" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/

gssh run 9ca7acb1b94c "cd /root/autodl-tmp && stdbuf -oL python3 -u r22_long.py > /tmp/r22.log 2>&1"
```

## 监控命令

```bash
gssh --json logs a0885c7b13f8 --tail 50 # R22 已完成 (历史)
gssh --json logs 0e6f316839ca --tail 50 # R23 (当前)
```
