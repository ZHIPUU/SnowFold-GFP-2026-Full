# Round 17 接手指南

## R17 状态
- ✅ A800 80GB 已可用
- ✅ HF ESMFold 评分校准
- ✅ R15 Top 6 真实分数 0.908

## R17 Top 6 (最终提交)
| Seq | Name | Score |
|:---:|:----|:-----:|
| 1 | R14_A_T01_013 | 0.9080 |
| 2 | R14_A_T02_037 | 0.9075 |
| 3 | R14_A_T01_020 | 0.8961 |
| 4 | R14_D_T02_033 | 0.8918 |
| 5 | R14_A_T01_023 | 0.8896 |
| 6 | R14_D_T02_039 | 0.8841 |

## 留给接手者
1. **永远用 HF `out.plddt[..., 1]` 评分**(0-1 尺度)
2. **R17 校准后, R14_A_T01_013 是真 Top 1**
3. **提交 CSV**: `D:\workspace\round17\submission_r17.csv`

## 快速复现
```bash
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r17_pipeline.py"
```
