# Round 16 接手指南

## 当前状态
- RTX 5080 16GB 已达到物理上限
- 已发现 R15 综合分 0.98 实际是 fair-esm 评分 bug
- 已切到云端 A800 80GB

## 留给接手者的关键事
1. **永远用 HuggingFace ESMFold 评分**(`out.plddt[..., 1]` 是 0-1 尺度)
2. **不要混用 fair-esm 和 transformers** 的 pLDDT 输出
3. **R17 已修正评分**,请优先以 R17 真实分数为准

## 接手快速命令
```bash
# 重现 R17 重评
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r17_pipeline.py"
```
