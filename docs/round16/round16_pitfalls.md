# Round 16 踩坑记录

> **核心教训**: **R15 的 "综合分 0.98" 实际是 fair-esm 0-100 尺度评分 bug**

---

## 坑 1: fair-esm vs HuggingFace pLDDT 尺度不一致 ⚠️ 严重

**症状**: 同一序列用两套 ESMFold 实现评分,**分数相差 100 倍**

**根因**:
- `fair-esm` 的 `out["plddt"]` 输出 0-100 整数(实际是中心化的预测区间)
- `transformers` 的 `out.plddt[..., 1]` 是 0-1 浮点(取 argmax 索引对应的 bin 中心)

**bug 触发条件**: R15 评分公式 `0.3 * pLDDT/100`,但实际加载的是 fair-esm 输出

**修复方法**:
```python
# fair-esm
plogits = out["plddt"]                # (L, 37)
probs = F.softmax(plogits, dim=-1)
centers = torch.linspace(0.5/37, 1-0.5/37, 37, device=logits.device)
plddt_01 = (probs * centers).sum(-1)  # 0-1 尺度
# 或直接用整数 / 100

# HuggingFace
plddt_01 = out.plddt[0, :, 1]        # 已经是 0-1 尺度
```

**教训**: 跨实现迁移评分函数时,**必须验证输出数值范围**。

---

## 坑 2: RTX 5080 16GB OOM 边界

**症状**: `num_recycles=16` 评估时随机 OOM

**根因**:
- ESMFold 每步 recycle 都会扩大中间激活
- 16GB VRAM 在 r=12 时接近上限
- r=16 + ESM2-3B 同时加载 → 必然 OOM

**解决方案**:
```python
model.trunk.set_chunk_size(64)  # 减半 chunk
with torch.cuda.amp.autocast():  # FP16
    out = model(tokens, num_recycles=16)
```

**教训**: 在 16GB 设备上,**chunk_size 和 recycles 必须协同调整**。

---

## 坑 3: MPNN 父代 PDB 必须用 ESMFold 预测而非晶体

**症状**: R11 用 2B3P 晶体(228aa) → MPNN 强制截短 → 序列变短

**根因**: 2B3P 晶体缺失两端残基(2 个 loop 区无电子密度)

**解决方案**: **永远用 ESMFold 预测的完整 PDB 作为 MPNN 输入**

---

## 坑 4: 网络 turbo 模式导致 git clone 超时

**症状**: `git clone https://github.com/...` 失败(timeout)

**根因**: `source /etc/network_turbo` 关闭后,github 直连超时

**解决方案**:
```bash
source /etc/network_turbo  # 加速 pip 但 git 不加速
# 或用 ghproxy.com
git clone https://ghproxy.com/https://github.com/...
```

---

## 坑 5: `gssh cp` 路径格式严格

**症状**: `gssh cp D:/file.txt 9ca7acb1b94c:/path/` 报"一端必须是 <session-id>:<path>"

**根因**: gssh 的 cp 解析器对 `D:/...` 形式的 Windows 路径在 PowerShell 下解析有歧义

**解决方案**:
```bash
gssh cp "D:\workspace\file.py" session-id:/remote/path/  # ✅ 双引号
```

---

*记录者: Trae AI Agent (Claude) | 最后更新: 2026-06-28*
