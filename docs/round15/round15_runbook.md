# Round 15 复现 Runbook

> **目标**: 通过多 recycles + ESM2 + 稳定性 4 维度共识选 Top 6

---

## 一、ESM2 似然计算

```python
import torch
import numpy as np
from transformers import AutoModelForMaskedLM, AutoTokenizer

# 加载 ESM2 650M
model_name = "facebook/esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name).cuda()
model.eval()

def compute_esm2_log_prob(seq):
    """计算 ESM2 似然（log probability）"""
    with torch.no_grad():
        tokens = tokenizer(seq, return_tensors="pt")["input_ids"].cuda()
        outputs = model(tokens)
        logits = outputs.logits
        log_probs = torch.log_softmax(logits, dim=-1)
        # 累加每个 token 的对数概率
        token_log_probs = [log_probs[0, i, tokens[0, i]].item() for i in range(1, tokens.size(1)-1)]
        mean_log_prob = float(np.mean(token_log_probs))
    return mean_log_prob, np.exp(-mean_log_prob)  # log_prob, perplexity
```

## 二、多 recycles 投票

```python
import statistics
from transformers import AutoTokenizer, EsmForProteinFolding

# 加载 ESMFold
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()

def evaluate_multi_recycles(seq, recycles_list=[4, 6, 8, 12]):
    """多 recycles 投票"""
    tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
    scores = {}
    for r in recycles_list:
        with torch.no_grad():
            output = model(tokens, num_recycles=r)
        plddt = compute_plddt(output)
        ptm = float(output.ptm.cpu().item())
        chromo = plddt[58:73].mean()
        score = 0.40 * ptm + 0.30 * (plddt.mean() / 100) + 0.30 * (chromo / 100)
        scores[r] = {"plddt": plddt.mean(), "ptm": ptm, "chromo": chromo, "score": score}
    avg_score = sum(s["score"] for s in scores.values()) / len(scores)
    score_std = statistics.stdev([s["score"] for s in scores.values()])
    return {
        "scores_per_r": scores,
        "avg_score": avg_score,
        "score_std": score_std,
        "best_score": max(s["score"] for s in scores.values()),
        "best_recycles": max(scores, key=lambda r: scores[r]["score"])
    }
```

## 三、综合评分

```python
def compute_final_score(candidate):
    """4 维度综合评分"""
    return (
        0.40 * normalize(candidate["avg_score"]) +
        0.30 * normalize(candidate["best_score"]) +
        0.15 * normalize_esm2(candidate["esm2_log_prob"]) +
        0.15 * normalize_stability(candidate["score_std"])
    )
```

## 四、关键参数

| 参数 | 推荐值 | 备注 |
|:-----|:-------|:-----|
| recycles_list | [4, 6, 8, 12] | 多 recycles 投票 |
| ESM2 模型 | facebook/esm2_t33_650M_UR50D | 650M 参数 |
| 权重 | 0.4/0.3/0.15/0.15 | 多模型共识 |
| final_score 阈值 | 0.92+ | Top 6 |

---

*Runbook 作者: Trae AI Agent (Claude)*
*最后更新: 2026-06-25*
