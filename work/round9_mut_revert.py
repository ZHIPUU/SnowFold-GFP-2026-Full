"""
Round 9: 从 MPNN_T01_014 回退突变构建低突变变体
================================================
策略: 对 57-mut 序列，回退 27 个"最不重要"的突变回 WT
保留 30 个"最重要"的突变，测试 pTM 保持情况
"""
import json, time, torch, numpy as np, random
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path(r"D:\生信\2026Protein Design")
R9 = ROOT / "work" / "round9"
R9.mkdir(parents=True, exist_ok=True)

WT = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# 加载 MPNN_T01_014 序列
with open(R9.parent / "round6" / "final_6_new_rules.json", encoding="utf-8") as f:
    data = json.load(f)
seq1 = [d["seq"] for d in data if "MPNN_T01_014" in d["name"]][0]

# 找出 57 个突变位点
mut_positions = []
mut_from = {}  # pos -> WT aa
mut_to = {}    # pos -> MPNN aa
for i, (a, b) in enumerate(zip(WT, seq1)):
    if a != b:
        mut_positions.append(i)
        mut_from[i] = a
        mut_to[i] = b

print(f"Total mutations: {len(mut_positions)}")

# 分类: 生色团/核心区 = "重要" (保留), 其他 = "不重要的" (回退)
# 生色团: Y66(65), G67(66), R96(95), E222(221)
# 核心 β-barrel: 58-72, 210-230 区域
chromo = set(range(57, 72)) | set(range(209, 230))  # 0-based

important = [p for p in mut_positions if p in chromo]
unimportant = [p for p in mut_positions if p not in chromo]

print(f"Important (keep): {len(important)}")
print(f"Unimportant (can revert): {len(unimportant)}")

# 目标: 保留 ~30 个突变
# 自动: 保留所有 important (8个), 再从 unimportant 中选 22 个保留
# 需要回退 27 个
n_keep_important = len(important)  # 8
n_keep_from_unimportant = 30 - n_keep_important  # 22
n_revert = len(unimportant) - n_keep_from_unimportant  # 27

print(f"\n保留: {n_keep_important} 重要 + {n_keep_from_unimportant} 不重要 = {n_keep_important + n_keep_from_unimportant}")
print(f"回退: {n_revert}")

# ============================================================
# 构建变体: 挑选不同的不重要位点子集回退
# ============================================================
def build_variant(keep_positions):
    """构建变体: 在 WT 序列上只应用 keep_positions 位置的突变"""
    seq = list(WT)
    for p in keep_positions:
        seq[p] = mut_to[p]
    return "".join(seq)

# 排序不重要位点 (按位置, 但也可以随机化)
unimportant_sorted = sorted(unimportant)

# 生成 50 个不同变体, 每次随机挑选保留的位点
random.seed(42)
variants = []
seen = set()

for _ in range(200):  # 生成多一点再筛选
    # 从不重要位点中随机选 n_keep_from_unimportant 个保留
    kept = random.sample(unimportant_sorted, n_keep_from_unimportant)
    keep_set = set(important) | set(kept)
    seq = build_variant(keep_set)
    
    if seq not in seen:
        seen.add(seq)
        n_muts = sum(1 for a, b in zip(seq, WT) if a != b)
        variants.append({
            "name": f"R9_rm_{len(variants)+1:03d}",
            "n_muts": n_muts,
            "seq": seq,
            "n_keep": len(keep_set),
        })
    
    if len(variants) >= 50:
        break

# 还保留了一些在 chromo 区外的不重要突变
# 但有些位置可能本身就对结构重要, 让我们也生成一些"激进回退"的变体
for n_extra_revert in [5, 10, 15]:
    n_keep_extra = max(15, n_keep_from_unimportant - n_extra_revert)
    for _ in range(20):
        kept = random.sample(unimportant_sorted, min(n_keep_extra, len(unimportant_sorted)))
        keep_set = set(important) | set(kept)
        seq = build_variant(keep_set)
        if seq not in seen:
            seen.add(seq)
            n_muts = sum(1 for a, b in zip(seq, WT) if a != b)
            variants.append({
                "name": f"R9_rm_{len(variants)+1:03d}",
                "n_muts": n_muts,
                "seq": seq,
                "n_keep": len(keep_set),
            })

print(f"\nVariants: {len(variants)}")
nmuts_list = [v["n_muts"] for v in variants]
print(f"Muts range: {min(nmuts_list)}-{max(nmuts_list)}")
for t in [20, 25, 30, 35]:
    print(f"  <{t}: {sum(1 for m in nmuts_list if m < t)}")

# 保存
with open(R9 / "variants_r9.json", "w", encoding="utf-8") as f:
    json.dump(variants, f, indent=2, ensure_ascii=False)

# ============================================================
# ESMFold 验证
# ============================================================
print(f"\n{'='*60}")
print(f"ESMFold validating {len(variants)} variants...")

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True
).cuda().eval()

t0 = time.time()
for i, v in enumerate(variants):
    with torch.no_grad():
        tok = tokenizer([v["seq"]], return_tensors="pt", add_special_tokens=False)
        tok = {k: v.cuda() for k, v in tok.items()}
        out = model(**tok)
    
    plddt = out["plddt"][0, 1:len(v["seq"])+1].cpu().numpy()
    plddt_r = plddt.mean(axis=1) * 100
    ptm = out["ptm"].item()
    r1, r2 = slice(57, 72), slice(209, 230)
    chromo_val = np.concatenate([plddt_r[r1], plddt_r[r2]]).mean() if len(plddt_r) > 230 else float(plddt_r.mean())
    
    v["plddt_mean"] = float(plddt_r.mean())
    v["plddt_chromo"] = float(chromo_val)
    v["ptm"] = float(ptm)
    v["score"] = round(ptm * 0.50 + plddt_r.mean() / 100 * 0.30 + chromo_val / 100 * 0.20 - v["n_muts"] / 200 * 0.05, 4)
    
    if (i+1) % 20 == 0:
        print(f"  [{i+1}/{len(variants)}] {time.time()-t0:.0f}s")

print(f"Done! {time.time()-t0:.0f}s")

# ============================================================
# 结果分析
# ============================================================
variants.sort(key=lambda x: -x["ptm"])

ref_ptm = 0.765  # MPNN_T01_014 pTM
threshold = ref_ptm - 0.05  # 0.715

good = [v for v in variants if v["ptm"] > threshold]

print(f"\n{'='*60}")
print(f"Ref pTM={ref_ptm}, threshold={threshold}")
print(f"pTM保持: {len(good)}/{len(variants)}")
print(f"\nTop-10 by pTM:")
print(f"{'#':<4} {'Name':<15} {'Muts':>4} {'pLDDT':>6} {'pTM':>6} {'Score':>7} {'Keep':>5}")
for i, v in enumerate(variants[:10], 1):
    flag = "✅" if v["ptm"] > threshold else ""
    print(f"{i:<4} {v['name']:<15} {v['n_muts']:>4} {v['plddt_mean']:>6.1f} {v['ptm']:>6.4f} {v['score']:>7.4f} {v['n_keep']:>5} {flag}")

# 最佳折中: 突变数少 + pTM 高
print(f"\n折中方案 (muts<30, sort by pTM):")
tradeoff = [v for v in variants if v["n_muts"] < 30]
tradeoff.sort(key=lambda x: -x["ptm"])
for v in tradeoff[:5]:
    print(f"  {v['name']}: muts={v['n_muts']}, pLDDT={v['plddt_mean']:.1f}, pTM={v['ptm']:.4f}")

# 保存
with open(R9 / "esmfold_r9.json", "w", encoding="utf-8") as f:
    json.dump(variants, f, indent=2, ensure_ascii=False)
print(f"\nSaved: {R9 / 'esmfold_r9.json'}")
