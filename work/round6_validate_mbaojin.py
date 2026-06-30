"""
验证 mBaoJin: ESMFold 是否系统性低估非-avGFP 骨架?

方法:
1. 用 ESMFold 折叠 mBaoJin WT (已知晶体结构 PDB 8QBJ)
2. 折叠 sfGFP WT (比赛 WT, 已知高 pLDDT) 作为对照
3. 折叠 mBaoJin 候选 M5
4. 对比 pLDDT 分布

如果 mBaoJin WT 也低 pLDDT → 模型偏差 → 信任文献
如果 mBaoJin WT 高但候选低 → 候选真的有问题
"""
import json, time, torch, numpy as np
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design")

# === mBaoJin WT (回退 M5_D230E: E230→D) ===
MBAOJIN_WT = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"
SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# === 加载 M5_mBaoJin_D230E 序列 ===
with open(ROOT / "work" / "round5" / "evolvepro_scored.json", encoding="utf-8") as f:
    all_c = json.load(f)
m5_seq = None
for c in all_c:
    if c.get("name") == "M5_mBaoJin_D230E":
        m5_seq = c["seq"]
        print(f"找到 M5: len={len(m5_seq)}")
        break

print(f"\nmBaoJin WT len={len(MBAOJIN_WT)}")
print(f"sfGFP len={len(SFGFP)}")
print(f"\n序列差异: MBAOJIN vs sfGFP")
for i, (a, b) in enumerate(zip(MBAOJIN_WT[:50], SFGFP[:50])):
    if a != b:
        print(f"  pos {i+1}: {a} vs {b}")

# === 加载 ESMFold ===
print("\n加载 ESMFold...")
import esm
esmfold = esm.pretrained.esmfold_v1().cuda().eval()
esmfold.set_chunk_size(64)  # 节省显存
print("ESMFold 已加载")

def fold_and_score(name, seq):
    print(f"\n--- 折叠 {name} (len={len(seq)}) ---")
    t0 = time.time()
    with torch.no_grad():
        output = esmfold.infer(seq)
    elapsed = time.time() - t0
    
    # pLDDT: 每个残基的预测局部置信度
    plddt = output["plddt"][0, 1:len(seq)+1].cpu().numpy()
    ptm = output["ptm"].item()
    mean_plddt = plddt.mean()
    
    # Chromophore 区域 (GFP 的 TYG/CRO 三联体在 ~65-67 位)
    # sfGFP chromophore: Y66-G67
    # mBaoJin: 从序列看 chromophore 可能在类似位置
    # 简单用中心 40% 残基估算
    mid_start = len(seq) // 3
    mid_end = 2 * len(seq) // 3
    chromo_plddt = plddt[mid_start:mid_end].mean()
    
    # 残基级 pLDDT 分布
    lt50 = (plddt < 50).sum()
    gt70 = (plddt > 70).sum()
    
    print(f"  耗时: {elapsed:.1f}s")
    print(f"  平均 pLDDT: {mean_plddt:.2f}")
    print(f"  中部 pLDDT: {chromo_plddt:.2f}")
    print(f"  pTM: {ptm:.4f}")
    print(f"  pLDDT < 50: {lt50}/{len(seq)} ({100*lt50/len(seq):.0f}%)")
    print(f"  pLDDT > 70: {gt70}/{len(seq)} ({100*gt70/len(seq):.0f}%)")
    
    return {
        "name": name, "length": len(seq),
        "plddt_mean": float(mean_plddt),
        "plddt_chromo_region": float(chromo_plddt),
        "plddt_lt50": int(lt50),
        "plddt_gt70": int(gt70),
        "ptm": float(ptm),
        "elapsed_s": float(elapsed),
    }

results = []
results.append(fold_and_score("sfGFP_WT", SFGFP))
results.append(fold_and_score("mBaoJin_WT", MBAOJIN_WT))
results.append(fold_and_score("M5_mBaoJin_D230E", m5_seq))

# === 对比分析 ===
print(f"\n{'='*80}")
print(f"对比分析: ESMFold 是否对 mBaoJin 有偏见?")
print(f"{'='*80}")
print(f"{'Name':<20} {'pLDDT':>7} {'pTM':>7} {'<50':>5} {'>70':>5} {'time':>6}")
print("-" * 60)
for r in results:
    print(f"{r['name']:<20} {r['plddt_mean']:>7.2f} {r['ptm']:>7.4f} "
          f"{r['plddt_lt50']:>5} {r['plddt_gt70']:>5} {r['elapsed_s']:>6.1f}s")

print(f"\n结论:")
sf_plddt = results[0]["plddt_mean"]
mb_plddt = results[1]["plddt_mean"]
m5_plddt = results[2]["plddt_mean"]
print(f"  sfGFP WT pLDDT = {sf_plddt:.1f}")
print(f"  mBaoJin WT pLDDT = {mb_plddt:.1f}")
print(f"  M5 candidate pLDDT = {m5_plddt:.1f}")

if mb_plddt < sf_plddt - 10:
    print(f"  >>> ESMFold 对 mBaoJin 骨架系统性低估（差 {sf_plddt-mb_plddt:.1f} 分）")
    print(f"  >>> 建议: 信任文献 Tm=92°C, 减少 pLDDT 罚分")
elif m5_plddt < mb_plddt - 5:
    print(f"  >>> mBaoJin 候选（M5）比 WT 差很多")
    print(f"  >>> 建议: 保持 pLDDT 罚分，使用最佳 mBaoJin 候选")
else:
    print(f"  >>> mBaoJin WT 和候选 pLDDT 相近")
    print(f"  >>> 综合判断中...")

# 保存
output_path = ROOT / "work" / "round6" / "mbaojin_validation.json"
output_path.parent.mkdir(parents=True, exist_ok=True)
with open(output_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\n结果已保存: {output_path}")
