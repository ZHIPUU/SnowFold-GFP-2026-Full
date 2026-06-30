"""Step 4: ESMFold 验证所有候选序列(pLDDT),使用 HuggingFace transformers。"""
import time
import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
CAND_FILE = ROOT / "candidates_round2_design.csv"
OUT_JSON = ROOT / "step4_esmfold_results.json"
OUT_CSV = ROOT / "step4_esmfold_results.csv"

# 加载 ESMFold (HuggingFace facebook/esmfold_v1)
from transformers import EsmForProteinFolding, AutoTokenizer
print("加载 ESMFold 模型 (facebook/esmfold_v1, ~14B params)...")
t0 = time.time()
try:
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1",
        torch_dtype=torch.float16,  # FP16 省 VRAM
        low_cpu_mem_usage=True,
    ).cuda().eval()
    print(f"  加载完成, 耗时 {time.time()-t0:.1f}s, GPU mem={torch.cuda.memory_allocated()/1e9:.2f} GB")
except Exception as e:
    print(f"❌ ESMFold 加载失败: {e}")
    raise

# 读候选
cand = pd.read_csv(CAND_FILE)
print(f"\n候选 {len(cand)} 条\n")

results = []
for i, row in cand.iterrows():
    name = row["name"]
    seq = row["seq"]
    print(f"[{i+1}/{len(cand)}] ESMFold: {name} (len={len(seq)})", flush=True)
    t0 = time.time()
    try:
        # Tokenize
        ids = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        ids = {k: v.cuda() for k, v in ids.items()}
        # Infer with attention_mask to suppress warning
        with torch.no_grad():
            output = model(
                input_ids=ids["input_ids"],
                attention_mask=ids.get("attention_mask", torch.ones_like(ids["input_ids"])),
                num_recycles=1,  # 1 recycle for speed
            )
        # output: PlDDT per residue, positions
        plddt = output.plddt[0].cpu().float().numpy()  # (L,)
        plddt_mean = float(plddt.mean())
        plddt_min = float(plddt.min())
        plddt_max = float(plldt.max())
        # positions: (L, 14, 3) atom14 coords (last frame)
        pos = output.positions[-1, 0].cpu().float().numpy()  # (L, 14, 3)
        ca = pos[:, 1, :]  # CA atoms
        # Rg
        rg = float(np.sqrt(((ca - ca.mean(axis=0))**2).sum(axis=1).mean()))
        # CA-CA min distance
        from scipy.spatial.distance import pdist
        dists = pdist(ca)
        min_dist = float(dists.min())
        has_clash = min_dist < 3.0
        elapsed = time.time() - t0
        print(f"  pLDDT mean={plddt_mean:.2f} min={plddt_min:.2f}, Rg={rg:.1f}Å, CA-CA min={min_dist:.1f}Å, 耗时 {elapsed:.1f}s", flush=True)
        if has_clash:
            print(f"  ⚠️ CA-CA 冲突 min={min_dist:.1f}Å")
        if plddt_mean < 60:
            print(f"  ⚠️ pLDDT 低 (<60), 可能未正确折叠")
        results.append({
            "id": int(row["id"]),
            "name": name,
            "scaffold": row["scaffold"],
            "length": len(seq),
            "plddt_mean": plddt_mean,
            "plddt_min": plddt_min,
            "plddt_max": plddt_max,
            "rg": rg,
            "ca_min_dist": min_dist,
            "has_clash": has_clash,
            "fold_time": elapsed,
            "status": "OK",
        })
        torch.cuda.empty_cache()
    except torch.cuda.OutOfMemoryError as e:
        print(f"  ❌ OOM: {e}")
        torch.cuda.empty_cache()
        results.append({
            "id": int(row["id"]),
            "name": name,
            "scaffold": row["scaffold"],
            "length": len(seq),
            "plddt_mean": float("nan"),
            "rg": float("nan"),
            "ca_min_dist": float("nan"),
            "has_clash": False,
            "fold_time": time.time() - t0,
            "status": "OOM",
        })
    except Exception as e:
        print(f"  ❌ 错误: {e}")
        torch.cuda.empty_cache()
        results.append({
            "id": int(row["id"]),
            "name": name,
            "scaffold": row["scaffold"],
            "length": len(seq),
            "plddt_mean": float("nan"),
            "rg": float("nan"),
            "ca_min_dist": float("nan"),
            "has_clash": False,
            "fold_time": time.time() - t0,
            "status": f"ERROR: {str(e)[:100]}",
        })

# 保存
df_res = pd.DataFrame(results).sort_values("plddt_mean", ascending=False)
df_res.to_csv(OUT_CSV, index=False)
with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print("\n=== ESMFold 全部完成 ===")
print(df_res[["id", "name", "plddt_mean", "rg", "ca_min_dist", "has_clash", "status"]].to_string(index=False))