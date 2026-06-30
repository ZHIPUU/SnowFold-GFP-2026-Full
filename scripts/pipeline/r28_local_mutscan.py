#!/usr/bin/env python3
"""R28 Local Mutation/Recombination Breadth Scan (~1h)

New method vs R27:
- No MPNN. Generate sequence variants directly from top candidates.
- Breadth directions:
  A: single conservative mutations
  B: double conservative mutations
  C: top-parent segment crossovers
  D: consensus-frequency sampling from top pool
  E: local back/forward swaps among R25/R26/R27 leaders

Runs after R27 finishes, uses local RTX 5080 with ESMFold r=8 chunk=64.
"""
import os, json, random, time, warnings, csv
from pathlib import Path
import torch
import numpy as np
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = Path("D:/workspace/round28_local")
WORK.mkdir(parents=True, exist_ok=True)
(WORK/"results").mkdir(exist_ok=True)
RECYCLES = 8
TARGET_N = 140
FIXED = {1, 65, 66, 67, 96, 222}
AA = "ACDEFGHIKLMNPQRSTVWY"
random.seed(2801)

GROUPS = [
    "AVLIMFWY",   # hydrophobic/aromatic-ish
    "STNQCY",     # polar
    "KRH",        # positive
    "DE",         # negative
    "GP",         # turn/small-special
]
AA_GROUP = {}
for g in GROUPS:
    for a in g:
        AA_GROUP[a] = g

def load_json(path):
    if Path(path).exists():
        try:
            return json.load(open(path, encoding="utf-8"))
        except Exception:
            return []
    return []

def best_pool():
    pool = []
    for p in [
        r"D:\workspace\round27_diverge\final_6_r27.json",
        r"D:\workspace\round27_diverge\all_passed.json",
        r"D:\workspace\round25\final_6_r25.json",
        r"D:\workspace\round26_local\final_6_r26.json",
        r"D:\workspace\round24\final_6_r24.json",
        r"D:\workspace\round22\final_6_r22.json",
    ]:
        data = load_json(p)
        if isinstance(data, list):
            pool.extend(data[:20])
    # Dedup by seq, sort by score
    seen = {}
    for x in pool:
        s = x.get("seq") or x.get("Sequence") or x.get("sequence")
        if not s: continue
        y = dict(x); y["seq"] = s
        if s not in seen or y.get("score", 0) > seen[s].get("score", 0):
            seen[s] = y
    pool = sorted(seen.values(), key=lambda x: x.get("score", 0), reverse=True)
    return pool

def conservative_alt(a):
    g = AA_GROUP.get(a, AA)
    opts = [x for x in g if x != a]
    return random.choice(opts or [x for x in AA if x != a])

def mutate(seq, k=1, conservative=True):
    arr = list(seq)
    positions = [i for i in range(len(arr)) if (i+1) not in FIXED]
    for pos in random.sample(positions, k=min(k, len(positions))):
        arr[pos] = conservative_alt(arr[pos]) if conservative else random.choice([x for x in AA if x != arr[pos]])
    return "".join(arr)

def crossover(a, b, cuts):
    out = []
    last = 0
    use_a = True
    for c in cuts + [min(len(a), len(b))]:
        out.append((a if use_a else b)[last:c])
        use_a = not use_a
        last = c
    s = "".join(out)
    if len(s) < len(a): s += a[len(s):]
    return s

def build_candidates(pool):
    if not pool:
        raise RuntimeError("No parent pool found")
    parents = [x["seq"] for x in pool[:8]]
    top = parents[0]
    candidates = []
    def add(seq, direction, parent_label):
        if not seq.startswith("M") or len(seq) < 220 or len(seq) > 250: return
        candidates.append({"seq": seq, "direction": direction, "parent": parent_label})

    # A: single conservative mutations around top 2
    for pi, p in enumerate(parents[:2]):
        for _ in range(35):
            add(mutate(p, 1, True), "A_single_conservative", f"pool_p{pi+1}")
    # B: double conservative
    for pi, p in enumerate(parents[:2]):
        for _ in range(25):
            add(mutate(p, 2, True), "B_double_conservative", f"pool_p{pi+1}")
    # C: crossovers among top candidates
    cuts_list = [[40], [80], [120], [160], [60,140], [90,180], [45,100,170]]
    for i in range(min(4, len(parents))):
        for j in range(i+1, min(6, len(parents))):
            for cuts in cuts_list[:4]:
                add(crossover(parents[i], parents[j], cuts), "C_segment_crossover", f"p{i+1}xp{j+1}")
                add(crossover(parents[j], parents[i], cuts), "C_segment_crossover", f"p{j+1}xp{i+1}")
    # D: consensus-frequency sampling
    L = min(len(p) for p in parents[:8])
    freqs = []
    for pos in range(L):
        counts = {}
        for p in parents[:8]:
            counts[p[pos]] = counts.get(p[pos], 0) + 1
        freqs.append(sorted(counts.items(), key=lambda x: x[1], reverse=True))
    for _ in range(35):
        arr = list(top[:L])
        for pos in range(L):
            if (pos+1) in FIXED: continue
            if random.random() < 0.035 and len(freqs[pos]) > 1:
                arr[pos] = random.choice([a for a,c in freqs[pos][1:]])
        add("".join(arr), "D_consensus_sampling", "top_pool")
    # E: mixed local random non-conservative small jumps
    for _ in range(25):
        add(mutate(random.choice(parents[:4]), random.choice([1,2,3]), False), "E_small_random_jump", "top4")

    # Dedup and trim
    seen = set(); out = []
    for c in candidates:
        if c["seq"] in seen: continue
        seen.add(c["seq"]); out.append(c)
    random.shuffle(out)
    return out[:TARGET_N]

print("Loading parent pool...", flush=True)
pool = best_pool()
print(f"Parent pool size={len(pool)}; top={pool[0].get('score')} parent={pool[0].get('parent')}", flush=True)
candidates = build_candidates(pool)
print(f"Generated {len(candidates)} unique direct variants", flush=True)
json.dump(candidates, open(WORK/"generated_candidates.json", "w"), indent=2)

print("Loading ESMFold chunk=64...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(64); model.eval()
print("Loaded.", flush=True)

def predict(seq):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k:v.cuda() for k,v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES)
    plddt = out.plddt[0,:,1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {"ptm":round(ptm,4),"plddt":round(gp,4),"chromo":round(cp,4),"score":round(score,4),"passes":ptm>0.6 and gp>0.6 and cp>0.55}

t0 = time.time()
passed = []
for i,c in enumerate(candidates):
    try:
        m = predict(c["seq"])
    except Exception:
        torch.cuda.empty_cache(); continue
    if m["passes"]:
        m.update(c); m["length"] = len(c["seq"]); m["recycles"] = RECYCLES
        passed.append(m)
    torch.cuda.empty_cache()
    if (i+1) % 20 == 0:
        top = max([x["score"] for x in passed], default=0)
        print(f"[{time.strftime('%H:%M:%S')}] {i+1}/{len(candidates)}, passed={len(passed)}, top={top:.4f}", flush=True)

passed.sort(key=lambda x:x["score"], reverse=True)
json.dump(passed, open(WORK/"all_passed.json", "w"), indent=2)
final6 = passed[:6]
json.dump(final6, open(WORK/"final_6_r28.json", "w"), indent=2)
with open(WORK/"submission_r28.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Team_Name", "Seq_ID", "Sequence"])
    for i,c in enumerate(final6): w.writerow(["SnowFold", i+1, c["seq"]])

print("\nTop 12:", flush=True)
for i,c in enumerate(passed[:12]):
    print(f"{i+1:2d}. score={c['score']:.4f} pTM={c['ptm']:.4f} chromo={c['chromo']:.3f} dir={c['direction']} parent={c['parent']}", flush=True)
print(f"Done: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
if passed:
    print(f"R28 Top1={passed[0]['score']:.4f} vs R25=0.9477 delta={passed[0]['score']-0.9477:+.4f}", flush=True)
