#!/usr/bin/env python3
"""R20 Local Consensus Scoring - Multi-recycles voting

对 R19 + R20 Top 6 (共 12 条) 用 r=4/8/12/16/20 多档评估
输出共识分数（加权平均）+ 单 recycles 详细
"""
import os, sys, json, time, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "D:/workspace/round20"
os.makedirs(WORK, exist_ok=True)

RECYCLES_LIST = [4, 8, 12, 16, 20]  # 5 档共识

print("Loading ESMFold (local RTX 5080)...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

def predict(seq, recycles):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {
        "ptm": round(ptm, 4), "plddt": round(gp, 4), "chromo": round(cp, 4),
        "score": round(score, 4),
        "passes": ptm > 0.60 and gp > 0.60 and cp > 0.55,
    }

if __name__ == "__main__":
    t0 = time.time()
    print(f"R20 Consensus Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Load R19 + R20 Top 6
    r19_top6 = json.load(open(r'D:\workspace\round19\final_6_r19.json'))
    r20_top6 = json.load(open(r'D:\workspace\round20\r20_top6.json'))

    # Mark origins
    for c in r19_top6: c['origin'] = 'R19'
    for c in r20_top6: c['origin'] = 'R20'

    # Combine unique sequences
    all_cands = {}
    for c in r19_top6 + r20_top6:
        if c['seq'] not in all_cands:
            all_cands[c['seq']] = c
    candidates = list(all_cands.values())
    print(f"Total unique candidates: {len(candidates)}", flush=True)

    # Multi-recycles scoring
    results = []
    for i, c in enumerate(candidates):
        seq = c['seq']
        print(f"\n[{i+1}/{len(candidates)}] origin={c['origin']} score={c['score']:.4f}", flush=True)
        per_r = {}
        for r in RECYCLES_LIST:
            t1 = time.time()
            m = predict(seq, r)
            t2 = time.time()
            per_r[r] = m
            print(f"  r={r}: pTM={m['ptm']:.4f} pLDDT={m['plddt']:.3f} chromo={m['chromo']:.3f} score={m['score']:.4f} ({t2-t1:.1f}s)", flush=True)
            torch.cuda.empty_cache()

        # Compute consensus score (weighted average favoring higher recycles)
        scores = np.array([per_r[r]['score'] for r in RECYCLES_LIST])
        weights = np.array([0.1, 0.15, 0.2, 0.25, 0.30])  # higher recycles = higher weight
        weights = weights / weights.sum()
        consensus_score = float(np.dot(scores, weights))

        # Best single r
        best_r = RECYCLES_LIST[np.argmax(scores)]
        best_score = float(scores.max())

        # Std (consistency indicator)
        score_std = float(scores.std())

        print(f"  Consensus (weighted): {consensus_score:.4f}", flush=True)
        print(f"  Best single: r={best_r} score={best_score:.4f}", flush=True)
        print(f"  Std: {score_std:.4f}", flush=True)

        result = dict(c)
        result['per_recycles'] = per_r
        result['consensus_score'] = round(consensus_score, 4)
        result['best_r'] = best_r
        result['best_score'] = round(best_score, 4)
        result['score_std'] = round(score_std, 4)
        results.append(result)

    # Sort by consensus
    results.sort(key=lambda x: x['consensus_score'], reverse=True)

    print(f"\n{'='*80}")
    print(f"{'Multi-Recycles Consensus Top 6':^80s}")
    print(f"{'='*80}")
    print(f"{'Origin':>7s} {'Orig':>7s} {'Consen':>7s} {'Best':>7s} {'Std':>7s} {'Best_r':>7s}")
    for i, r in enumerate(results):
        print(f"{r['origin']:>7s} {r['score']:>7.4f} {r['consensus_score']:>7.4f} {r['best_score']:>7.4f} {r['score_std']:>7.4f} {r['best_r']:>7d}")

    # Save
    json.dump(results, open(os.path.join(WORK, 'consensus_results.json'), 'w'), indent=2)

    # Generate submission CSV
    import csv
    csv_path = os.path.join(WORK, 'submission_r20_consensus.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Team_Name', 'Seq_ID', 'Sequence'])
        for i, r in enumerate(results[:6]):
            w.writerow(['SnowFold', i+1, r['seq']])
    print(f"\nSubmission: {csv_path}", flush=True)

    # Compare R20 single r=8 vs consensus
    r20_only = [r for r in results if r['origin'] == 'R20']
    r20_only.sort(key=lambda x: x['consensus_score'], reverse=True)
    print(f"\nR20 only consensus (by consensus_score):")
    for i, r in enumerate(r20_only):
        print(f"  {i+1}. {r['parent'][:25]:25s} R20_score={r['score']:.4f} Consensus={r['consensus_score']:.4f}", flush=True)

    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)