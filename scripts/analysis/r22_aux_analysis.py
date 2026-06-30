#!/usr/bin/env python3
"""R22 Auxiliary Analysis: ESM2-650M likelihood + sequence analysis for R20 Top 6

Runs on CPU (not GPU) - parallel with R22 long-run on A800.
"""
import os, sys, json, time, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, AutoModelForMaskedLM

WORK = "/root/autodl-tmp/r22_aux"
os.makedirs(WORK, exist_ok=True)

print("Loading ESM2-650M (CPU)...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True)
model = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True).eval()
print("Loaded.", flush=True)

def compute_log_likelihood(seq):
    """Compute mean log-likelihood per residue (lower = less natural)"""
    tokens = tok(seq, return_tensors="pt")
    input_ids = tokens["input_ids"]
    with torch.no_grad():
        outputs = model(input_ids)
        logits = outputs.logits  # (1, L, vocab_size)
    # Shift: predict token i+1 from tokens 0..i
    log_probs = F.log_softmax(logits, dim=-1)
    target = input_ids[:, 1:]   # (1, L-1)
    pred_logp = log_probs[:, :-1, :].gather(2, target.unsqueeze(-1)).squeeze(-1)
    return float(pred_logp.mean().item())

def sequence_features(seq):
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return {
        "len": len(seq),
        "M_start": seq[0] == "M",
        "charge": seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E"),
        "hydrophobic": sum(seq.count(a) for a in "AILMFVW") / len(seq),
        "polar": sum(seq.count(a) for a in "NQSTY") / len(seq),
        "Gly_count": seq.count("G"),
        "Pro_count": seq.count("P"),
        "Cys_count": seq.count("C"),
        "aa_violation": sum(1 for a in seq if a not in aa),
    }

if __name__ == "__main__":
    t0 = time.time()
    print(f"R22 aux start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    
    # Load R20 Top 6
    with open("/root/autodl-tmp/r22/results/r20_top6.json") as f:
        top6 = json.load(f)
    
    # Load sfGFP WT for comparison
    SFGFP_WT = "MASKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    
    # Compute ESM2 log-likelihood for sfGFP WT
    print("\nComputing sfGFP WT log-likelihood...", flush=True)
    wt_logp = compute_log_likelihood(SFGFP_WT)
    print(f"  sfGFP WT logP: {wt_logp:.4f}", flush=True)
    
    # Compute for each R20 Top 6
    results = []
    for i, c in enumerate(top6):
        seq = c["seq"]
        print(f"\n[{i+1}/6] {c['parent'][:30]} score={c['score']:.4f}", flush=True)
        logp = compute_log_likelihood(seq)
        delta = logp - wt_logp
        print(f"  ESM2-650M logP: {logp:.4f} (delta vs WT: {delta:+.4f})", flush=True)
        feats = sequence_features(seq)
        print(f"  Features: {feats}", flush=True)
        result = dict(c)
        result["esm2_logp"] = round(logp, 4)
        result["esm2_delta_vs_wt"] = round(delta, 4)
        result["features"] = feats
        results.append(result)
    
    # Combined score: 50% ESMFold + 50% ESM2 normalized
    esmf_scores = [c["score"] for c in results]
    esm2_scores = [c["esm2_logp"] for c in results]
    
    # Min-max normalize
    def minmax(xs):
        mn, mx = min(xs), max(xs)
        if mx - mn < 1e-9: return [0.5]*len(xs)
        return [(x - mn)/(mx - mn) for x in xs]
    
    n_esmf = minmax(esmf_scores)
    n_esm2 = minmax(esm2_scores)
    for i, c in enumerate(results):
        c["combined_score"] = round(0.5*n_esmf[i] + 0.5*n_esm2[i], 4)
    
    # Sort by combined
    results.sort(key=lambda x: x["combined_score"], reverse=True)
    
    # Save
    json.dump(results, open(os.path.join(WORK, "aux_results.json"), "w"), indent=2)
    json.dump({"wt_logp": wt_logp}, open(os.path.join(WORK, "wt_logp.json"), "w"))
    
    print(f"\n=== R20 Top 6 + ESM2-650M likelihood ===", flush=True)
    print(f"{'Seq':>4s} {'ESMFold':>7s} {'ESM2_logP':>10s} {'Δ_vs_WT':>8s} {'Combined':>9s}", flush=True)
    for i, c in enumerate(results):
        print(f"{i+1:>4d} {c['score']:>7.4f} {c['esm2_logp']:>10.4f} {c['esm2_delta_vs_wt']:>+8.4f} {c['combined_score']:>9.4f}", flush=True)
    
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)