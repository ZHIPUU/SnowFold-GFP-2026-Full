#!/usr/bin/env python3
"""R20 Pure CPU Sequence Analysis - No external models

Computes comprehensive sequence features for R20 Top 6 + R22 candidates.
Uses only Python stdlib + numpy.
"""
import os, sys, json, time, csv
from collections import Counter
import numpy as np

WORK_REMOTE = "/root/autodl-tmp/r20_analysis"
os.makedirs(WORK_REMOTE, exist_ok=True)

# Reference sequences
SFGFP_WT = "MASKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
avGFP_WT = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

# Amino acid properties
AA_HYDROPHOBIC = set("AILMFWVP")
AA_AROMATIC = set("FWY")
AA_POLAR = set("STNQ")
AA_POS = set("RKH")
AA_NEG = set("DE")
AA_TINY = set("AGC")
AA_ALL = "ACDEFGHIKLMNPQRSTVWY"
VALID_AA = set(AA_ALL)

def composition(seq):
    c = Counter(seq)
    L = len(seq)
    return {a: c.get(a, 0) / L for a in AA_ALL}

def features(seq):
    L = len(seq)
    c = Counter(seq)
    return {
        "length": L,
        "M_start": seq[0] == "M",
        "hydrophobic_pct": 100 * sum(c[a] for a in AA_HYDROPHOBIC) / L,
        "aromatic_pct": 100 * sum(c[a] for a in AA_AROMATIC) / L,
        "polar_pct": 100 * sum(c[a] for a in AA_POLAR) / L,
        "positive_charge": sum(c[a] for a in AA_POS),
        "negative_charge": sum(c[a] for a in AA_NEG),
        "net_charge": sum(c[a] for a in AA_POS) - sum(c[a] for a in AA_NEG),
        "tiny_pct": 100 * sum(c[a] for a in AA_TINY) / L,
        "cys_count": c.get("C", 0),  # disulfide potential
        "pro_count": c.get("P", 0),  # structural constraints
        "gly_count": c.get("G", 0),  # flexibility
        "invalid_count": sum(1 for a in seq if a not in VALID_AA),
    }

def hamming(s1, s2):
    """Hamming distance (assumes aligned to min length)"""
    L = min(len(s1), len(s2))
    return sum(1 for a, b in zip(s1[:L], s2[:L]) if a != b)

def identity_aligned(s1, s2):
    L = min(len(s1), len(s2))
    matches = sum(1 for a, b in zip(s1[:L], s2[:L]) if a == b)
    return 100 * matches / L

def mutation_summary(s, ref):
    """List mutations by position"""
    L = min(len(s), len(ref))
    muts = []
    for i in range(L):
        if s[i] != ref[i]:
            muts.append((i+1, ref[i], s[i]))
    return muts

def mutation_signature(s, ref):
    """Single-letter mutation signature like 'S30R,N105T'"""
    return ",".join([f"{ref[i]}{i+1}{s[i]}" for i in range(min(len(s), len(ref))) if s[i] != ref[i]])

def chromophore_check(s):
    """Check 5 core residues at positions 65, 66, 67, 96, 222 (0-indexed: 64, 65, 66, 95, 221)"""
    if len(s) < 223:
        return None
    return {
        "pos65": s[64],  # T (should be T or S)
        "pos66": s[65],  # Y (must be Y)
        "pos67": s[66],  # G (must be G)
        "pos96": s[95],  # R (should be R or K)
        "pos222": s[221], # E (should be E or D)
    }

if __name__ == "__main__":
    t0 = time.time()
    print(f"Pure CPU analysis start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    
    # Load R20 Top 6
    with open("/root/autodl-tmp/r22/results/r20_top6.json") as f:
        r20_top6 = json.load(f)
    
    # Load R22 passed candidates
    passed = []
    for path in ["/root/autodl-tmp/r22/results/phase1_progress.json",
                 "/root/autodl-tmp/r22/results/phase2_progress.json"]:
        if os.path.exists(path):
            with open(path) as f:
                passed.extend(json.load(f))
    seen = set()
    unique = []
    for c in passed:
        if c["seq"] not in seen:
            seen.add(c["seq"])
            unique.append(c)
    passed = unique
    
    print(f"R20 Top 6: {len(r20_top6)}", flush=True)
    print(f"R20+R22 passed unique: {len(passed)}", flush=True)
    
    # Analyze R20 Top 6
    results = []
    for i, c in enumerate(r20_top6):
        seq = c["seq"]
        f = features(seq)
        h_sf = hamming(seq, SFGFP_WT)
        h_av = hamming(seq, avGFP_WT)
        id_sf = identity_aligned(seq, SFGFP_WT)
        id_av = identity_aligned(seq, avGFP_WT)
        muts_sf = mutation_summary(seq, SFGFP_WT)
        mut_sig = mutation_signature(seq, SFGFP_WT)
        chromo = chromophore_check(seq)
        results.append({
            "rank": i+1,
            "parent": c["parent"][:35],
            "esmfold_score": c["score"],
            "ptm": c["ptm"],
            "plddt": c["plddt"],
            "chromo": c["chromo"],
            **f,
            "vs_sfGFP_hamming": h_sf,
            "vs_sfGFP_identity_pct": round(id_sf, 2),
            "vs_avGFP_hamming": h_av,
            "vs_avGFP_identity_pct": round(id_av, 2),
            "mutations_vs_sfGFP": mut_sig,
            "mut_count_vs_sfGFP": len(muts_sf),
            "chromophore_5core": chromo,
        })
    
    # Print
    print("\n" + "="*100)
    print("R20 Top 6 Pure CPU Sequence Analysis")
    print("="*100)
    for r in results:
        print(f"\nSeq {r['rank']}: score={r['esmfold_score']:.4f} parent={r['parent']}", flush=True)
        print(f"  Length: {r['length']}, Hydrophobic: {r['hydrophobic_pct']:.1f}%, Net charge: {r['net_charge']:+d}", flush=True)
        print(f"  vs sfGFP: hamming={r['vs_sfGFP_hamming']} ({100*r['vs_sfGFP_hamming']/r['length']:.1f}%), identity={r['vs_sfGFP_identity_pct']:.1f}%", flush=True)
        print(f"  vs avGFP: hamming={r['vs_avGFP_hamming']}, identity={r['vs_avGFP_identity_pct']:.1f}%", flush=True)
        print(f"  Chromophore 5-core: {r['chromophore_5core']}", flush=True)
        print(f"  Mutations vs sfGFP: {r['mutations_vs_sfGFP'][:150]}{'...' if len(r['mutations_vs_sfGFP'])>150 else ''}", flush=True)
    
    # Summary
    print("\n" + "="*100)
    print("Summary")
    print("="*100)
    h_sf = [r["vs_sfGFP_hamming"] for r in results]
    h_av = [r["vs_avGFP_hamming"] for r in results]
    muts = [r["mut_count_vs_sfGFP"] for r in results]
    print(f"vs sfGFP:  hamming min={min(h_sf)}, max={max(h_sf)}, avg={sum(h_sf)/6:.1f}", flush=True)
    print(f"vs avGFP:  hamming min={min(h_av)}, max={max(h_av)}, avg={sum(h_av)/6:.1f}", flush=True)
    print(f"Mutations vs sfGFP: min={min(muts)}, max={max(muts)}, avg={sum(muts)/6:.1f}", flush=True)
    
    # Common mutations
    all_muts = []
    for r in results:
        all_muts.extend(r["mutations_vs_sfGFP"].split(",") if r["mutations_vs_sfGFP"] else [])
    common = Counter(all_muts).most_common(10)
    print(f"\nTop 10 most common mutations vs sfGFP:", flush=True)
    for mut, count in common:
        print(f"  {mut}: {count}/6", flush=True)
    
    # Save CSV
    csv_path = f"{WORK_REMOTE}/r20_analysis.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print(f"\nCSV: {csv_path}", flush=True)
    
    json.dump(results, open(f"{WORK_REMOTE}/r20_analysis.json", "w"), indent=2)
    
    # Build frequency matrix across all passed candidates
    print("\n" + "="*100)
    print("All R20+R22 unique passed candidates - aggregate analysis")
    print("="*100)
    all_seqs = [c["seq"] for c in passed]
    all_freqs = np.zeros((len(AA_ALL), len(all_seqs)))
    for j, seq in enumerate(all_seqs):
        c = Counter(seq)
        for i, a in enumerate(AA_ALL):
            all_freqs[i, j] = c.get(a, 0) / len(seq)
    mean_freqs = all_freqs.mean(axis=1)
    std_freqs = all_freqs.std(axis=1)
    
    print(f"\n{'AA':>3s} {'Mean%':>6s} {'Std%':>6s} {'sfGFP%':>7s}", flush=True)
    for i, a in enumerate(AA_ALL):
        sf_pct = SFGFP_WT.count(a) / len(SFGFP_WT) * 100
        print(f"{a:>3s} {mean_freqs[i]*100:>5.2f}% {std_freqs[i]*100:>5.2f}% {sf_pct:>6.2f}%", flush=True)
    
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)