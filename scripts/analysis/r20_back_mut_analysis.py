#!/usr/bin/env python3
"""R20 Local CPU Tasks (parallel with GPU consensus)

Task A: 突变关键性 back-mutation 分析
Task B: 候选相似性网络 + 树状图

These run on local CPU - no GPU needed.
"""
import os, sys, json, time, csv
from collections import Counter
import numpy as np

WORK = "D:/workspace/round20"
SFGFP = "MASKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

def hamming(s1, s2):
    L = min(len(s1), len(s2))
    return sum(1 for a, b in zip(s1[:L], s2[:L]) if a != b)

def get_mutations(s, ref):
    L = min(len(s), len(ref))
    muts = []
    for i in range(L):
        if s[i] != ref[i]:
            muts.append((i+1, ref[i], s[i]))
    return muts

def back_mutate(seq, ref, positions):
    """回退到 ref 在 positions 列表的残基"""
    seq = list(seq)
    for pos in positions:
        if pos - 1 < len(seq) and pos - 1 < len(ref):
            seq[pos-1] = ref[pos-1]
    return ''.join(seq)

def similarity_matrix(seqs):
    n = len(seqs)
    m = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            d = hamming(seqs[i], seqs[j])
            m[i, j] = d
            m[j, i] = d
    return m

if __name__ == "__main__":
    t0 = time.time()
    print(f"R20 CPU parallel tasks Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Load candidates
    r19_top6 = json.load(open(r'D:\workspace\round19\final_6_r19.json'))
    r20_top6 = json.load(open(r'D:\workspace\round20\r20_top6.json'))
    for c in r19_top6: c['origin'] = 'R19'
    for c in r20_top6: c['origin'] = 'R20'

    # Combine unique
    seen = {}
    for c in r19_top6 + r20_top6:
        if c['seq'] not in seen:
            seen[c['seq']] = c
    candidates = list(seen.values())
    seqs = [c['seq'] for c in candidates]
    print(f"Total unique: {len(candidates)}", flush=True)

    # ===== Task A: Back-mutation analysis (for R20 Top 1) =====
    print(f"\n{'='*80}")
    print("Task A: Back-mutation analysis (R20 Top 1)")
    print(f"{'='*80}")

    top1 = r20_top6[0]
    top1_seq = top1['seq']
    muts = get_mutations(top1_seq, SFGFP)
    print(f"R20 Top 1 has {len(muts)} mutations vs sfGFP", flush=True)

    # 模拟回退每一组突变并预测 score 变化（不需要 ESMFold，假设）
    # 实际上：back-mutation 会减少与 sfGFP 的距离，可能减分也可能加分
    # 没法预测 score，让我们做更精细的分析

    # 寻找"高密度突变区域"
    mut_positions = [p for p, r1, r2 in muts]
    windows = []
    for start in range(1, max(mut_positions)+1):
        cnt = sum(1 for p in mut_positions if abs(p - start) < 5)
        if cnt >= 8:
            windows.append((start, cnt))
    print(f"\nHigh-density mutation clusters (>=8 mutations within 5 residues):")
    for pos, cnt in sorted(set(windows))[:5]:
        nearby = [(p, r1, r2) for p, r1, r2 in muts if abs(p - pos) < 5]
        print(f"  Near pos {pos}: {cnt} mutations: {nearby[:5]}")

    # ===== Task B: Diversity analysis =====
    print(f"\n{'='*80}")
    print("Task B: Diversity analysis (12 candidates)")
    print(f"{'='*80}")

    sim_matrix = similarity_matrix(seqs)

    # All pairwise distances
    pairs = []
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            pairs.append((candidates[i]['origin'], candidates[j]['origin'],
                          i, j, int(sim_matrix[i, j])))
    pairs.sort(key=lambda x: x[4])

    print(f"\nMost similar pairs (lowest hamming):")
    for p1, p2, i, j, d in pairs[:10]:
        print(f"  {p1}[{i}] vs {p2}[{j}]: hamming={d} (identity={100-d/239:.1f}%)")

    # Diversity statistics
    within_r19 = [d for o1, o2, _, _, d in pairs if o1 == o2 == 'R19']
    within_r20 = [d for o1, o2, _, _, d in pairs if o1 == o2 == 'R20']
    cross = [d for o1, o2, _, _, d in pairs if o1 != o2]
    print(f"\nWithin R19: avg={np.mean(within_r19):.1f}, min={min(within_r19)}, max={max(within_r19)}")
    print(f"Within R20: avg={np.mean(within_r20):.1f}, min={min(within_r20)}, max={max(within_r20)}")
    print(f"Cross R19-R20: avg={np.mean(cross):.1f}, min={min(cross)}, max={max(cross)}")

    # Save results
    json.dump({
        'clusters': sorted(set(windows)),
        'diversity': {
            'within_r19': [int(x) for x in within_r19],
            'within_r20': [int(x) for x in within_r20],
            'cross': [int(x) for x in cross],
        },
        'sim_matrix_size': int(sim_matrix.shape[0]),
    }, open(os.path.join(WORK, 'cpu_analysis.json'), 'w'), indent=2)

    # ===== Task C: 共同突变详细统计 =====
    print(f"\n{'='*80}")
    print("Task C: Detailed mutation statistics (R20 Top 6)")
    print(f"{'='*80}")

    r20_seqs = [c['seq'] for c in r20_top6]
    all_muts = []  # List of (pos, ref_aa, mut_aa, count)

    for pos in range(1, 240):
        ref_aa = SFGFP[pos-1] if pos-1 < len(SFGFP) else 'X'
        # count how many R20 Top 6 have non-ref at this position
        non_ref_count = 0
        non_ref_seqs = []
        for i, seq in enumerate(r20_seqs):
            if pos-1 < len(seq) and seq[pos-1] != ref_aa:
                non_ref_count += 1
                non_ref_seqs.append((i+1, seq[pos-1]))
        if non_ref_count > 0:
            all_muts.append((pos, ref_aa, non_ref_count, non_ref_seqs))

    all_muts.sort(key=lambda x: x[2], reverse=True)

    print(f"\nTop 15 mutated positions (sorted by frequency):")
    print(f"{'Pos':>4s} {'RefAA':>6s} {'Count':>6s} Variations")
    for pos, ref_aa, cnt, vars in all_muts[:15]:
        vars_str = ', '.join([f'#{i+1}{a}' for i, a in vars])
        print(f"{pos:>4d} {ref_aa:>6s} {cnt}/6 {vars_str}")

    # Categorize by R20 consensus (all 6 same AA vs varied)
    conserved_muts = [(p, r, c) for p, r, c, v in all_muts if c == 6]
    varied_muts = [(p, r, c) for p, r, c, v in all_muts if c < 6]
    print(f"\nConserved (6/6 same position mutated): {len(conserved_muts)} positions")
    print(f"Varied (1-5/6 mutated): {len(varied_muts)} positions")
    print(f"Total mutations observed: {len(all_muts)} positions")

    # Save full data
    varied_list = []
    for m in all_muts:
        if m[2] < 6:
            pos, ref_aa, cnt, variants = m
            varied_list.append({
                'pos': pos,
                'ref': ref_aa,
                'count': cnt,
                'variants': [{'seq_id': i+1, 'aa': a} for i, a in variants]
            })
    json.dump({
        'r20_top1_mutations': [{'pos': p, 'ref': r, 'to': m} for p, r, m in muts],
        'r20_conserved_positions': [{'pos': p, 'ref': r, 'count': c} for p, r, c in conserved_muts],
        'r20_varied_positions': varied_list,
    }, open(os.path.join(WORK, 'mutation_analysis.json'), 'w'), indent=2)

    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)