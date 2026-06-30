#!/usr/bin/env python3
"""R22 Long-Run Pipeline (12-14 hours total)

Phase 1: R20 finalize - 评估已生成的 2000 候选 (2 父代 × 1000)
Phase 2: R21 large-scale - 用 R20 Top 6 父代, 4 温度 × 150 候选 × 6 = 3600 候选
Phase 3: Top 50 排序 + 多 recycles 共识 + r=20 重算 + 综合排名
"""
import os, sys, json, time, glob, subprocess, copy, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding, AutoModel

WORK = "/root/autodl-tmp/r22"
R20_WORK = "/root/autodl-tmp/r20"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "pdbs_precise", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

# === Config ===
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
TOP_K_PRECISE = 20
FIXED = [1, 65, 66, 67, 96, 222]
TEMPS_R21 = [0.1, 0.2, 0.5, 1.0]
NUM_SEQ_PER_TEMP_R21 = 150  # 4*150*6 = 3600
BATCH = 25

print("Loading ESMFold...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

def predict(seq, recycles=RECYCLES_SCREEN):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {"ptm": round(ptm,4), "plddt": round(gp,4), "chromo": round(cp,4),
            "score": round(score,4), "passes": ptm>0.60 and gp>0.60 and cp>0.55}

def predict_and_save_pdb(seq, pdb_path, recycles=RECYCLES_SCREEN):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    metrics = {"ptm":round(ptm,4), "plddt":round(gp,4), "chromo":round(cp,4),
               "score":round(score,4), "passes":ptm>0.60 and gp>0.60 and cp>0.55}
    positions = out.positions[-1][0].cpu().numpy()
    with open(pdb_path, "w") as f:
        f.write("REMARK  R22\n"); aidx=1
        for i in range(len(seq)):
            rn = aa3.get(seq[i], 'ALA')
            for j, an in enumerate(["N","CA","C","O"]):
                x,y,z = positions[i,j]
                f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt[i]*100:6.2f}\n")
                aidx+=1
        f.write("END\n")
    return metrics

def parse_fa(paths):
    seqs = []
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as f:
            n, s = "", ""
            for l in f:
                l = l.strip()
                if l.startswith(">"):
                    if s: seqs.append({"name":n,"seq":s})
                    n, s = l[1:], ""
                else:
                    if l: s += l
            if s: seqs.append({"name":n,"seq":s})
    return seqs

def list_fa_files(outdir):
    fa_files = []
    seqs_dir = os.path.join(outdir, "seqs")
    if not os.path.isdir(seqs_dir): return []
    for root, dirs, files in os.walk(seqs_dir):
        for f in files:
            if f.endswith(".fa"):
                fa_files.append(os.path.join(root, f))
    return sorted(fa_files)

def run_mpnn(pdb_path, name, num_seq, temps):
    outdir = os.path.join(WORK, "mpnn_out", name)
    os.makedirs(outdir, exist_ok=True)
    flag = os.path.join(outdir, "done.flag")
    if os.path.isfile(flag):
        files = list_fa_files(outdir)
        print(f"  [MPNN cached] {len(files)} files", flush=True)
        return files
    fixed_json = os.path.join(outdir, "fixed.jsonl")
    pdb_key = os.path.basename(pdb_path).replace(".pdb", "")
    with open(fixed_json, "w") as f:
        f.write(json.dumps({pdb_key: {"A": FIXED}}) + "\n")
    cmd = [sys.executable, os.path.join(MPNN, "protein_mpnn_run.py"),
           "--pdb_path", pdb_path, "--pdb_path_chains", "A",
           "--path_to_model_weights", os.path.join(MPNN, "vanilla_model_weights"),
           "--fixed_positions_jsonl", fixed_json, "--out_folder", outdir,
           "--num_seq_per_target", str(num_seq), "--batch_size", str(BATCH),
           "--sampling_temp", " ".join(str(t) for t in temps),
           "--seed", "42", "--suppress_print", "1"]
    print(f"  [MPNN] {name}...", end="", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    files = list_fa_files(outdir)
    if not files:
        print(f" retrying...", end="", flush=True)
        subprocess.run(cmd, timeout=1800)
        files = list_fa_files(outdir)
    if files:
        open(flag, "w").close()
        print(f" {len(files)} files", flush=True)
    else:
        print(f" FAILED (rc={r.returncode})", flush=True)
    return files

def screen_candidates(seqs, parent_label, recycles=RECYCLES_SCREEN):
    results = []
    n_total = len(seqs)
    print(f"  [screen] {n_total} candidates @ r={recycles} (parent: {parent_label})", flush=True)
    for i, e in enumerate(seqs):
        s = e["seq"]
        if not s.startswith("M") or len(s) < 220 or len(s) > 250:
            continue
        try:
            m = predict(s, recycles)
        except Exception as ex:
            torch.cuda.empty_cache()
            continue
        if not m["passes"]:
            continue
        m["name"] = e["name"][:60]
        m["seq"] = s
        m["length"] = len(s)
        m["parent"] = parent_label
        m["recycles"] = recycles
        results.append(m)
        torch.cuda.empty_cache()
        if (i+1) % 100 == 0:
            print(f"    [{time.strftime('%H:%M:%S')}] {i+1}/{n_total}, {len(results)} passed, VRAM={torch.cuda.memory_allocated()/1024**3:.1f}GB", flush=True)
    results.sort(key=lambda x: x["score"], reverse=True)
    print(f"  Result: {len(results)}/{n_total} passed ({100*len(results)/max(1,n_total):.1f}%)", flush=True)
    return results

if __name__ == "__main__":
    t0 = time.time()
    print(f"="*70, flush=True)
    print(f"R22 LONG Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"Phase 1: R20 finalize (2000 候选)")
    print(f"Phase 2: R21 large-scale (3600 候选)")
    print(f"Phase 3: Top 50 + r=20 重算")
    print(f"="*70, flush=True)
    
    all_results = []
    
    # ==================== Phase 1: R20 finalize ====================
    print(f"\n[Phase 1] R20 finalize - 评估 R20 已生成的 2000 候选", flush=True)
    R19_TOP6 = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))
    R19_BY_NAME = {p["name"]: p for p in R19_TOP6}
    
    r20_mpnn_dirs = sorted(glob.glob(os.path.join(R20_WORK, "mpnn_out", "*")))
    print(f"Found {len(r20_mpnn_dirs)} R20 MPNN outputs", flush=True)
    
    for mpnn_dir in r20_mpnn_dirs:
        parent_name = os.path.basename(mpnn_dir)
        parent = R19_BY_NAME.get(parent_name)
        if parent is None:
            print(f"  WARN: skip {parent_name}", flush=True)
            continue
        ps = parent["seq"]
        fa_files = list_fa_files(mpnn_dir)
        all_seqs = parse_fa(fa_files)
        filt = [x for x in all_seqs if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"\n[R20 finalize] {parent_name[:40]}: {len(filt)} candidates", flush=True)
        r = screen_candidates(filt, parent_name, RECYCLES_SCREEN)
        all_results.extend(r)
        # Save progress
        json.dump(all_results, open(os.path.join(WORK, "results", "phase1_progress.json"), "w"), indent=2)
    
    print(f"\n[Phase 1 DONE] Total passed: {len(all_results)}", flush=True)
    all_results.sort(key=lambda x: x["score"], reverse=True)
    
    # Take Top 6 from R20 as R21 parents
    R20_TOP6 = all_results[:6]
    json.dump(R20_TOP6, open(os.path.join(WORK, "results", "r20_top6.json"), "w"), indent=2)
    print(f"[Phase 1] R20 Top 6 saved. Top1={R20_TOP6[0]['score']:.4f}", flush=True)
    
    # ==================== Phase 2: R21 large-scale ====================
    print(f"\n[Phase 2] R21 large-scale - 用 R20 Top 6 父代, 4 温度 × 150 = 600/父代", flush=True)
    
    for i, parent in enumerate(R20_TOP6):
        pn = parent.get("name", f"p{i}")[:30]
        ps = parent["seq"]
        print(f"\n[R21] Parent {i+1}/6: {pn} ({len(ps)}aa)", flush=True)
        
        # ESMFold PDB
        pdb = os.path.join(WORK, "pdbs", f"r21_p{i+1}.pdb")
        if not os.path.isfile(pdb):
            print(f"  ESMFold r=8 → PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb, RECYCLES_SCREEN)
            print("done", flush=True)
        
        # MPNN
        ffs = run_mpnn(pdb, f"r21_p{i+1}", NUM_SEQ_PER_TEMP_R21, TEMPS_R21)
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"  filtered: {len(filt)}", flush=True)
        if not filt:
            continue
        
        # ESMFold evaluate
        r = screen_candidates(filt, f"r21_p{i+1}", RECYCLES_SCREEN)
        all_results.extend(r)
        # Save progress
        json.dump(all_results, open(os.path.join(WORK, "results", "phase2_progress.json"), "w"), indent=2)
    
    print(f"\n[Phase 2 DONE] Total passed (R20+R21): {len(all_results)}", flush=True)
    all_results.sort(key=lambda x: x["score"], reverse=True)
    json.dump(all_results, open(os.path.join(WORK, "results", "all_passed.json"), "w"), indent=2)
    
    # ==================== Phase 3: Top 50 + r=20 重算 ====================
    print(f"\n[Phase 3] Top 50 + r=20 high-precision recount", flush=True)
    top50 = all_results[:50]
    print(f"\nTop 10 after screening:", flush=True)
    for i, c in enumerate(top50[:10]):
        print(f"  {i+1:2d}. parent={c['parent'][:25]:25s} score={c['score']:.4f} pTM={c['ptm']:.4f}", flush=True)
    
    top20 = top50[:TOP_K_PRECISE]
    precise = []
    for i, c in enumerate(top20):
        nm = f"top{i+1:02d}_p{c['parent'][:10]}_s{c['score']:.4f}"
        pdb_precise = os.path.join(WORK, "pdbs_precise", f"{nm}.pdb")
        m_precise = predict_and_save_pdb(c["seq"], pdb_precise, RECYCLES_PRECISE)
        m_precise["name"] = c["name"]
        m_precise["seq"] = c["seq"]
        m_precise["parent"] = c["parent"]
        ptm_diff = abs(m_precise["ptm"] - c["ptm"])
        if ptm_diff > 0.05:
            if m_precise["score"] < c["score"]:
                chosen = m_precise
                chosen["note"] = f"r=20 chosen (diff={ptm_diff:.3f})"
            else:
                chosen = c
                chosen["note"] = f"r=8 kept (diff={ptm_diff:.3f})"
        else:
            chosen = m_precise
            chosen["note"] = f"r=20 OK (diff={ptm_diff:.3f})"
        chosen["recycles"] = RECYCLES_PRECISE
        chosen["ptm_diff"] = round(ptm_diff, 4)
        precise.append(chosen)
        print(f"  {i+1:2d}. r=20: score={chosen['score']:.4f} pTM={chosen['ptm']:.4f} (diff={ptm_diff:.4f}, {chosen['note']})", flush=True)
    
    pool = precise + top50[TOP_K_PRECISE:]
    pool.sort(key=lambda x: x["score"], reverse=True)
    final6 = pool[:6]
    
    import csv
    csv_path = os.path.join(WORK, "submission_r22.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for i, c in enumerate(final6):
            w.writerow(["SnowFold", i+1, c["seq"]])
            print(f"  FINAL {i+1}: score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f} r={c.get('recycles',8)}", flush=True)
    
    json.dump(final6, open(os.path.join(WORK, "final_6_r22.json"), "w"), indent=2)
    
    # Compare with R19
    r19_top = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))[0]["score"]
    r22_top = final6[0]["score"]
    print(f"\n=== R19 vs R22 ===", flush=True)
    print(f"  R19 Top 1: {r19_top:.4f}", flush=True)
    print(f"  R22 Top 1: {r22_top:.4f}", flush=True)
    delta = r22_top - r19_top
    print(f"  Δ: {delta:+.4f} ({(delta/r19_top)*100:+.2f}%)", flush=True)
    print(f"  Status: {'✅ BREAKTHROUGH' if delta > 0 else '⚠️ NO IMPROVEMENT'}", flush=True)
    
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"Submit: {csv_path}", flush=True)
