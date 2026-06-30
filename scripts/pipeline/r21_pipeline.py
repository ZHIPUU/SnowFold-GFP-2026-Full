#!/usr/bin/env python3
"""R21 Pipeline - Use R20 Top 6 as parents, more aggressive MPNN

Strategy:
- Parents: R20 Top 6 (the best from R20's 1000 candidates)
- 4 temperatures, 250/temp = 1000/parent
- Top 50 → r=20 recount (project rule 5.3)
- Skip the WT failed path (R19 verified)
"""
import os, sys, json, time, glob, subprocess, copy, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r21"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "pdbs_precise", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

NUM_SEQ_PER_TEMP = 250
TEMPS = [0.1, 0.2, 0.5, 1.0]
FIXED = [1, 65, 66, 67, 96, 222]  # M + 5 chromophore
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
TOP_K_SCREEN = 50
TOP_K_PRECISE = 20
BATCH = 25

# Load R20 Top 6 (will be created after R20 completes)
R20_TOP6_PATH = "/root/autodl-tmp/r20/final_6_r20.json"
if os.path.exists(R20_TOP6_PATH):
    R20_TOP6 = json.load(open(R20_TOP6_PATH))
    print(f"Loaded R20 Top 6 ({len(R20_TOP6)} candidates)", flush=True)
else:
    # Fallback to R19 if R20 not done
    R20_TOP6 = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))
    print(f"WARN: R20 not done, using R19 Top 6", flush=True)

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
        f.write("REMARK  R21\n"); aidx=1
        for i in range(len(seq)):
            rn = aa3.get(seq[i], 'ALA')
            for j, an in enumerate(["N","CA","C","O"]):
                x,y,z = positions[i,j]
                f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt[i]*100:6.2f}\n")
                aidx+=1
        f.write("END\n")
    return metrics

def list_fa_files(outdir):
    fa_files = []
    seqs_dir = os.path.join(outdir, "seqs")
    if not os.path.isdir(seqs_dir):
        return []
    for root, dirs, files in os.walk(seqs_dir):
        for f in files:
            if f.endswith(".fa"):
                fa_files.append(os.path.join(root, f))
    return sorted(fa_files)

def run_mpnn(pdb_path, name):
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
           "--fixed_positions_jsonl", fixed_json,
           "--out_folder", outdir,
           "--num_seq_per_target", str(NUM_SEQ_PER_TEMP),
           "--batch_size", str(BATCH),
           "--sampling_temp", " ".join(str(t) for t in TEMPS),
           "--seed", "42", "--suppress_print", "1"]
    print(f"  [MPNN] {name}...", end="", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    files = list_fa_files(outdir)
    if not files:
        print(f" retrying...", end="", flush=True)
        subprocess.run(cmd, timeout=600)
        files = list_fa_files(outdir)
    if files:
        open(flag, "w").close()
        print(f" {len(files)} files", flush=True)
    else:
        print(f" FAILED (rc={r.returncode})", flush=True)
    return files

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

def screen_candidates(seqs, parent):
    results = []
    n_total = len(seqs)
    print(f"  [screen] {n_total} candidates @ r={RECYCLES_SCREEN}", flush=True)
    for i, e in enumerate(seqs):
        s = e["seq"]; n = e.get("name", f"{parent}_{i}")
        if not s.startswith("M") or len(s) < 220 or len(s) > 250: continue
        try:
            m = predict(s, RECYCLES_SCREEN)
        except Exception as ex:
            torch.cuda.empty_cache()
            continue
        if not m["passes"]: continue
        m["name"] = n; m["seq"] = s; m["length"] = len(s); m["parent"] = parent
        m["recycles"] = RECYCLES_SCREEN
        results.append(m)
        torch.cuda.empty_cache()
        if (i+1) % 50 == 0:
            print(f"    [{time.strftime('%H:%M:%S')}] screened {i+1}/{n_total}, {len(results)} passed, VRAM={torch.cuda.memory_allocated()/1024**3:.1f}GB", flush=True)
    results.sort(key=lambda x: x["score"], reverse=True)
    return results

def high_precision_recount(seq_data, name):
    seq = seq_data["seq"]
    pdb_precise = os.path.join(WORK, "pdbs_precise", f"{name}.pdb")
    m_precise = predict_and_save_pdb(seq, pdb_precise, RECYCLES_PRECISE)
    m_precise["name"] = name
    m_precise["seq"] = seq
    m_precise["length"] = len(seq)
    ptm_diff = abs(m_precise["ptm"] - seq_data["ptm"])
    if ptm_diff > 0.05:
        if m_precise["score"] < seq_data["score"]:
            chosen = m_precise
            chosen["note"] = f"r=20 chosen (diff={ptm_diff:.3f})"
        else:
            chosen = seq_data
            chosen["note"] = f"r=8 kept (diff={ptm_diff:.3f})"
    else:
        chosen = m_precise
        chosen["note"] = f"r=20 (diff={ptm_diff:.3f})"
    chosen["recycles"] = RECYCLES_PRECISE
    chosen["ptm_diff"] = round(ptm_diff, 4)
    return chosen

if __name__ == "__main__":
    t0 = time.time()
    print(f"="*70, flush=True)
    print(f"R21 Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"="*70, flush=True)
    
    all_results = []
    for parent in R20_TOP6:
        pn = parent.get("name", "unknown")[:30]; ps = parent["seq"]
        print(f"\n[{time.strftime('%H:%M:%S')}] Parent: {pn} ({len(ps)}aa)", flush=True)
        
        pdb = os.path.join(WORK, "pdbs", f"r21_p{i}.pdb")
        if not os.path.isfile(pdb):
            print(f"  ESMFold r=8 → PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb, RECYCLES_SCREEN)
            print("done", flush=True)
        
        ffs = run_mpnn(pdb, f"r21_p{i}")
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"  filtered: {len(filt)}", flush=True)
        if not filt: continue
        r = screen_candidates(filt, f"r21_p{i}")
        all_results.extend(r)
    
    all_results.sort(key=lambda x: x["score"], reverse=True)
    json.dump(all_results, open(os.path.join(WORK, "results", "all.json"), "w"), indent=2)
    
    print(f"\n=== Top 50 ===", flush=True)
    top50 = all_results[:TOP_K_SCREEN]
    for i, c in enumerate(top50[:10]):
        print(f"  {i+1}. score={c['score']:.4f} pTM={c['ptm']:.4f}", flush=True)
    
    print(f"\n=== Top {TOP_K_PRECISE} high-precision (r={RECYCLES_PRECISE}) ===", flush=True)
    top20 = top50[:TOP_K_PRECISE]
    precise = []
    for i, c in enumerate(top20):
        nm = f"top{i+1:02d}_score{c['score']:.4f}"
        m = high_precision_recount(c, nm)
        precise.append(m)
        print(f"  {i+1}. score={m['score']:.4f} r={m.get('recycles',8)} diff={m['ptm_diff']}", flush=True)
    
    pool = precise + top50[TOP_K_PRECISE:]
    pool.sort(key=lambda x: x["score"], reverse=True)
    final6 = pool[:6]
    
    import csv
    csv_path = os.path.join(WORK, "submission_r21.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for i, c in enumerate(final6):
            w.writerow(["SnowFold", i+1, c["seq"]])
            print(f"  {i+1}. score={c['score']:.4f} pTM={c['ptm']:.4f}", flush=True)
    
    json.dump(final6, open(os.path.join(WORK, "final_6_r21.json"), "w"), indent=2)
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
