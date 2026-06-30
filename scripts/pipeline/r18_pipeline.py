#!/usr/bin/env python3
"""R18 Pipeline - A800 full MPNN exploration + correct ESMFold scoring"""
import os, sys, json, time, glob, subprocess, copy
import numpy as np
import torch
import warnings; warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r18"
MPNN = "/root/autodl-tmp/ProteinMPNN"
NUM_SEQ = 50
BATCH = 10
TEMPS = [0.1, 0.15, 0.2]
FIXED = [65, 66, 67, 96, 222]
RECYCLES = 8

os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

# === R17 Top 6 parents (from correct scoring) ===
R17_TOP6 = json.load(open("/root/autodl-tmp/r17/final_6_r17.json"))

# === Model loading ===
print("Loading ESMFold...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

def predict(seq, recycles=RECYCLES):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {"ptm":round(ptm,4),"plddt":round(gp,2),"chromo":round(cp,2),
            "score":round(score,4),"passes":ptm>0.6 and gp>0.6 and cp>0.55}

def predict_and_save_pdb(seq, pdb_path):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    metrics = {"ptm":round(ptm,4),"plddt":round(gp,2),"chromo":round(cp,2),
               "score":round(score,4),"passes":ptm>0.6 and gp>0.6 and cp>0.55}
    
    positions = out.positions[-1][0].cpu().numpy()
    with open(pdb_path, "w") as f:
        f.write("REMARK  R18 pipeline\n"); aidx=1
        for i in range(len(seq)):
            rn = aa3.get(seq[i], 'ALA')
            for j, an in enumerate(["N","CA","C","O"]):
                x,y,z = positions[i,j]
                f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt[i]*100:6.2f}\n")
                aidx+=1
        f.write("END\n")
    return metrics

def run_mpnn(pdb_path, name):
    outdir = os.path.join(WORK, "mpnn_out", name)
    os.makedirs(outdir, exist_ok=True)
    flag = os.path.join(outdir, "done.flag")
    if os.path.isfile(flag):
        files = sorted(glob.glob(os.path.join(outdir, "seqs", "*.fa")))
        print(f"  MPNN cached: {len(files)} files", flush=True); return files
    
    fixed = os.path.join(outdir, "fixed.jsonl")
    pdb_key = os.path.basename(pdb_path).replace(".pdb", "")
    with open(fixed, "w") as f:
        f.write(json.dumps({pdb_key: {"A": FIXED}}) + "\n")
    
    cmd = [sys.executable, os.path.join(MPNN, "protein_mpnn_run.py"),
           "--pdb_path", pdb_path, "--pdb_path_chains", "A",
           "--path_to_model_weights", os.path.join(MPNN, "vanilla_model_weights"),
           "--fixed_positions_jsonl", fixed,
           "--out_folder", outdir,
           "--num_seq_per_target", str(NUM_SEQ),
           "--batch_size", str(BATCH),
           "--sampling_temp", " ".join(str(t) for t in TEMPS),
           "--seed", "42", "--suppress_print", "1"]
    print(f"  MPNN {name}...", end="", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    files = sorted(glob.glob(os.path.join(outdir, "seqs", "*.fa")))
    if not files:
        # Try running without capture_output to see full error
        print(f" no output files, retrying...", end="", flush=True)
        r2 = subprocess.run(cmd, timeout=600)
        files = sorted(glob.glob(os.path.join(outdir, "seqs", "*.fa")))
        if not files:
            print(f" still failed (rc={r.returncode})", flush=True)
            return []
    print(f" {len(files)} files", flush=True)
    if files: open(flag, "w").close()
    return files

def parse_fa(paths):
    seqs = []
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as f:
            n,s = "",""
            for l in f:
                l=l.strip()
                if l.startswith(">"):
                    if s: seqs.append({"name":n,"seq":s})
                    n,s = l[1:],""
                else: s+=l
            if s: seqs.append({"name":n,"seq":s})
    return seqs

def evaluate(seqs, parent):
    results = []
    for i,e in enumerate(seqs):
        s,n = e["seq"], e.get("name",f"{parent}_{i}")
        if not s.startswith("M") or len(s)<220 or len(s)>250: continue
        try:
            m = predict(s)
        except Exception as e:
            torch.cuda.empty_cache(); continue
        if not m["passes"]: continue
        m["name"] = n; m["seq"] = s; m["length"] = len(s); m["parent"] = parent
        results.append(m)
        if (i+1)%50==0: print(f"  {i+1}/{len(seqs)}, {len(results)} passed", flush=True)
    results.sort(key=lambda x:x["score"], reverse=True)
    n=len(results)
    print(f"  {parent}: {n}/{len(seqs)} passed, Top={results[0]['score']:.4f}" if n else f"  {parent}: 0 passed", flush=True)
    return results

if __name__ == "__main__":
    t0 = time.time()
    print(f"R18 Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    
    all_results = []
    for parent in R17_TOP6:
        pn = parent["name"]; ps = parent["seq"]
        print(f"\n--- {pn} ({len(ps)}aa) ---", flush=True)
        
        # Step 1: ESMFold → PDB
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print("  ESMFold+PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb)
            print("done", flush=True)
        
        # Step 2: MPNN
        ffs = run_mpnn(pdb, pn)
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"]!=ps and x["seq"].startswith("M")]
        print(f"  {len(filt)} to evaluate", flush=True)
        if not filt: continue
        
        # Step 3: Evaluate
        r = evaluate(filt, pn)
        all_results.extend(r)
        if r:
            json.dump(r[:30], open(os.path.join(WORK,"results",f"{pn}_top.json"),"w"), indent=2)
    
    all_results.sort(key=lambda x:x["score"], reverse=True)
    json.dump(all_results, open(os.path.join(WORK,"results","all.json"),"w"), indent=2)
    
    print(f"\n=== FINAL Top 6 ===", flush=True)
    print(f"Total passed: {len(all_results)}", flush=True)
    top6 = all_results[:6] if all_results else R17_TOP6[:6]
    
    import csv
    csv_path = os.path.join(WORK, "submission_r18.csv")
    with open(csv_path,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["Team_Name","Seq_ID","Sequence"])
        for i,c in enumerate(top6):
            w.writerow(["SnowFold",i+1,c["seq"]])
            print(f"  {i+1}. {c['name'][:28]:<30s} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.2f} chromo={c['chromo']:.2f} score={c['score']:.4f}", flush=True)
    
    json.dump(top6, open(os.path.join(WORK,"final_6_r18.json"),"w"), indent=2)
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"Submit: {csv_path}", flush=True)
