#!/usr/bin/env python3
"""R23 Server - DIVERGENT EXPLORATION on A800

基于 R20 Top 3 父代, 5 温度 × 150 候选 = 750 候选/父代 × 3 = 2250 候选
高温度比例 (T=0.3/0.5/0.7/1.0) 更大 → 多样性

R22 已经达到 0.9430. R23 是发散性:
- 用 R20 (而非 R22) Top 3 父代,避免重复
- 高温度更激进,跳出 R22 的局部最优
"""
import os, sys, json, time, glob, subprocess, copy, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r23"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs","mpnn_out","results"]:
    os.makedirs(os.path.join(WORK,d), exist_ok=True)

NUM_SEQ_PER_TEMP = 150
TEMPS = [0.1, 0.3, 0.5, 0.7, 1.0]
FIXED = [1, 65, 66, 67, 96, 222]
RECYCLES_SCREEN = 8
BATCH = 25

print("Loading ESMFold...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

def predict(seq):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES_SCREEN)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {"ptm": round(ptm,4),"plddt": round(gp,4),"chromo": round(cp,4),
            "score": round(score,4),"passes": ptm>0.60 and gp>0.60 and cp>0.55}

def predict_and_save_pdb(seq, pdb_path):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES_SCREEN)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    metrics = {"ptm":round(ptm,4),"plddt":round(gp,4),"chromo":round(cp,4),
               "score":round(score,4),"passes":ptm>0.60 and gp>0.60 and cp>0.55}
    positions = out.positions[-1][0].cpu().numpy()
    with open(pdb_path, "w") as f:
        f.write("REMARK  R23\n"); aidx=1
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
    if not os.path.isdir(seqs_dir): return []
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
           "--fixed_positions_jsonl", fixed_json, "--out_folder", outdir,
           "--num_seq_per_target", str(NUM_SEQ_PER_TEMP), "--batch_size", str(BATCH),
           "--sampling_temp", " ".join(str(t) for t in TEMPS),
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

def screen(seqs, parent, recycles=RECYCLES_SCREEN):
    results = []
    n_total = len(seqs)
    print(f"  [screen] {n_total} candidates @ r={recycles}", flush=True)
    for i, e in enumerate(seqs):
        s = e["seq"]
        if not s.startswith("M") or len(s) < 220 or len(s) > 250: continue
        try:
            m = predict(s)
        except Exception as ex:
            torch.cuda.empty_cache()
            continue
        if not m["passes"]: continue
        m["name"] = e["name"][:50]
        m["seq"] = s
        m["length"] = len(s)
        m["parent"] = parent
        m["recycles"] = recycles
        results.append(m)
        torch.cuda.empty_cache()
        if (i+1) % 100 == 0:
            print(f"    [{time.strftime('%H:%M:%S')}] {i+1}/{n_total}, {len(results)} passed", flush=True)
    results.sort(key=lambda x: x["score"], reverse=True)
    return results

if __name__ == "__main__":
    t0 = time.time()
    print(f"R23 Server Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"Strategy: 3 R20 Top 父代 × 5 温度 × 150 候选 = 2250 候选", flush=True)
    print(f"高温度比例更大: {TEMPS}", flush=True)
    print(f"Fixed: {FIXED}", flush=True)

    r20_top6 = json.load(open("/root/autodl-tmp/r22/results/r20_top6.json"))
    r20_top3 = r20_top6[:3]
    print(f"\nR20 Top 3 父代:", flush=True)
    for i, p in enumerate(r20_top3):
        print(f"  P{i+1}: score={p['score']:.4f}", flush=True)

    all_passed = []
    for i, parent in enumerate(r20_top3):
        ps = parent["seq"]
        pn = f"r23_p{i+1}"

        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print(f"\n[Parent {i+1}/3] ESMFold r={RECYCLES_SCREEN} -> PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb)
            print("done", flush=True)

        ffs = run_mpnn(pdb, pn)
        if not ffs:
            print(f"  MPNN failed for {pn}", flush=True)
            continue
        all_seqs = parse_fa(ffs)
        filt = [x for x in all_seqs if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"  filtered: {len(filt)}", flush=True)
        if not filt: continue

        r = screen(filt, pn, RECYCLES_SCREEN)
        print(f"  Parent {pn}: {len(r)}/{len(filt)} passed", flush=True)
        all_passed.extend(r)
        json.dump(all_passed, open(os.path.join(WORK,"results","progress.json"),"w"), indent=2)

    all_passed.sort(key=lambda x: x["score"], reverse=True)
    json.dump(all_passed, open(os.path.join(WORK,"results","all_passed.json"),"w"), indent=2)
    print(f"\nTotal passed: {len(all_passed)}", flush=True)
    print(f"Top 10:", flush=True)
    for j, c in enumerate(all_passed[:10]):
        print(f"  {j+1:2d}. parent={c['parent']} score={c['score']:.4f} pTM={c['ptm']:.4f}", flush=True)

    final6 = all_passed[:6]
    import csv
    csv_path = os.path.join(WORK, "submission_r23.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for j, c in enumerate(final6):
            w.writerow(["SnowFold", j+1, c["seq"]])
    json.dump(final6, open(os.path.join(WORK,"final_6_r23.json"),"w"), indent=2)

    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"Submit: {csv_path}", flush=True)