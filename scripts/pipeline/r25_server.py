#!/usr/bin/env python3
"""R25 Server - 链式自动执行 (R24 完成后由 watcher 启动)

Strategy: 用 R24 Top 6 (预期 0.9447 series) 作父代
- 温度: [0.05, 0.1, 0.2, 0.4, 0.6] (中低温度, 精细优化)
- 候选: 1000/父代 × 6 = 6000
- Phase 2: Top 20 用 r=20 重算
"""
import os, sys, json, time, glob, subprocess, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r25"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "pdbs_precise", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

NUM_SEQ_PER_TEMP = 200
TEMPS = [0.05, 0.1, 0.2, 0.4, 0.6]  # 中低温度精细优化
FIXED = [1, 65, 66, 67, 96, 222]
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
TOP_K = 50
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
    return {"ptm": round(ptm,4),"plddt": round(gp,4),"chromo": round(cp,4),
            "score": round(score,4),"passes": ptm>0.60 and gp>0.60 and cp>0.55}

def predict_and_save_pdb(seq, pdb_path, recycles=RECYCLES_SCREEN):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    metrics = {"ptm":round(ptm,4),"plddt":round(gp,4),"chromo":round(cp,4),
               "score":round(score,4),"passes":ptm>0.60 and gp>0.60 and cp>0.55}
    positions = out.positions[-1][0].cpu().numpy()
    with open(pdb_path, "w") as f:
        f.write("REMARK  R25\n"); aidx=1
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
            m = predict(s, recycles)
        except Exception:
            torch.cuda.empty_cache(); continue
        if not m["passes"]: continue
        m["name"] = e["name"][:50]; m["seq"] = s; m["length"] = len(s)
        m["parent"] = parent; m["recycles"] = recycles
        results.append(m)
        torch.cuda.empty_cache()
        if (i+1) % 100 == 0:
            print(f"    [{time.strftime('%H:%M:%S')}] {i+1}/{n_total}, {len(results)} passed", flush=True)
    results.sort(key=lambda x: x["score"], reverse=True)
    return results

if __name__ == "__main__":
    t0 = time.time()
    print(f"R25 Server Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)

    # Load R24 Top 6 (fallback to progress.json if final not ready)
    r24_final = "/root/autodl-tmp/r24/final_6_r24.json"
    r24_prog = "/root/autodl-tmp/r24/results/all_passed.json"
    if os.path.exists(r24_final):
        parents = json.load(open(r24_final))
        print(f"Loaded R24 final_6 ({len(parents)} parents)", flush=True)
    elif os.path.exists(r24_prog):
        parents = json.load(open(r24_prog))[:6]
        print(f"WARN: R24 final not found, using progress.json top 6", flush=True)
    else:
        print("ERROR: No R24 results found! Falling back to R22 Top 6", flush=True)
        parents = json.load(open("/root/autodl-tmp/r22/final_6_r22.json"))

    print(f"\nR25 父代 (from R24):", flush=True)
    for i, p in enumerate(parents[:6]):
        print(f"  P{i+1}: score={p['score']:.4f}", flush=True)

    all_passed = []
    for i, parent in enumerate(parents[:6]):
        ps = parent["seq"]
        pn = f"r25_p{i+1}"
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print(f"\n[Parent {i+1}/6] ESMFold r=8 -> PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb, RECYCLES_SCREEN)
            print("done", flush=True)
        ffs = run_mpnn(pdb, pn)
        if not ffs:
            print(f"  MPNN failed for {pn}", flush=True); continue
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
    for j, c in enumerate(all_passed[:10]):
        print(f"  {j+1:2d}. parent={c['parent']} score={c['score']:.4f} pTM={c['ptm']:.4f}", flush=True)

    # Phase 2: r=20 recount top 20
    print(f"\n=== Phase 2: r=20 recount top 20 ===", flush=True)
    top20 = all_passed[:20]
    precise = []
    for i, c in enumerate(top20):
        nm = f"top{i+1:02d}_{c['parent'][:10]}_s{c['score']:.4f}"
        pdb_p = os.path.join(WORK, "pdbs_precise", f"{nm}.pdb")
        m = predict_and_save_pdb(c["seq"], pdb_p, RECYCLES_PRECISE)
        m["name"]=c["name"]; m["seq"]=c["seq"]; m["parent"]=c["parent"]
        ptm_diff = abs(m["ptm"]-c["ptm"])
        chosen = m if ptm_diff<=0.05 or m["score"]>=c["score"] else c
        chosen["recycles"]=RECYCLES_PRECISE; chosen["ptm_diff"]=round(ptm_diff,4)
        precise.append(chosen)
        print(f"  {i+1:2d}. score={chosen['score']:.4f} pTM={chosen['ptm']:.4f} (diff={ptm_diff:.4f})", flush=True)

    pool = precise + all_passed[20:TOP_K]
    pool.sort(key=lambda x: x["score"], reverse=True)
    final6 = pool[:6]
    import csv
    csv_path = os.path.join(WORK, "submission_r25.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Team_Name","Seq_ID","Sequence"])
        for j, c in enumerate(final6):
            w.writerow(["SnowFold", j+1, c["seq"]])
    json.dump(final6, open(os.path.join(WORK,"final_6_r25.json"),"w"), indent=2)
    print(f"\n=== R25 FINAL Top 6 ===", flush=True)
    for j, c in enumerate(final6):
        print(f"  {j+1}. score={c['score']:.4f} pTM={c['ptm']:.4f} chromo={c['chromo']:.3f}", flush=True)
    r24_top1 = parents[0]['score']
    print(f"\nR25 Top 1: {final6[0]['score']:.4f} | R24 Top 1: {r24_top1:.4f} | Delta: {final6[0]['score']-r24_top1:+.4f}", flush=True)
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)