#!/usr/bin/env python3
"""R27 Multi-Direction Divergence Exploration (~2h local)

5 directions, ~100-200 candidates each, total ~1200 candidates:
  A: 极低温微调 [0.01,0.05] × R25 Top2, 100/parent = 400
  B: 超高温跳跃 [1.0,1.5] × R25 Top2, 100/parent = 400  
  C: 极少固定 [66,67,222] only × R25 Top1, T=0.1, 200
  D: 超宽固定 [1,2,65,66,67,96,203,222] × R25 Top1, T=0.2, 200
  E: 全谱温度 [0.01..1.5] × R25 Top1, 20/temp × 8 = 160
"""
import os, sys, json, time, subprocess, warnings, csv, hashlib
import torch, numpy as np
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "D:/workspace/round27_diverge"
MPNN_DIR = "C:/proteinmpnn_r10"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs","mpnn_out","results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

RECYCLES = 8
BATCH = 4
aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY","ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

print("Loading ESMFold chunk=64...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(64); model.eval()
print("Loaded.", flush=True)

r25 = json.load(open(r'D:\workspace\round25\final_6_r25.json'))

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

def save_pdb(seq, path):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k:v.cuda() for k,v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES)
    plddt = out.plddt[0,:,1].cpu().numpy()
    pos = out.positions[-1][0].cpu().numpy()
    with open(path,"w") as f:
        f.write("REMARK R27\n"); aidx=1
        for i,a in enumerate(seq):
            rn=aa3.get(a,"ALA")
            for j,an in enumerate(["N","CA","C","O"]):
                x,y,z=pos[i,j]
                f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt[i]*100:6.2f}\n")
                aidx+=1
        f.write("END\n")

def list_fa(outdir):
    seqs_dir=os.path.join(outdir,"seqs"); res=[]
    if not os.path.isdir(seqs_dir): return []
    for root,dirs,files in os.walk(seqs_dir):
        for f in files:
            if f.endswith(".fa"): res.append(os.path.join(root,f))
    return sorted(res)

def parse_fa(paths):
    out=[]
    for p in paths:
        n=""; s=""
        with open(p,encoding="utf-8",errors="replace") as f:
            for line in f:
                line=line.strip()
                if line.startswith(">"):
                    if s: out.append({"name":n,"seq":s})
                    n=line[1:]; s=""
                elif line: s+=line
            if s: out.append({"name":n,"seq":s})
    return out

def run_mpnn(pdb, name, temps, num_seq, fixed):
    pdb_mpnn = pdb.replace("\\","/")
    outdir = os.path.join(WORK, "mpnn_out", name)
    os.makedirs(outdir, exist_ok=True)
    fixed_json = os.path.join(outdir, "fixed.jsonl")
    key = os.path.basename(pdb).replace(".pdb","")
    with open(fixed_json,"w") as f:
        f.write(json.dumps({key:{"A":fixed}})+"\n")
    cmd = [sys.executable, os.path.join(MPNN_DIR,"protein_mpnn_run.py"),
           "--pdb_path", pdb_mpnn, "--pdb_path_chains", "A",
           "--path_to_model_weights", os.path.join(MPNN_DIR,"vanilla_model_weights"),
           "--fixed_positions_jsonl", fixed_json, "--out_folder", outdir,
           "--num_seq_per_target", str(num_seq), "--batch_size", str(BATCH),
           "--sampling_temp", " ".join(map(str,temps)),
           "--seed", str(np.random.randint(0,99999)), "--suppress_print", "1"]
    print(f"  [MPNN] {name} temps={temps} n={num_seq} fixed={fixed}...", end="", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    files = list_fa(outdir)
    print(f" {len(files)} files", flush=True)
    if not files:
        print(f"  STDERR: {r.stderr[-300:]}", flush=True)
    return files

def screen(seqs, direction, parent_label):
    passed = []
    n = len(seqs)
    print(f"  [screen] {n} candidates @ r={RECYCLES} ({direction})", flush=True)
    for j, e in enumerate(seqs):
        s = e["seq"]
        if not s.startswith("M") or len(s)<220 or len(s)>250: continue
        try: m = predict(s)
        except Exception:
            torch.cuda.empty_cache(); continue
        if not m["passes"]: continue
        m["seq"]=s; m["name"]=e["name"][:50]; m["parent"]=parent_label
        m["direction"]=direction; m["length"]=len(s); m["recycles"]=RECYCLES
        passed.append(m)
        torch.cuda.empty_cache()
        if (j+1)%50==0:
            print(f"    [{time.strftime('%H:%M:%S')}] {j+1}/{n}, {len(passed)} passed", flush=True)
    passed.sort(key=lambda x:x["score"], reverse=True)
    return passed

if __name__ == "__main__":
    t0 = time.time()
    print(f"R27 Divergence Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    all_passed = []

    # === Direction A: 极低温微调 ===
    print(f"\n{'='*60}\nDirection A: 极低温微调 [0.01, 0.05]\n{'='*60}", flush=True)
    for i in range(2):
        parent = r25[i]
        pn = f"A_p{i+1}"
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print(f"  PDB for {pn}...", end="", flush=True)
            save_pdb(parent["seq"], pdb)
            print(" done", flush=True)
        files = run_mpnn(pdb, pn, [0.01, 0.05], 50, [1,65,66,67,96,222])
        if files:
            seqs = parse_fa(files)
            filt = [x for x in seqs if x["seq"]!=parent["seq"] and x["seq"].startswith("M")]
            r = screen(filt, "A_ultra_low", pn)
            print(f"  {pn}: {len(r)}/{len(filt)} passed, top={r[0]['score']:.4f}" if r else f"  {pn}: 0 passed", flush=True)
            all_passed.extend(r)

    # === Direction B: 超高温跳跃 ===
    print(f"\n{'='*60}\nDirection B: 超高温跳跃 [1.0, 1.5]\n{'='*60}", flush=True)
    for i in range(2):
        parent = r25[i]
        pn = f"B_p{i+1}"
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            save_pdb(parent["seq"], pdb)
        files = run_mpnn(pdb, pn, [1.0, 1.5], 50, [1,65,66,67,96,222])
        if files:
            seqs = parse_fa(files)
            filt = [x for x in seqs if x["seq"]!=parent["seq"] and x["seq"].startswith("M")]
            r = screen(filt, "B_ultra_high", pn)
            print(f"  {pn}: {len(r)}/{len(filt)} passed, top={r[0]['score']:.4f}" if r else f"  {pn}: 0 passed", flush=True)
            all_passed.extend(r)

    # === Direction C: 极少固定 (只固定 chromophore 核心 3 位) ===
    print(f"\n{'='*60}\nDirection C: 极少固定 [66,67,222] only\n{'='*60}", flush=True)
    parent = r25[0]
    pn = "C_minimal_fixed"
    pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
    if not os.path.isfile(pdb):
        save_pdb(parent["seq"], pdb)
    files = run_mpnn(pdb, pn, [0.1], 50, [66, 67, 222])
    if files:
        seqs = parse_fa(files)
        filt = [x for x in seqs if x["seq"].startswith("M")]
        r = screen(filt, "C_minimal_fixed", pn)
        print(f"  {pn}: {len(r)}/{len(filt)} passed, top={r[0]['score']:.4f}" if r else f"  {pn}: 0 passed", flush=True)
        all_passed.extend(r)

    # === Direction D: 超宽固定 (固定更多关键位点) ===
    print(f"\n{'='*60}\nDirection D: 超宽固定 [1,2,65,66,67,96,203,222]\n{'='*60}", flush=True)
    parent = r25[0]
    pn = "D_wide_fixed"
    pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
    if not os.path.isfile(pdb):
        save_pdb(parent["seq"], pdb)
    files = run_mpnn(pdb, pn, [0.2], 50, [1, 2, 65, 66, 67, 96, 203, 222])
    if files:
        seqs = parse_fa(files)
        filt = [x for x in seqs if x["seq"]!=parent["seq"] and x["seq"].startswith("M")]
        r = screen(filt, "D_wide_fixed", pn)
        print(f"  {pn}: {len(r)}/{len(filt)} passed, top={r[0]['score']:.4f}" if r else f"  {pn}: 0 passed", flush=True)
        all_passed.extend(r)

    # === Direction E: 全谱温度 ===
    print(f"\n{'='*60}\nDirection E: 全谱温度 [0.01..1.5]\n{'='*60}", flush=True)
    parent = r25[0]
    pn = "E_full_spectrum"
    pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
    if not os.path.isfile(pdb):
        save_pdb(parent["seq"], pdb)
    all_temps = [0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5]
    files = run_mpnn(pdb, pn, all_temps, 10, [1, 65, 66, 67, 96, 222])
    if files:
        seqs = parse_fa(files)
        filt = [x for x in seqs if x["seq"]!=parent["seq"] and x["seq"].startswith("M")]
        r = screen(filt, "E_full_spectrum", pn)
        print(f"  {pn}: {len(r)}/{len(filt)} passed, top={r[0]['score']:.4f}" if r else f"  {pn}: 0 passed", flush=True)
        all_passed.extend(r)

    # === Summary ===
    all_passed.sort(key=lambda x: x["score"], reverse=True)
    json.dump(all_passed, open(os.path.join(WORK, "all_passed.json"), "w"), indent=2)

    print(f"\n{'='*60}\nSUMMARY: {len(all_passed)} total passed\n{'='*60}", flush=True)
    by_dir = {}
    for c in all_passed:
        d = c.get("direction","?")
        by_dir.setdefault(d, []).append(c)
    for d, items in sorted(by_dir.items()):
        print(f"  {d}: {len(items)} passed, top={items[0]['score']:.4f}", flush=True)

    print(f"\nTop 15 overall:", flush=True)
    for i, c in enumerate(all_passed[:15]):
        print(f"  {i+1:2d}. score={c['score']:.4f} pTM={c['ptm']:.4f} chromo={c['chromo']:.3f} dir={c['direction']} parent={c['parent']}", flush=True)

    final6 = all_passed[:6]
    json.dump(final6, open(os.path.join(WORK, "final_6_r27.json"), "w"), indent=2)
    with open(os.path.join(WORK, "submission_r27.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Team_Name","Seq_ID","Sequence"])
        for i, c in enumerate(final6):
            w.writerow(["SnowFold", i+1, c["seq"]])

    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    if all_passed:
        print(f"R27 Top 1: {all_passed[0]['score']:.4f} | R25 Top 1: 0.9477 | Delta: {all_passed[0]['score']-0.9477:+.4f}", flush=True)