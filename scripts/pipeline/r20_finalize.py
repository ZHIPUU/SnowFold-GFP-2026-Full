#!/usr/bin/env python3
"""R20 Finalize - Use existing MPNN outputs from R20 (2 parents), no new MPNN

Evaluates 2000 candidates (2 parents x 1000), selects Top 50, recounts Top 20 with r=20
"""
import os, sys, json, time, glob, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r20"
for d in ["pdbs_fin", "pdbs_precise_fin", "results_fin"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
TOP_K_PRECISE = 20

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
        f.write("REMARK  R20_FIN\n"); aidx=1
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

if __name__ == "__main__":
    t0 = time.time()
    print(f"R20 Finalize Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    
    # Get all MPNN outputs from R20
    mpnn_dirs = sorted(glob.glob(os.path.join(WORK, "mpnn_out", "*")))
    print(f"Found {len(mpnn_dirs)} MPNN outputs:", flush=True)
    
    all_passed = []
    
    # Load R19 Top 6 parents
    R19_TOP6 = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))
    R19_BY_NAME = {p["name"]: p for p in R19_TOP6}
    
    for mpnn_dir in mpnn_dirs:
        parent_name = os.path.basename(mpnn_dir)
        parent = R19_BY_NAME.get(parent_name)
        if parent is None:
            print(f"  WARN: No parent {parent_name} found in R19", flush=True)
            continue
        ps = parent["seq"]
        
        # Parse .fa files
        fa_files = []
        for root, dirs, files in os.walk(os.path.join(mpnn_dir, "seqs")):
            for f in files:
                if f.endswith(".fa"):
                    fa_files.append(os.path.join(root, f))
        all_seqs = parse_fa(fa_files)
        filt = [x for x in all_seqs if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"\n[{time.strftime('%H:%M:%S')}] {parent_name[:40]}: {len(filt)} candidates", flush=True)
        
        # ESMFold evaluate
        passed = []
        for i, e in enumerate(filt):
            s = e["seq"]
            try:
                m = predict(s, RECYCLES_SCREEN)
            except Exception as ex:
                torch.cuda.empty_cache()
                continue
            if not m["passes"]:
                continue
            m["name"] = e["name"][:60]
            m["seq"] = s
            m["length"] = len(s)
            m["parent"] = parent_name
            m["recycles"] = RECYCLES_SCREEN
            passed.append(m)
            torch.cuda.empty_cache()
            if (i+1) % 100 == 0:
                print(f"    [{time.strftime('%H:%M:%S')}] {i+1}/{len(filt)}, {len(passed)} passed, VRAM={torch.cuda.memory_allocated()/1024**3:.1f}GB", flush=True)
        
        passed.sort(key=lambda x: x["score"], reverse=True)
        print(f"  Total passed: {len(passed)}/{len(filt)} ({100*len(passed)/max(1,len(filt)):.1f}%)", flush=True)
        all_passed.extend(passed)
    
    # Global sort
    all_passed.sort(key=lambda x: x["score"], reverse=True)
    json.dump(all_passed, open(os.path.join(WORK, "results_fin", "all_passed.json"), "w"), indent=2)
    
    # Top 50
    top50 = all_passed[:50]
    print(f"\n=== Top 50 (after r={RECYCLES_SCREEN} screen) ===", flush=True)
    for i, c in enumerate(top50):
        print(f"  {i+1:2d}. parent={c['parent'][:25]:25s} score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f}", flush=True)
    
    json.dump(top50, open(os.path.join(WORK, "results_fin", "top50.json"), "w"), indent=2)
    
    # Top 20 high-precision recount (rule 5.3)
    print(f"\n=== Top {TOP_K_PRECISE} high-precision (r={RECYCLES_PRECISE}) ===", flush=True)
    top20 = top50[:TOP_K_PRECISE]
    precise = []
    for i, c in enumerate(top20):
        nm = f"top{i+1:02d}_p{c['parent'][:10]}_s{c['score']:.4f}"
        pdb_precise = os.path.join(WORK, "pdbs_precise_fin", f"{nm}.pdb")
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
        print(f"  {i+1:2d}. r={RECYCLES_PRECISE}: score={chosen['score']:.4f} pTM={chosen['ptm']:.4f} (diff={ptm_diff:.4f}, {chosen['note']})", flush=True)
    
    # Final Top 6
    pool = precise + top50[TOP_K_PRECISE:]
    pool.sort(key=lambda x: x["score"], reverse=True)
    final6 = pool[:6]
    
    import csv
    csv_path = os.path.join(WORK, "submission_r20.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for i, c in enumerate(final6):
            w.writerow(["SnowFold", i+1, c["seq"]])
            print(f"  FINAL {i+1}: score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f} r={c.get('recycles',8)}", flush=True)
    
    json.dump(final6, open(os.path.join(WORK, "final_6_r20.json"), "w"), indent=2)
    
    # Compare with R19
    r19_top = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))[0]["score"]
    r20_top = final6[0]["score"]
    print(f"\n=== R19 vs R20 ===", flush=True)
    print(f"  R19 Top 1: {r19_top:.4f}", flush=True)
    print(f"  R20 Top 1: {r20_top:.4f}", flush=True)
    delta = r20_top - r19_top
    print(f"  Δ: {delta:+.4f} ({(delta/r19_top)*100:+.2f}%)", flush=True)
    print(f"  Status: {'✅ BREAKTHROUGH' if delta > 0 else '⚠️ NO IMPROVEMENT'}", flush=True)
    
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"Submit: {csv_path}", flush=True)
