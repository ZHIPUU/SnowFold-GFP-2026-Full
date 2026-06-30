#!/usr/bin/env python3
"""R20 Pipeline - A800 strict rules-compliant batch + high-prec r=20 recount

Strict project rules:
- 5.2: ~1000 sequences/batch, keep Top 50
- 5.2: 4 temperatures (0.1, 0.2, 0.5, 1.0), >=200/temp
- 5.3: Top 20 MUST high-precision recount (num_recycles>=20)
- 5.3: If pTM diff > 0.05, take min
- Score: 0.4*pTM + 0.3*pLDDT/100 + 0.3*chromo_pLDDT/100
"""
import os, sys, json, time, glob, subprocess, copy, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r20"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "pdbs_precise", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

# === Strict rules-compliant config ===
NUM_SEQ_PER_TEMP = 250
TEMPS = [0.1, 0.2, 0.5, 1.0]      # 4 档温度
FIXED = [1, 65, 66, 67, 96, 222]  # position 1 (M) + 5 chromophore核心
RECYCLES_SCREEN = 8                # 筛选 r=8
RECYCLES_PRECISE = 20              # 高精度 r=20 (项目规则 5.3)
TOP_K_SCREEN = 50                  # 规则 5.2: 保留 Top 50
TOP_K_PRECISE = 20                 # 规则 5.3: Top 20 高精度重算
BATCH = 25

# === R19 Top 6 parents ===
R19_TOP6 = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))

# === Load ESMFold ===
print("Loading ESMFold...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

def predict(seq, recycles=RECYCLES_SCREEN):
    """Returns {ptm, plddt(0-1), chromo(0-1), score (0-1 scale)}"""
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()  # 0-1
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean()); cp = float(plddt[57:72].mean())
    # 项目规则公式: 0.4*pTM + 0.3*pLDDT + 0.3*chromo (pLDDT 0-1)
    score = 0.40*ptm + 0.30*gp + 0.30*cp
    return {
        "ptm": round(ptm, 4), 
        "plddt": round(gp, 4), 
        "chromo": round(cp, 4),
        "score": round(score, 4),
        "passes": ptm > 0.60 and gp > 0.60 and cp > 0.55
    }

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
        f.write("REMARK  R20\n"); aidx=1
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

def list_fa_files(outdir):
    """Robust .fa file listing — glob has issues with names containing ',' or spaces"""
    fa_files = []
    seqs_dir = os.path.join(outdir, "seqs")
    if not os.path.isdir(seqs_dir):
        return []
    # Walk manually
    for root, dirs, files in os.walk(seqs_dir):
        for f in files:
            if f.endswith(".fa"):
                fa_files.append(os.path.join(root, f))
    return sorted(fa_files)

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
    """Round 1: r=8 screening, keep Top 50"""
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
        # Free GPU memory after each prediction
        torch.cuda.empty_cache()
        if (i+1) % 25 == 0:
            print(f"    [{time.strftime('%H:%M:%S')}] screened {i+1}/{n_total}, {len(results)} passed, VRAM={torch.cuda.memory_allocated()/1024**3:.1f}GB", flush=True)
    results.sort(key=lambda x: x["score"], reverse=True)
    return results

def high_precision_recount(seq_data, name):
    """Rule 5.3: Top 20 with r=20"""
    seq = seq_data["seq"]
    pdb_precise = os.path.join(WORK, "pdbs_precise", f"{name}.pdb")
    if os.path.isfile(pdb_precise):
        # Already done, just re-evaluate
        pass
    # Predict at r=20 and save PDB
    m_precise = predict_and_save_pdb(seq, pdb_precise, RECYCLES_PRECISE)
    m_precise["name"] = name
    m_precise["seq"] = seq
    m_precise["length"] = len(seq)
    
    # Compare with r=8 result
    ptm_diff = abs(m_precise["ptm"] - seq_data["ptm"])
    
    # Rule 5.3: if pTM diff > 0.05, take the lower
    if ptm_diff > 0.05:
        # Take min score
        if m_precise["score"] < seq_data["score"]:
            chosen = m_precise
            chosen["note"] = f"r=20 chosen (r=8 was higher, diff={ptm_diff:.3f})"
        else:
            chosen = seq_data
            chosen["note"] = f"r=8 kept (r=20 was lower, diff={ptm_diff:.3f})"
    else:
        # Use r=20 as authoritative
        chosen = m_precise
        chosen["note"] = f"r=20 (diff={ptm_diff:.3f} < 0.05)"
    
    chosen["recycles"] = RECYCLES_PRECISE
    chosen["ptm_diff"] = round(ptm_diff, 4)
    return chosen

if __name__ == "__main__":
    t0 = time.time()
    print(f"="*70, flush=True)
    print(f"R20 Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"Config: {len(TEMPS)} temps={TEMPS}, {NUM_SEQ_PER_TEMP}/temp = {NUM_SEQ_PER_TEMP*len(TEMPS)}/parent", flush=True)
    print(f"Rules: Top {TOP_K_SCREEN} after r={RECYCLES_SCREEN} screen, Top {TOP_K_PRECISE} with r={RECYCLES_PRECISE}", flush=True)
    print(f"="*70, flush=True)
    
    all_results_screen = []  # 全局 Top 50 候选 (跨父代)
    
    # === Phase 1: ESMFold PDB + MPNN ===
    for parent in R19_TOP6:
        pn = parent["name"]; ps = parent["seq"]
        print(f"\n[{time.strftime('%H:%M:%S')}] Parent: {pn} ({len(ps)}aa)", flush=True)
        
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print(f"  ESMFold r=8 → PDB...", end=" ", flush=True)
            predict_and_save_pdb(ps, pdb, RECYCLES_SCREEN)
            print("done", flush=True)
        
        ffs = run_mpnn(pdb, pn)
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"] != ps and x["seq"].startswith("M")]
        print(f"  filtered: {len(filt)} (excl parent & non-M)", flush=True)
        if not filt: continue
        
        r = screen_candidates(filt, pn)
        all_results_screen.extend(r)
        print(f"  Parent {pn}: {len(r)} passed screen", flush=True)
    
    # === Phase 2: Cross-parent Top 50 ===
    print(f"\n[{time.strftime('%H:%M:%S')}] === Top {TOP_K_SCREEN} after screening ===", flush=True)
    all_results_screen.sort(key=lambda x: x["score"], reverse=True)
    top50 = all_results_screen[:TOP_K_SCREEN]
    print(f"Global Top {TOP_K_SCREEN}:", flush=True)
    for i, c in enumerate(top50[:10]):
        print(f"  {i+1}. {c['parent'][:20]:20s} score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f}", flush=True)
    if len(top50) > 10:
        print(f"  ... and {len(top50)-10} more", flush=True)
    
    json.dump(top50, open(os.path.join(WORK, "results", "top50_screen.json"), "w"), indent=2)
    
    # === Phase 3: Top 20 High-Precision Recount (r=20) ===
    print(f"\n[{time.strftime('%H:%M:%S')}] === Top {TOP_K_PRECISE} High-Precision Recount (r={RECYCLES_PRECISE}) ===", flush=True)
    top20 = top50[:TOP_K_PRECISE]
    precise_results = []
    for i, c in enumerate(top20):
        nm = f"top{i+1:02d}_{c['parent'][:15]}_{c['score']:.4f}"
        print(f"  [{i+1}/{len(top20)}] {nm}...", end=" ", flush=True)
        m = high_precision_recount(c, nm)
        precise_results.append(m)
        print(f"score={m['score']:.4f} (pTM={m['ptm']:.4f}, diff={m['ptm_diff']:.4f}, note={m['note']})", flush=True)
    
    json.dump(precise_results, open(os.path.join(WORK, "results", "top20_precise.json"), "w"), indent=2)
    
    # === Phase 4: Final Top 6 ===
    print(f"\n[{time.strftime('%H:%M:%S')}] === FINAL Top 6 (after r=20 recount) ===", flush=True)
    # Use precise_results for top20, fall back to top50[20:] for 21-50
    pool = precise_results + top50[TOP_K_PRECISE:]
    pool.sort(key=lambda x: x["score"], reverse=True)
    final6 = pool[:6]
    
    import csv
    csv_path = os.path.join(WORK, "submission_r20.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for i, c in enumerate(final6):
            w.writerow(["SnowFold", i+1, c["seq"]])
            print(f"  {i+1}. parent={c['parent'][:20]:20s} score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f} r={c.get('recycles',8)}", flush=True)
    
    json.dump(final6, open(os.path.join(WORK, "final_6_r20.json"), "w"), indent=2)
    
    # === Compare with R19 ===
    r19_top = R19_TOP6[0]["score"]
    r20_top = final6[0]["score"]
    print(f"\n[{time.strftime('%H:%M:%S')}] === R19 vs R20 ===", flush=True)
    print(f"  R19 Top 1: {r19_top:.4f}", flush=True)
    print(f"  R20 Top 1: {r20_top:.4f}", flush=True)
    delta = r20_top - r19_top
    print(f"  Δ: {delta:+.4f} ({(delta/r19_top)*100:+.2f}%)", flush=True)
    print(f"  Status: {'✅ BREAKTHROUGH' if delta > 0 else '⚠️ NO IMPROVEMENT'}", flush=True)
    
    print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print(f"Submit: {csv_path}", flush=True)
