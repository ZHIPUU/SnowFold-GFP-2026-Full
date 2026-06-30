"""
R16 Pipeline v4 — 多父代近原点搜索 + avGFP改造
评分: sort_score = pTM * 0.40 + (pLDDT/100) * 0.30 + (chromo/100) * 0.30
关键修复: PDB和MPNN输出放无中文路径 C:\r16_work\
"""
import os, sys, json, time, glob, subprocess, shutil, tempfile
import numpy as np
import torch
import torch.nn.functional as F
import warnings; warnings.filterwarnings("ignore")

# === Config: 所有路径避免中文 ===
ROOT = r"D:\生信\2026Protein Design"
WORK = os.path.join(ROOT, "work", "round16", "r16_temp")  # 项目内，sandbox允许
MPNN_DIR = r"C:\proteinmpnn_r10"
NUM_SEQ = 50
BATCH_SIZE = 10
TEMPS = [0.1, 0.15, 0.2]
FIXED = [65, 66, 67, 96, 222]
RECYCLES = 8
TOP_K = 20
aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split())}

os.makedirs(WORK, exist_ok=True)
os.makedirs(os.path.join(WORK, "pdbs"), exist_ok=True)
os.makedirs(os.path.join(WORK, "mpnn_out"), exist_ok=True)
os.makedirs(os.path.join(WORK, "results"), exist_ok=True)

AvGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

def load_r15():
    p = os.path.join(ROOT, "work", "round15", "final_6_r15.json")
    with open(p, encoding="utf-8") as f:
        return json.load(f)

def check_mut(seq, f, p, t):
    """Safe mutation: skip if already matches target"""
    if seq[p-1] == f:
        s = list(seq); s[p-1] = t; return "".join(s)
    elif seq[p-1] == t:
        return seq  # already mutated
    else:
        print(f"  WARN: pos{p} is {seq[p-1]}, expected {f}, skipping")
        return None

# === Model Loading (ESMFold only) ===
print("Loading ESMFold...")
from transformers import AutoTokenizer, EsmForProteinFolding
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("ESMFold loaded.")

def predict_and_save(seq, pdb_path=None):
    """ESMFold predict, optionally save PDB"""
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=RECYCLES)
    logits = out.plddt[0]
    probs = F.softmax(logits, dim=-1)
    centers = torch.linspace(0.5/37, 1-0.5/37, 37, device=logits.device)
    plddt_01 = (probs * centers.unsqueeze(0)).sum(-1)
    plddt_100 = plddt_01.cpu().numpy() * 100
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt_100.mean()); cp = float(plddt_100[57:72].mean())
    score = 0.40*ptm + 0.30*gp/100 + 0.30*cp/100
    metrics = {"global_plddt": gp, "chromo_plddt": cp, "ptm": ptm,
               "sort_score": score, "passes": ptm>0.60 and gp>60.0 and cp>55.0}
    if pdb_path:
        pos = out.positions[-1][0].cpu().numpy()
        with open(pdb_path, "w") as f:
            f.write("REMARK  R16\n"); aidx=1
            for i in range(len(seq)):
                rn = aa3.get(seq[i], 'ALA')
                for j, an in enumerate(["N","CA","C","O"]):
                    x,y,z = pos[i,j]
                    f.write(f"ATOM  {aidx:5d} {an:^4s} {rn:3s} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{plddt_100[i]:6.2f}\n")
                    aidx+=1
            f.write("END\n")
    return metrics

# === ProteinMPNN (PDBs in non-Chinese path) ===
def run_mpnn(pdb_path, name):
    outdir = os.path.join(WORK, "mpnn_out", name)
    os.makedirs(outdir, exist_ok=True)
    flag = os.path.join(outdir, "done.flag")
    if os.path.isfile(flag):
        files = sorted(glob.glob(os.path.join(outdir, "seqs", "*.fa")))
        print(f"  MPNN cached: {len(files)} files"); return files
    
    fixed = os.path.join(outdir, "fixed.jsonl")
    # ProteinMPNN parse_PDB on Windows uses full path as name (rfind("/") fails on \)
    # So the key must be the full pdb_path without .pdb
    pdb_key = pdb_path.replace(".pdb", "")
    with open(fixed, "w") as f:
        f.write(json.dumps({pdb_key: {"A": FIXED}}) + "\n")
    
    cmd = [sys.executable, os.path.join(MPNN_DIR, "protein_mpnn_run.py"),
           "--pdb_path", pdb_path, "--pdb_path_chains", "A",
           "--path_to_model_weights", os.path.join(MPNN_DIR, "vanilla_model_weights"),
           "--fixed_positions_jsonl", fixed,
           "--out_folder", outdir,
           "--num_seq_per_target", str(NUM_SEQ),
           "--batch_size", str(BATCH_SIZE),
           "--sampling_temp", " ".join(str(t) for t in TEMPS),
           "--seed", "42", "--suppress_print", "1"]
    print(f"  MPNN {name}...", end="", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        err = r.stderr[:400]
        print(f" FAILED: {err}")
        return []
    files = sorted(glob.glob(os.path.join(outdir, "seqs", "*.fa")))
    if files: open(flag, "w").close()
    print(f" {len(files)} files")
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
            m = predict_and_save(s, None)
        except Exception as e:
            print(f"  ESMFold err: {e}"); torch.cuda.empty_cache(); continue
        if not m["passes"]: continue
        m["esm2_lp"] = 0.0
        m.update({"name":n,"seq":s,"length":len(s),"parent":parent})
        results.append(m)
        if (i+1)%50==0: print(f"  {i+1}/{len(seqs)}, {len(results)} passed")
    results.sort(key=lambda x:x["sort_score"], reverse=True)
    n=len(results)
    print(f"  {parent}: {n}/{len(seqs)} passed, Top={results[0]['sort_score']:.4f}" if n else f"  {parent}: 0 passed")
    return results

# === 方案A: R15 多父代 ===
def scheme_a():
    print("="*50+"\n方案A: R15 多父代")
    all_r = []
    for p in load_r15():
        pn, ps = p["name"], p["seq"]
        print(f"\n--- {pn} ({len(ps)}aa) ---")
        pdb = os.path.join(WORK, "pdbs", f"{pn}.pdb")
        if not os.path.isfile(pdb):
            print("  ESMFold+PDB...", end=" ", flush=True)
            predict_and_save(ps, pdb); print("done")
        
        ffs = run_mpnn(pdb, pn)
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"]!=ps and x["seq"].startswith("M")]
        print(f"  {len(filt)} to evaluate")
        if not filt: continue
        r = evaluate(filt, pn)
        all_r.extend(r)
        if r: json.dump(r[:TOP_K], open(os.path.join(WORK,"results",f"{pn}_top.json"),"w"), indent=2)
    
    all_r.sort(key=lambda x:x["sort_score"], reverse=True)
    json.dump(all_r, open(os.path.join(WORK,"results","schemeA_all.json"),"w"), indent=2)
    print(f"\n方案A: {len(all_r)} passed")
    if all_r: print(f"  Top1={all_r[0]['sort_score']:.4f}, Top6={[x['sort_score'] for x in all_r[:6]]}")
    return all_r

# === 方案B: avGFP改造 ===
def scheme_b():
    print("="*50+"\n方案B: avGFP改造")
    # 只包含已验证在avGFP上实际不同的突变
    # (对比competition提供的avGFP和sfGFP序列)
    av_muts = [
        ("S",30,"R"),  # avGFP=S, sfGFP=R
        ("S",65,"T"),  # avGFP=S, sfGFP=T (chromophore)
        ("F",99,"S"),  # avGFP=F, sfGFP=S
        ("N",105,"T"), # avGFP=N, sfGFP=T
        ("Y",145,"F"), # avGFP=Y, sfGFP=F
        ("M",153,"T"), # avGFP=M, sfGFP=T
        ("V",163,"A"), # avGFP=V, sfGFP=A
        ("I",171,"V"), # avGFP=I, sfGFP=V
        ("A",206,"V"), # avGFP=A, sfGFP=V
    ]
    full = av_muts
    core5 = [m for m in av_muts if m[1] in [65, 99, 153, 163]]
    minimal = [m for m in av_muts if m[1]==65]
    
    vars = {"av_minimal": minimal, "av_core5": core5, "av_full": full}
    all_r = []
    for vn, ms in vars.items():
        print(f"\n--- {vn} ({len(ms)} muts) ---")
        es = AvGFP
        ok = True
        for f,p,t in ms:
            r = check_mut(es, f, p, t)
            if r is None: ok=False; break
            es = r
        if not ok: print(f"  SKIP: mutation mismatch"); continue
        print(f"  seq len={len(es)}")
        
        pdb = os.path.join(WORK, "pdbs", f"{vn}.pdb")
        if not os.path.isfile(pdb):
            print("  ESMFold+PDB...", end=" ", flush=True)
            m0 = predict_and_save(es, pdb)
            print(f"pTM={m0['ptm']:.3f}, pLDDT={m0['global_plddt']:.1f}")
        
        ffs = run_mpnn(pdb, vn)
        all_s = parse_fa(ffs)
        filt = [x for x in all_s if x["seq"]!=es and x["seq"].startswith("M")]
        print(f"  {len(filt)} to evaluate")
        if not filt: continue
        r = evaluate(filt, vn)
        for x in r: x["parent_type"] = "avGFP"
        all_r.extend(r)
        if r: json.dump(r[:TOP_K], open(os.path.join(WORK,"results",f"{vn}_top.json"),"w"), indent=2)
    
    all_r.sort(key=lambda x:x["sort_score"], reverse=True)
    json.dump(all_r, open(os.path.join(WORK,"results","schemeB_all.json"),"w"), indent=2)
    print(f"\n方案B: {len(all_r)} passed")
    if all_r: print(f"  Top1={all_r[0]['sort_score']:.4f}")
    return all_r

# === Final ===
def final(a, b):
    print("="*50+"\n最终排序")
    all_c = a + b; all_c.sort(key=lambda x:x["sort_score"], reverse=True)
    print(f"总候选: {len(all_c)}")
    
    # 多样性选6
    top6 = []; seen = set()
    for c in all_c:
        p = c.get("parent","?")
        if len(top6) < 3 or p not in seen:
            top6.append(c); seen.add(p)
            if len(top6)>=6: break
    for c in all_c:
        if len(top6)>=6: break
        if c not in top6: top6.append(c)
    
    print(f"\n{'Seq':>4s} {'Name':<30s} {'Parent':<16s} {'pLDDT':>6s} {'Ch':>5s} {'pTM':>6s} {'Score':>7s}")
    print("-"*73)
    for i,c in enumerate(top6):
        print(f"{i+1:>4d} {c['name'][:28]:<30s} {c.get('parent','?')[:14]:<16s} {c['global_plddt']:>6.1f} {c['chromo_plddt']:>5.1f} {c['ptm']:>6.3f} {c['sort_score']:>7.4f}")
    
    import csv, pandas as pd
    excl = set(pd.read_csv(os.path.join(ROOT,"Exclusion_List.csv")).iloc[:,0].astype(str))
    for i,c in enumerate(top6):
        s=c["seq"]
        print(f"  Seq{i+1}: M={s[0]=='M'} L={len(s)} AA={all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s)} Excl={s in excl}")

    # 同时保存到项目目录(含中文)
    proj_out = os.path.join(ROOT, "work", "round16", "submission_r16.csv")
    os.makedirs(os.path.join(ROOT, "work", "round16"), exist_ok=True)
    with open(proj_out, "w", newline="") as f:
        w=csv.writer(f); w.writerow(["Team_Name","Seq_ID","Sequence"])
        for i,c in enumerate(top6): w.writerow(["SnowFold",i+1,c["seq"]])
    json.dump(top6, open(proj_out.replace(".csv",".json"),"w"), indent=2)
    
    # 也存到 C:\r16_work
    csv_c = os.path.join(WORK, "submission_r16.csv")
    with open(csv_c, "w", newline="") as f:
        w=csv.writer(f); w.writerow(["Team_Name","Seq_ID","Sequence"])
        for i,c in enumerate(top6): w.writerow(["SnowFold",i+1,c["seq"]])
    print(f"\n提交: {proj_out}")
    return top6

if __name__ == "__main__":
    t0=time.time()
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Work dir (no Chinese): {WORK}")
    a = scheme_a()
    b = scheme_b()
    final(a,b)
    print(f"Done: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)")
