#!/usr/bin/env python3
"""R17 Pipeline - A800 full exploration w/ HuggingFace ESMFold"""
import os, sys, json, time, glob, subprocess, copy
import numpy as np
import torch
import warnings; warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

WORK = "/root/autodl-tmp/r17"
RECYCLES = 12
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs","mpnn_out","results"]:
    os.makedirs(os.path.join(WORK,d), exist_ok=True)

# Model Loading
print("Loading ESMFold (HuggingFace)...", flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print("Loaded.", flush=True)

def predict(seq, recycles=RECYCLES):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, num_recycles=recycles)
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = float(out.ptm.cpu().item())
    gp = float(plddt.mean())
    cp = float(plddt[57:72].mean())
    score = 0.40*ptm + 0.30*gp + 0.30*cp  # pLDDT is 0-1 scale on HF
    return {"ptm":round(ptm,4),"plddt":round(gp,2),"chromo":round(cp,2),
            "score":round(score,4),"passes":ptm>0.6 and gp>0.6 and cp>0.55}

# R15 Top 6 sequences
R15 = [
    ("R14_A_T02_037","MTIPGETLLAGVVPVRVNLDGDVNGKKFKVVGEGEGDATKGELRLTFEVTEGELPLDPLLLSYILTYPLSIFRKMPKDNPLRPFYLACLPEGYVIERTLDFKGEGTLNVTSETYFDGDTLVSNITLKGTGFKEGGKLLTKKVKEIRLVGDLTITPDEEKKGVKLTYTLELTFEDGSTSTADVKELIYPKGKGPETLPEPQTLTVDLTLTAVPEVEGDRFRFEHRDPPGVEPPPLSELSK"),
    ("R14_D_T02_033","MTIPGEELLSGVVPVKVNLDGDVNGQKFKVKGEGTGDATTGKLSLEFEVTEGTLPLDPLLLSYILTYPLSIFRKMPEDDPLRPFLLACLPEGYVREITLDFEGEGTLKVKSETHFEGDTLVSNITLKGTDFKEGGKLLTKKIKSIRLVGKRTITPDEERHGVNLTYTLELTWEDGSTSTAKVKELIYPIGEGPEELPEPQEITIDEVFKAKPEVNGNKFKYEYQDEPGVTPPPLSELSK"),
    ("R14_A_T01_013","MTIPGETLLSGKVPVKVNLDGDVNGEKFKVEGEGTGDATKGELKLEFKVTEGELPLDWLLLSYILTYPLSIFKKMPADHPLKPFYLCCLPEGYIIERTLDFEGEGTLKVTSKTYFDGDTLVSEITLKGTDFKEGGKLLTKEIKEIVEEGELTITPDEEKHGVKRTYTLKLTFKDGSTLTAKVDELIYPNGKGPEKLPEAQKWTINVVYKAKPKVEGDKAKVEWKDPEGVTPPPLSELSK"),
    ("R14_A_T01_020","MTIPGDTLLAGVVPVKVKLDGDVNGEKFKVEGEGEGDATTGRLKLKFEVTEGELPLDPLLLTYILTYPLRIFRKMPADDPLKPFLKACLPEGYVVERTLDFEGEGTLKVKSKTYFDGDTLVSEIELKGTDFKEGGLLLTKKVKEIRVTGELTISPDEEKHGVNLTYTLVLTFEDGSTATAKVKELIYPIGKGPETLPAPQTLNVDFVFEAKPKVEGNKFEFEWEDPPGVEPPPLSELSK"),
    ("R14_D_T02_039","MTIPGETLLNGVVPVKVNLDGDVNGEKFKVTGEGEGDATKGELKLDFEVTEGKLPLDPLLLSFILTYPLRIFRKMPADDPLKPFYLACLPEGYVVETTLDFENEGTLKVTSKTYFEGDTLVSDITLKGTDFKKGGKLLTKKVKEIRLTGDLTVTPDEAKKGVKLTYTLVLTFEDGSTSTAKVEQLIYPNGKGPETLPAAQTLTIDWVFTAKPKVEGNKFHYEYKDPEGVEPPPLSELSK"),
    ("R14_A_T01_023","MTIPGETLLSGVVPVKVNLDGDVNGNKFKVKGEGEGDATKGYLKLKFEVTEGELPLDPLLLTFILTYPLSIFRKMPADHPLRPFLLACLPEGYVVERTFDFEGEGTLKVRSDTYFDGDTLVSDITLKGTDFKEGGKLLTKKVESIKLKGELTITPDEAKHGVKLTYTLELTFTDGSTSTAKVEELIYPIGKGPATLPAAQTWTIDFTFTSKPKVEGNKFEYEYQDTPGVEPPPLSELSK"),
]

# New directions
SF = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
AV = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# avGFP + sfGFP hybrid (combining avGFP backbone with sfGFP key mutations verified on this avGFP)
av_sf_hybrid = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
sf_S30R = "MSKGEELFTGVVPILVELDGDVNGHKFRVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
sf_TGP_CTERM = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGIDYGMDELYK"

NEW = [
    ("av_sf_hybrid", av_sf_hybrid),
    ("sf_S30R", sf_S30R),
    ("sf_TGP_Cterm", sf_TGP_CTERM),
    ("sfGFP_WT", SF),
    ("avGFP_WT", AV),
]

# Run all evaluations
t0 = time.time()
print(f"R17 Start: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
all_results = []

print("\n=== R15 Top 6 Re-evaluation ===", flush=True)
for name, seq in R15:
    m = predict(seq)
    m["name"] = name; m["seq"] = seq
    all_results.append(m)
    print(f"  {name}: pTM={m['ptm']:.4f} pLDDT={m['plddt']:.2f} chromo={m['chromo']:.2f} score={m['score']:.4f}", flush=True)

print("\n=== New Directions ===", flush=True)
for name, seq in NEW:
    m = predict(seq)
    m["name"] = name; m["seq"] = seq
    all_results.append(m)
    print(f"  {name}: pTM={m['ptm']:.4f} pLDDT={m['plddt']:.2f} chromo={m['chromo']:.2f} score={m['score']:.4f}", flush=True)

# Sort and save
all_results.sort(key=lambda x: x["score"], reverse=True)
with open(os.path.join(WORK,"results","all.json"),"w") as f:
    json.dump(all_results, f, indent=2)

# Final Top 6
print(f"\n=== FINAL Top 6 (survival: pTM>0.6, pLDDT>0.6, chromo>0.55) ===", flush=True)
passed = [c for c in all_results if c["passes"]]
print(f"Passed survival: {len(passed)}/{len(all_results)}", flush=True)

top6 = passed[:6] if len(passed) >= 6 else all_results[:6]
import csv
csv_path = os.path.join(WORK,"submission_r17.csv")
with open(csv_path,"w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["Team_Name","Seq_ID","Sequence"])
    for i,c in enumerate(top6):
        w.writerow(["SnowFold",i+1,c["seq"]])
        print(f"  {i+1}. {c['name'][:28]:<30s} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.2f} chromo={c['chromo']:.2f} score={c['score']:.4f}", flush=True)

json.dump(top6, open(os.path.join(WORK,"final_6_r17.json"),"w"), indent=2)

# Also save to workspace
import shutil
try:
    shutil.copy(csv_path, os.path.join(WORK,"..", "submission_r17.csv"))
except:
    pass

print(f"\nDone: {time.strftime('%Y-%m-%d %H:%M:%S')} ({(time.time()-t0)/60:.1f} min)", flush=True)
print(f"Submit file: {csv_path}", flush=True)
