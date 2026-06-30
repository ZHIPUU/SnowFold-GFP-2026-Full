#!/usr/bin/env python3
"""R22 Pipeline - Overnight long-running task (10-13h on A800)

Strategy:
- Parents: R19 Top 6 (proven effective)
- MPNN: 5 temperatures x 200 seqs = 1000 per parent
- Total: 6 parents x 1000 = 6000 candidates
- ESMFold r=8 evaluation
- Top 50 + multi-recycles consensus (r=8/16/24)
- ESM2-650M likelihood scoring (zero-shot, no training)
- Final Top 6 with combined ranking
"""
import os, sys, json, time, glob, subprocess, copy, warnings
import numpy as np
import torch
import torch.nn.functional as F
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding, AutoModelForMaskedLM

WORK = "/root/autodl-tmp/r22"
MPNN = "/root/autodl-tmp/ProteinMPNN"
os.makedirs(WORK, exist_ok=True)
for d in ["pdbs", "pdbs_precise", "mpnn_out", "results"]:
    os.makedirs(os.path.join(WORK, d), exist_ok=True)

NUM_SEQ_PER_TEMP = 200
TEMPS = [0.1, 0.2, 0.3, 0.5, 1.0]  # 5 温度
FIXED = [1, 65, 66, 67, 96, 222]  # M + 5 chromophore
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
RECYCLES_CONSENSUS = [16, 24, 32]  # 多档投票
TOP_K_SCREEN = 50
TOP_K_PRECISE = 20
BATCH = 25

# R19 Top 6 parents (proven)
R19_TOP6 = json.load(open("/root/autodl-tmp/r19/final_6_r19.json"))
print(f"R22 parents: {len(R19_TOP6)} from R19 Top 6", flush=True)

# === ESMFold (primary model) ===
print("Loading ESMFold...", flush=True)
t0 = time.time()
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
print(f"ESMFold loaded ({time.time()-t0:.0f}s)", flush=True)

# === ESM2-650M for likelihood scoring ===
print("Loading ESM2-650M...", flush=True)
t1 = time.time()
esm2_tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True)
esm2_model = AutoModelForMaskedLM.from_pretrained("facebook/esm2_t33_650M_UR50D", local_files_only=True).cuda().eval()
print(f"ESM2-650M loaded ({time.time()-t1:.0f}s)", flush=True)

aa3 = {a:b for a,b in zip("ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP