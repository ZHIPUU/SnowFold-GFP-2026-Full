"""Test fair-esm ESMFold vs HuggingFace ESMFold pLDDT values"""
import warnings, time, json
warnings.filterwarnings("ignore")
import esm
import torch
import torch.nn.functional as F

print("Loading fair-esm ESMFold v1...")
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
model.set_chunk_size(128)
print("Loaded.")

sequences = {
    "sfGFP_WT": "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "avGFP_WT": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
}

# Load R15 Top 1
with open("/root/autodl-tmp/r15_top1_seq.txt", "r") as f:
    r15_seq = f.read().strip()
sequences["R15_Top1"] = r15_seq

results = {}
for name, seq in sequences.items():
    print(f"\n--- {name} (len={len(seq)}) ---")
    t0 = time.time()
    with torch.no_grad():
        out = model.infer(seq, num_recycles=8)
    elapsed = time.time() - t0
    
    ptm = float(out["ptm"].cpu().item())
    
    # pLDDT from binned logits (same method as before)
    plddt_logits = out["plddt"]  # (L, 37)
    probs = F.softmax(plddt_logits, dim=-1)
    centers = torch.linspace(0.5/37, 1-0.5/37, 37, device=plddt_logits.device)
    plddt_01 = (probs * centers.unsqueeze(0)).sum(-1)
    plddt_100 = (plddt_01 * 100).cpu().numpy()
    
    global_plddt = float(plddt_100.mean())
    chromo_plddt = float(plddt_100[57:72].mean())
    score = 0.40 * ptm + 0.30 * global_plddt / 100 + 0.30 * chromo_plddt / 100
    
    print(f"  pTM={ptm:.4f}, Global pLDDT={global_plddt:.1f}, Chromo pLDDT={chromo_plddt:.1f}")
    print(f"  Score={score:.4f}, Time={elapsed:.1f}s")
    
    results[name] = {
        "ptm": round(ptm, 4),
        "global_plddt": round(global_plddt, 1),
        "chromo_plddt": round(chromo_plddt, 1),
        "score": round(score, 4),
        "time_s": round(elapsed, 1)
    }

print("\n\n=== SUMMARY ===")
print(json.dumps(results, indent=2))
