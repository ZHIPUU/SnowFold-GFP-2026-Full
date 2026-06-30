#!/usr/bin/env python3
"""R16 广度探索：在 A800 上用 ESMFold 预测 R15 候选 + avGFP/sfGFP 对照 + 新方向序列"""
import warnings, torch, time, json, sys
warnings.filterwarnings("ignore")
from transformers import AutoTokenizer, EsmForProteinFolding

def load_model():
    print("Loading ESMFold...", flush=True)
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).cuda()
    model.trunk.set_chunk_size(128)
    model.eval()
    print("ESMFold loaded.", flush=True)
    return tok, model

def predict(tok, model, seq, num_recycles=8):
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.cuda() for k, v in inputs.items()}
    t0 = time.time()
    with torch.no_grad():
        out = model(**inputs, num_recycles=num_recycles)
    t1 = time.time()
    plddt = out.plddt[0, :, 1].cpu().numpy()
    ptm = out.ptm.item()
    global_plddt = float(plddt.mean())
    chromo_plddt = float(plddt[57:72].mean())
    score = ptm * 0.40 + global_plddt / 100 * 0.30 + chromo_plddt / 100 * 0.30
    return {
        "seq_len": len(seq),
        "ptm": round(ptm, 4),
        "global_plddt": round(global_plddt, 2),
        "chromo_plddt": round(chromo_plddt, 2),
        "sort_score": round(score, 4),
        "time_s": round(t1 - t0, 1),
        "num_recycles": num_recycles,
    }

# R15 Top 6 候选
R15_TOP6 = [
    ("R14_A_T02_037", "MTIPGETLLAGVVPVRVNLDGDVNGKKFKVVGEGEGDATKGELRLTFEVTEGELPLDPLLLSYILTYPLSIFRKMPKDNPLRPFYLACLPEGYVIERTLDFKGEGTLNVTSETYFDGDTLVSNITLKGTGFKEGGKLLTKKVKEIRLVGDLTITPDEEKKGVKLTYTLELTFEDGSTSTADVKELIYPKGKGPETLPEPQTLTVDLTLTAVPEVEGDRFRFEHRDPPGVEPPPLSELSK"),
    ("R14_D_T02_033", "MTIPGEELLSGVVPVKVNLDGDVNGQKFKVKGEGTGDATTGKLSLEFEVTEGTLPLDPLLLSYILTYPLSIFRKMPEDDPLRPFLLACLPEGYVREITLDFEGEGTLKVKSETHFEGDTLVSNITLKGTDFKEGGKLLTKKIKSIRLVGKRTITPDEERHGVNLTYTLELTWEDGSTSTAKVKELIYPIGEGPEELPEPQEITIDEVFKAKPEVNGNKFKYEYQDEPGVTPPPLSELSK"),
    ("R14_A_T01_013", "MTIPGETLLSGKVPVKVNLDGDVNGEKFKVEGEGTGDATKGELKLEFKVTEGELPLDWLLLSYILTYPLSIFKKMPADHPLKPFYLCCLPEGYIIERTLDFEGEGTLKVTSKTYFDGDTLVSEITLKGTDFKEGGKLLTKEIKEIVEEGELTITPDEEKHGVKRTYTLKLTFKDGSTLTAKVDELIYPNGKGPEKLPEAQKWTINVVYKAKPKVEGDKAKVEWKDPEGVTPPPLSELSK"),
]

# 对照组
CONTROLS = [
    ("sfGFP_WT", "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"),
    ("avGFP_WT", "MSKGEELFTGVPVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"),
]

# 广度探索：新方向序列
# 方向1: avGFP + sfGFP 关键突变混合 (R15 handoff 建议)
# 方向2: C端延伸稳定化 (R15 指出 C端 pLDDT 偏低)
# 方向3: 二硫键引入 (R15 handoff P1 建议)
EXPLORATION = [
    # 方向1: avGFP 骨架 + sfGFP 的 6 个关键稳定突变 (S30R, Y39N, S72A, F99S, M153T, V163A)
    ("avGFP_sf6mut", "MSKGEELFTGVPVPILVELDGDVNGHKFSVRGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"),
    # 方向2: sfGFP + C端 GITH->GIDY 突变 (延长稳定螺旋)
    ("sfGFP_cstable", "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGIDYGMDELYK"),
    # 方向3: R14_A_T01_013 + 引入二硫键 (C48-C70)
    ("R14A013_disulfide", "MTIPGETLLSGKVPVKVNLDGDVNGEKFKVEGEGTGDATKGELKLEFKVTEGELPLDWLLLSYILCYPLSIFKKMPADHPLKPFYLCCLPEGYIIERTLDFEGEGTLKVTSKTYFDGDTLVSEITLKGTDFKEGGKLLTKEIKEIVEEGELTITPDEEKHGVKRTYTLKLTFKDGSTLTAKVDELIYPNGKGPEKLPEAQKWTINVVYKAKPKVEGDKAKVEWKDPEGVTPPPLSELSK"),
]

if __name__ == "__main__":
    tok, model = load_model()
    results = {}
    
    all_seqs = [("control", CONTROLS), ("r15_top3", R15_TOP6[:3]), ("exploration", EXPLORATION)]
    
    for group_name, seqs in all_seqs:
        print(f"\n=== {group_name} ===", flush=True)
        results[group_name] = []
        for name, seq in seqs:
            print(f"  Predicting {name} ({len(seq)}aa)...", flush=True)
            r8 = predict(tok, model, seq, num_recycles=8)
            r12 = predict(tok, model, seq, num_recycles=12)
            entry = {"name": name, "seq": seq, "r8": r8, "r12": r12}
            results[group_name].append(entry)
            print(f"    r8: score={r8['sort_score']} pTM={r8['ptm']} pLDDT={r8['global_plddt']} chromo={r8['chromo_plddt']} ({r8['time_s']}s)", flush=True)
            print(f"    r12: score={r12['sort_score']} pTM={r12['ptm']} pLDDT={r12['global_plddt']} chromo={r12['chromo_plddt']} ({r12['time_s']}s)", flush=True)
    
    # 保存结果
    with open("/root/autodl-tmp/r16_explore_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\n=== Results saved to /root/autodl-tmp/r16_explore_results.json ===", flush=True)
    
    # 汇总
    print("\n=== SUMMARY ===", flush=True)
    for group_name, entries in results.items():
        for e in entries:
            r8 = e["r8"]
            r12 = e["r12"]
            best = max(r8["sort_score"], r12["sort_score"])
            best_r = 8 if r8["sort_score"] >= r12["sort_score"] else 12
            print(f"  {e['name']:25s} best_score={best:.4f} (r={best_r}) pTM={max(r8['ptm'],r12['ptm']):.4f} chromo={max(r8['chromo_plddt'],r12['chromo_plddt']):.1f}", flush=True)
    
    print("\nDONE", flush=True)
