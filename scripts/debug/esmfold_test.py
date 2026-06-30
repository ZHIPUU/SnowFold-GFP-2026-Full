import torch
from transformers import AutoTokenizer, EsmForProteinFolding
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", local_files_only=True)
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True, local_files_only=True).cuda()
model.trunk.set_chunk_size(128); model.eval()
seq = "MTIPGETLLSGVVPVKVNLDGDVNGKKFKVKGEGEGDAEKGVLKVTFTVTEGELPLDPLLLSYILTYPLRIFKKMPKDDPLKPFLLSCLPSGYVIERELDFEGEGTLKVTSKVYFEGDTLVNEVTLKGSGFKEGSKLLTKKVASIRLVGDLTITPDEEKKGVKVTYTLELTFEDGSTSTADVKELIYPKGKGPETLPEPQTLTVDLTLTAVPEVEGDRFRFEHRDPPGVEPPPLSELSY"
print(f"seq len: {len(seq)}")
inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
inputs = {k: v.cuda() for k, v in inputs.items()}
with torch.no_grad():
    out = model(**inputs, num_recycles=8)
print(f"pTM: {out.ptm.item():.4f}")
plddt = out.plddt[0,:,1].cpu().numpy()
print(f"pLDDT mean: {plddt.mean():.3f}")
print(f"VRAM max: {torch.cuda.max_memory_allocated()/1024**3:.1f} GB")
