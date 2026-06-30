#!/bin/bash
cd /home/a/geoevo_work
export PATH="$HOME/miniconda3/bin:$PATH"
source "$HOME/miniconda3/bin/activate" geoevo

echo "=== Download ESM2-3B via HF mirror ==="
export HF_ENDPOINT=https://hf-mirror.com

python << 'PYEOF'
import os
from huggingface_hub import hf_hub_download

target = os.path.expanduser("~/.cache/torch/hub/checkpoints/esm2_t36_3B_UR50D.pt")
if os.path.exists(target):
    print(f"Already exists: {os.path.getsize(target)/1e6:.0f}MB")
else:
    print("Downloading esm2_t36_3B_UR50D.pt from HF mirror...")
    pt_path = hf_hub_download(
        repo_id="facebook/esm2_t36_3B_UR50D",
        filename="esm2_t36_3B_UR50D.pt",
        local_dir=os.path.expanduser("~/.cache/torch/hub/checkpoints"),
        local_files_only=False
    )
    if pt_path != target:
        os.rename(pt_path, target)
    print(f"Downloaded: {os.path.getsize(target)/1e6:.0f}MB")
PYEOF

echo "=== Verify ==="
ls -la ~/.cache/torch/hub/checkpoints/esm2_t36_3B_UR50D.pt
