#!/bin/bash
cd /home/a/geoevo_work

echo "=== 1. TOS ==="
export PATH="$HOME/miniconda3/bin:$PATH"
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

echo "=== 2. Create env ==="
conda create -n geoevo python=3.11 -y -c conda-forge 2>&1 | tail -3

echo "=== 3. Activate ==="
source "$HOME/miniconda3/bin/activate" geoevo 2>/dev/null || \
eval "$(conda shell.bash hook)" && conda activate geoevo

echo "=== 4. Install torch via conda ==="
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y 2>&1 | tail -5

echo "=== 5. Install other deps via conda ==="
conda install -c conda-forge biopython pandas tqdm -y 2>&1 | tail -3

echo "=== 6. Install pip deps ==="
pip install fair-esm torch_geometric -q 2>&1 | tail -3

echo "=== 7. Test ==="
python -c "import torch; print('PyTorch:', torch.__version__); print('CUDA:', torch.cuda.is_available())" 2>&1
python -c "from geoevobuilder import sequence_design; print('GeoEvoBuilder: OK')" 2>&1

echo "=== 8. Copy PDB + Run GeoEvoBuilder ==="
cp /mnt/d/生信/2026Protein\ Design/work/round4/pdbs/2B3P.pdb ./1GFL.pdb 2>/dev/null
python run_GeoEvoBuilder.py -iP ./ -i 1GFL.pdb --chainID A --SM 50 --ST 0.1 --n 30 -oP ./output/ 2>&1

echo "=== DONE ==="
ls ./output/ 2>/dev/null && cat ./output/*.fa 2>/dev/null | head -100
