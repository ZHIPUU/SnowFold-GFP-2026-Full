#!/bin/bash
cd /home/a/geoevo_work
export PATH="$HOME/miniconda3/bin:$PATH"
source "$HOME/miniconda3/bin/activate" geoevo 2>/dev/null || eval "$(conda shell.bash hook)" && conda activate geoevo

echo "=== Fix torch ==="
pip uninstall torch torchvision torchaudio -y 2>/dev/null || true
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu -q 2>&1 | tail -3

echo "=== Install TF ==="
pip install tensorflow-cpu -q 2>&1 | tail -3

echo "=== Test ==="
python -c "import torch; print('Torch:', torch.__version__)" 2>&1
python -c "import tensorflow; print('TF:', tensorflow.__version__)" 2>&1
python -c "from geoevobuilder import sequence_design; print('GeoEvoBuilder: OK')" 2>&1

echo "=== Run GeoEvoBuilder ==="
cp /mnt/d/生信/2026Protein\ Design/work/round4/pdbs/2B3P.pdb ./1GFL.pdb 2>/dev/null
python run_GeoEvoBuilder.py -iP ./ -i 1GFL.pdb --chainID A --SM 50 --ST 0.1 --n 30 -oP ./output/ 2>&1
echo "=== DONE ==="
ls ./output/ 2>/dev/null
