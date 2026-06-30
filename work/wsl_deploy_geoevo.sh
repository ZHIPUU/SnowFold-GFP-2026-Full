#!/bin/bash
# WSL2 部署 GeoEvoBuilder + 运行 GFP 设计
set -e

echo "=== 1. 创建 WSL 工作目录 ==="
mkdir -p ~/geoevo && cd ~/geoevo

echo "=== 2. 从 Windows 复制文件 ==="
cp -r /mnt/d/生信/2026Protein\ Design/work/round5/GeoEvoBuilder/* ./
cp /mnt/d/生信/2026Protein\ Design/work/GeoEvoBuilder/geoevobuilder/params/Se.pt ./geoevobuilder/params/Se.pt
cp /mnt/d/生信/2026Protein\ Design/work/round5/GeoEvoBuilder/geoevobuilder/Se.pt ./geoevobuilder/Se.pt 2>/dev/null || true

echo "=== 3. 安装 Miniconda (Python 3.11) ==="
if ! command -v conda &> /dev/null; then
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p ~/miniconda3
    eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    conda init
fi
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"

echo "=== 4. 创建 conda env ==="
conda create -n geoevobuilder python=3.11 -y
conda activate geoevobuilder

echo "=== 5. 安装系统依赖 ==="
conda install -c salilab dssp -y 2>/dev/null || echo "dssp 安装跳过"
pip install uv -q

echo "=== 6. 安装 Python 依赖 ==="
uv pip install --system biopython fair-esm pandas tensorflow-cpu tqdm -q
uv pip install --system torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128 -q 2>/dev/null || \
pip install torch --index-url https://download.pytorch.org/whl/cu128 -q
uv pip install --system torch_geometric -q 2>/dev/null || pip install torch_geometric -q

echo "=== 7. 准备 GFP backbone (PDB 1GFL) ==="
# 1GFL 是 GeoEvoBuilder 论文中 GFP 设计的模板
if [ ! -f "1GFL.pdb" ]; then
    wget -q https://files.rcsb.org/download/1GFL.pdb -O 1GFL.pdb || \
    cp /mnt/d/生信/2026Protein\ Design/work/round4/pdbs/2B3P.pdb 1GFL.pdb
fi

echo "=== 8. 运行 GeoEvoBuilder GFP 设计 ==="
python run_GeoEvoBuilder.py \
    -iP ./ \
    -i 1GFL.pdb \
    --chainID A \
    --SM 50 \
    --ST 0.1 \
    --n 50 \
    --Fixed "" \
    -oP ./output/

echo "=== 9. 结果复制回 Windows ==="
cp -r ./output/ /mnt/d/生信/2026Protein\ Design/work/round8/geoevo_output/ 2>/dev/null || \
mkdir -p /mnt/d/生信/2026Protein\ Design/work/round8/geoevo_output && \
cp -r ./output/* /mnt/d/生信/2026Protein\ Design/work/round8/geoevo_output/ 2>/dev/null

echo "=== 完成 ==="
ls -la ./output/
