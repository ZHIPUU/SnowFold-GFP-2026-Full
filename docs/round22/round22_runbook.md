# Round 22 Runbook

## 复现步骤

### Step 1: 上传文件到服务器

```bash
gssh cp "D:\workspace\r22_long.py" 9ca7acb1b94c:/root/autodl-tmp/r22_long.py
gssh cp "C:\proteinmpnn_r10\protein_mpnn_run.py" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/
gssh cp "C:\proteinmpnn_r10\protein_mpnn_utils.py" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/
gssh cp "C:\proteinmpnn_r10\vanilla_model_weights\v_48_020.pt" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/vanilla_model_weights/
```

### Step 2: 启动 R22 长线任务

```bash
gssh run 9ca7acb1b94c "cd /root/autodl-tmp && stdbuf -oL -eL python3 -u r22_long.py > /tmp/r22.log 2>&1"
```

### Step 3: 监控

```bash
# 看进度
gssh --json exec 9ca7acb1b94c "tail -50 /tmp/r22.log"

# 看 GPU/进程
gssh --json exec 9ca7acb1b94c "nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader && ps -eo pid,etime,pcpu,cmd | grep r22 | grep -v grep"
```

### Step 4: 下载结果

```bash
mkdir -p D:/workspace/round22
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r22/submission_r22.csv "D:/workspace/round22/"
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r22/final_6_r22.json "D:/workspace/round22/"
```

### Step 5: 合规检查

```bash
python D:/workspace/round22/analyze_r22.py
```

## 关键参数

```python
FIXED = [1, 65, 66, 67, 96, 222]  # ← 关键! position 1 (M) + 5 chromophore
TEMPS = [0.1, 0.2, 0.5, 1.0]      # 4 档
NUM_SEQ_PER_TEMP = 250            # Phase 1: 1000 候选/父代
NUM_SEQ_PER_TEMP_R21 = 150       # Phase 2: 600 候选/父代
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20            # 项目规则 5.3
TOP_K_PRECISE = 20
BATCH = 25
```

## 流程

```
Phase 1: R20 finalize        (4.5h, 2000 candidates)
       ↓
Phase 2: R21 large MPNN      (7h, 3600 candidates)
       ↓
Phase 3: Top 50 + r=20 recount (30min)
       ↓
Final Top 6 (0.9430)
```

## 预计总时间

~12 小时（Phase 1+2+3）

---

*最后更新: 2026-06-29*
