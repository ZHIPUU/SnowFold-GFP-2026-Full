# Round 20 Runbook

## 复现步骤

### Step 1: 上传脚本
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
# 看最新日志
gssh --json exec 9ca7acb1b94c "tail -30 /tmp/r22.log"

# 检查 GPU + 进程
gssh --json exec 9ca7acb1b94c "nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader && ps -eo pid,etime,pcpu,cmd | grep r22 | grep -v grep"
```

### Step 4: 下载结果
```bash
# Phase 1 完成时下载 R20 Top 6
mkdir -p D:/workspace/round20
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r22/results/r20_top6.json "D:/workspace/round20/"

# 全跑完下载 R22 Top 6
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r22/final_6_r22.json "D:/workspace/round22/"
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r22/submission_r22.csv "D:/workspace/round22/"
```

## 关键参数

```python
FIXED = [1, 65, 66, 67, 96, 222]  # ← 关键: position 1 (M)
TEMPS = [0.1, 0.2, 0.5, 1.0]
NUM_SEQ_PER_TEMP = 250  # 1000/父代
RECYCLES_SCREEN = 8
RECYCLES_PRECISE = 20
TOP_K_PRECISE = 20
BATCH = 25
```

## 预计总时间

- Phase 1: 4.5h (R20 finalize 2000 候选)
- Phase 2: 7h (R21 6 父代 × 1000)
- Phase 3: 1h (Top 50 + r=20 重算)
- **总计: ~12.5h**

---

*最后更新: 2026-06-29*
