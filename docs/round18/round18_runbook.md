# Round 18 Runbook

## 复现步骤

### Step 1: 准备 ProteinMPNN 权重
```bash
gssh cp "C:\proteinmpnn_r10\protein_mpnn_run.py" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/
gssh cp "C:\proteinmpnn_r10\protein_mpnn_utils.py" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/
gssh cp "C:\proteinmpnn_r10\vanilla_model_weights\v_48_020.pt" 9ca7acb1b94c:/root/autodl-tmp/ProteinMPNN/vanilla_model_weights/
```

### Step 2: 跑 R18
```bash
gssh cp "D:\workspace\r18_pipeline.py" 9ca7acb1b94c:/root/autodl-tmp/r18_pipeline.py
gssh exec 9ca7acb1b94c "cd /root/autodl-tmp && python3 r18_pipeline.py"
```

### Step 3: 下载结果
```bash
gssh cp 9ca7acb1b94c:/root/autodl-tmp/r18/submission_r18.csv "D:\workspace\round18\"
```

## 关键参数
- 父代: R17 Top 6 (来自 `final_6_r17.json`)
- 固定残基: [65, 66, 67, 96, 222]
- 温度: 0.1 / 0.15 / 0.2
- 候选/温度: 50

---

*最后更新: 2026-06-28*
