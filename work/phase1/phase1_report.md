# Phase 1.2 数据深度分析报告

## 1. 数据规模
- 总突变记录: 141572 条
- 4 种 GFP 母体: avGFP / amacGFP / cgreGFP / ppluGFP
- 单点突变(已验证): 241 条

## 2. WT 基线 (CFPS 初始亮度, log10 尺度)
- avGFP: 3.719 (相对 = 1.0)
- amacGFP: 3.971 (相对 = 1.0)
- cgreGFP: 4.497 (相对 = 1.0)
- ppluGFP: 4.226 (相对 = 1.0)

## 3. 突变数 vs 亮度 (以 avGFP 为例)
- n_mut=0: count=1, mean=3.719, max=3.719
- n_mut=1: count=1084, mean=3.431, max=4.114
- n_mut=2: count=12777, mean=3.346, max=4.108
- n_mut=3: count=12336, mean=3.013, max=4.123
- n_mut=4: count=9387, mean=2.555, max=4.107
- n_mut=5: count=6825, mean=2.111, max=4.014
- n_mut=6: count=4298, mean=1.798, max=3.954
- n_mut=7: count=2526, mean=1.600, max=3.908

## 4. 各母体 Top-10 Hotspot (按 mean delta_brightness 排序, 至少 3 次观测)

### avGFP

- pos 220.0 (L): n=3.0, mean_d=+0.031, max_d=+0.074
- pos 5.0 (E): n=4.0, mean_d=+0.020, max_d=+0.066
- pos 226.0 (A): n=5.0, mean_d=-0.036, max_d=+0.012
- pos 194.0 (L): n=5.0, mean_d=-0.073, max_d=+0.034
- pos 83.0 (F): n=4.0, mean_d=-0.139, max_d=+0.045
- pos 102.0 (D): n=3.0, mean_d=-0.173, max_d=+0.028
- pos 11.0 (V): n=3.0, mean_d=-0.186, max_d=-0.070
- pos 49.0 (T): n=3.0, mean_d=-0.228, max_d=-0.032
- pos 62.0 (T): n=5.0, mean_d=-0.320, max_d=+0.048
- pos 99.0 (F): n=5.0, mean_d=-0.723, max_d=-0.128

### amacGFP

- pos 174.0 (G): n=6.0, mean_d=+0.064, max_d=+0.084
- pos 5.0 (E): n=6.0, mean_d=-0.000, max_d=+0.039
- pos 158.0 (N): n=8.0, mean_d=-0.008, max_d=+0.016
- pos 173.0 (G): n=6.0, mean_d=-0.009, max_d=+0.019
- pos 49.0 (T): n=5.0, mean_d=-0.027, max_d=+0.008
- pos 83.0 (F): n=4.0, mean_d=-0.442, max_d=-0.048
- pos 223.0 (F): n=7.0, mean_d=-0.514, max_d=-0.106
- pos 102.0 (D): n=7.0, mean_d=-0.628, max_d=-0.056
- pos 62.0 (T): n=4.0, mean_d=-0.847, max_d=+0.075
- pos 99.0 (F): n=4.0, mean_d=-0.905, max_d=-0.061

### cgreGFP

- pos 214.0 (D): n=6.0, mean_d=-0.001, max_d=+0.022
- pos 193.0 (E): n=7.0, mean_d=-0.009, max_d=+0.016
- pos 41.0 (T): n=5.0, mean_d=-0.037, max_d=-0.013
- pos 175.0 (G): n=5.0, mean_d=-0.061, max_d=-0.013
- pos 52.0 (T): n=5.0, mean_d=-0.081, max_d=-0.009
- pos 213.0 (D): n=7.0, mean_d=-0.101, max_d=-0.063
- pos 174.0 (G): n=7.0, mean_d=-0.158, max_d=-0.049
- pos 160.0 (G): n=6.0, mean_d=-0.381, max_d=-0.056
- pos 147.0 (P): n=5.0, mean_d=-0.488, max_d=-0.043
- pos 31.0 (I): n=5.0, mean_d=-0.585, max_d=+0.006

### ppluGFP

- pos 195.0 (R): n=5.0, mean_d=+0.075, max_d=+0.109
- pos 116.0 (V): n=5.0, mean_d=-0.012, max_d=+0.014
- pos 79.0 (N): n=8.0, mean_d=-0.017, max_d=+0.004
- pos 188.0 (G): n=6.0, mean_d=-0.030, max_d=-0.008
- pos 132.0 (I): n=8.0, mean_d=-0.037, max_d=+0.034
- pos 23.0 (G): n=6.0, mean_d=-0.043, max_d=-0.014
- pos 163.0 (G): n=6.0, mean_d=-0.055, max_d=-0.012
- pos 81.0 (G): n=4.0, mean_d=-0.058, max_d=-0.019
- pos 51.0 (L): n=3.0, mean_d=-0.067, max_d=+0.001
- pos 165.0 (Y): n=4.0, mean_d=-0.070, max_d=+0.000