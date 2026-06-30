"""检查 sfGFP 经典突变在 avGFP 加性模型里的效应"""
import json
import pickle
from pathlib import Path

with open(r'D:\生信\2026Protein Design\work\phase1\model_report_v2.json') as f:
    rep = json.load(f)

# 找 avGFP 的所有 effects
import sys
sys.path.insert(0, str(Path(r'D:\生信\2026Protein Design\work\phase1')))
with open(r'D:\生信\2026Protein Design\work\phase1\additive_models_v2.pkl', 'rb') as f:
    models = pickle.load(f)

av_effs = models['avGFP']['effects']

# 经典 sfGFP 突变 (在 avGFP 上下文)
# avGFP: S65 S72 M153 V163 T203 S202 A206 I171 V163 ...
classical = [
    ('S', 65, 'T'),  # S65T
    ('S', 72, 'A'),  # S72A
    ('M', 153, 'T'), # M153T
    ('V', 163, 'A'), # V163A
    ('T', 203, 'I'), # T203I
    ('T', 203, 'Y'), # T203Y
    ('S', 202, 'D'), # S202D
    ('A', 206, 'V'), # A206V
    ('A', 206, 'K'), # A206K
    ('I', 171, 'V'), # I171V
    ('N', 149, 'Y'), # N149Y
    ('F', 99, 'S'),  # F99S
    ('N', 105, 'K'), # N105K
    ('K', 166, 'T'), # K166T
    ('I', 167, 'T'), # I167T
    ('K', 79, 'R'),  # K79R
    ('S', 30, 'R'),  # S30R
    ('E', 172, 'G'), # E172G
]
print('=== avGFP 加性模型对 sfGFP 经典突变的预测 ===')
print('(正值 = 提升亮度, 负值 = 降低亮度)')
for wt_aa, pos, new_aa in classical:
    # effects key 是 (pos, wt_aa, new_aa) 或 (pos, new_aa)
    key_variants = [(pos, wt_aa, new_aa), (pos, new_aa)]
    eff = None
    for k in key_variants:
        if k in av_effs:
            eff = av_effs[k]
            break
    if eff is None:
        print(f'  {wt_aa}{pos}{new_aa}: NOT IN MODEL')
    else:
        sign = '+' if eff > 0 else ''
        print(f'  {wt_aa}{pos}{new_aa}: {sign}{eff:.3f}')

# 输出所有 positive 效应的前 30
print('\n=== avGFP top 30 positive mutations ===')
sorted_effs = sorted(av_effs.items(), key=lambda x: -x[1])
for key, e in sorted_effs[:30]:
    if len(key) == 3:
        pos, wt_aa, new_aa = key
    else:
        pos, new_aa = key
    print(f'  {wt_aa}{pos}{new_aa}: +{e:.3f}')

print('\n=== avGFP top 30 negative mutations ===')
for key, e in sorted_effs[-30:]:
    if len(key) == 3:
        pos, wt_aa, new_aa = key
    else:
        pos, new_aa = key
    print(f'  {wt_aa}{pos}{new_aa}: {e:.3f}')