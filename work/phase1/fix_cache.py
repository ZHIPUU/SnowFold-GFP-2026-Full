"""快速修复: 把 phase1_cache 重新保存为 dict-only 格式"""
import pickle
import sys
from pathlib import Path

# dummy parse_mut 注入到 __main__, 这样 pickle.load 不会报错
def parse_mut(s):
    return []

sys.modules['__main__'].parse_mut = parse_mut

PHASE1 = Path(r"D:\生信\2026Protein Design\work\phase1")

with open(PHASE1 / "phase1_cache.pkl", "rb") as f:
    cache = pickle.load(f)

# 去掉 parse_mut (不要函数引用)
cache_clean = {k: v for k, v in cache.items() if k != "parse_mut"}
print("Saving without parse_mut:", list(cache_clean.keys()))

with open(PHASE1 / "phase1_cache_v2.pkl", "wb") as f:
    pickle.dump(cache_clean, f)
print("Done")
