"""稳健下载 LigandMPNN 模型权重 (Python urllib + retry)"""
import urllib.request
import os
import time
from pathlib import Path

OUT_DIR = Path(r"D:\生信\2026Protein Design\work\round5\LigandMPNN\model_params")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 关键文件 (Round 5 优先用 ligandmpnn + soluble + proteinmpnn baseline)
FILES = [
    "ligandmpnn_v_32_010_25.pt",  # 主力, 中等噪声
    "ligandmpnn_v_32_020_25.pt",  # 高噪声 = 多样性
    "ligandmpnn_v_32_005_25.pt",  # 低噪声 = 保守
    "solublempnn_v_48_020.pt",    # 可溶性 backbone
    "proteinmpnn_v_48_020.pt",    # baseline 对比
]

BASE = "https://files.ipd.uw.edu/pub/ligandmpnn/"

def download(url, out, retries=3, chunk=8192):
    """带重试和进度的下载"""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=120) as resp:
                total = int(resp.headers.get('Content-Length', 0))
                downloaded = 0
                t0 = time.time()
                with open(out, 'wb') as f:
                    while True:
                        buf = resp.read(chunk)
                        if not buf:
                            break
                        f.write(buf)
                        downloaded += len(buf)
                        if downloaded % (chunk * 100) == 0:
                            pct = downloaded * 100 / total if total else 0
                            speed = downloaded / (time.time() - t0) / 1024
                            print(f"  ... {downloaded/1024/1024:.1f}MB ({pct:.0f}%) {speed:.0f} KB/s", end="\r")
            print(f"  ✓ {downloaded/1024/1024:.2f}MB")
            return True
        except Exception as e:
            print(f"  attempt {attempt+1} failed: {str(e)[:80]}")
            if attempt < retries - 1:
                time.sleep(3)
    return False


print(f"下载到 {OUT_DIR}\n")
for fname in FILES:
    out = OUT_DIR / fname
    if out.exists() and out.stat().st_size > 5_000_000:  # > 5 MB
        print(f"[skip] {fname} ({out.stat().st_size/1024/1024:.2f}MB)")
        continue
    print(f"[download] {fname}")
    url = BASE + fname
    success = download(url, out)
    if not success:
        print(f"  ✗ failed after retries")
    print()

print("\n=== 最终状态 ===")
for f in sorted(OUT_DIR.glob("*.pt")):
    print(f"  {f.name}: {f.stat().st_size/1024/1024:.2f}MB")
