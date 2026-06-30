"""Launch step2_train.py in background, properly detached."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
SCRIPT = ROOT / "step2_train.py"
OUT_LOG = ROOT / "step2_stdout.log"
ERR_LOG = ROOT / "step2_stderr.log"
PY = Path(r"C:\Python314\python.exe")

# 0x00000008 = DETACHED_PROCESS, 0x00000200 = CREATE_NEW_PROCESS_GROUP
DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200

with open(OUT_LOG, "wb") as out, open(ERR_LOG, "wb") as err:
    p = subprocess.Popen(
        [str(PY), "-u", str(SCRIPT)],
        cwd=str(ROOT),
        stdout=out,
        stderr=err,
        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
        close_fds=True,
    )

print(f"Launched step2_train.py as PID {p.pid}")
print(f"  stdout -> {OUT_LOG}")
print(f"  stderr -> {ERR_LOG}")
print(f"  cwd    -> {ROOT}")