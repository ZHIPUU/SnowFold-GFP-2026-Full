"""Launch stepC epistasis training in background."""
import subprocess
from pathlib import Path

ROOT = Path(r"D:\生信\2026Protein Design\work\round2")
SCRIPT = ROOT / "stepC_epistasis.py"
OUT_LOG = ROOT / "stepC_stdout.log"
ERR_LOG = ROOT / "stepC_stderr.log"
PY = Path(r"C:\Python314\python.exe")

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
print(f"Launched stepC as PID {p.pid}")