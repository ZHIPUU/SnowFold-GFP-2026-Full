import os
from pathlib import Path

BASE_DIR = Path(os.getenv("ATLAS_BASE_DIR", "/app"))
DB_PATH = Path(os.getenv("ATLAS_DB", str(BASE_DIR / "data" / "design.db")))
SOURCE_ROOTS = [Path(p) for p in os.getenv(
    "ATLAS_SOURCE_ROOTS",
    "D:/workspace;D:/生信/2026Protein Design"
).split(";") if p]
ARTIFACT_DIR = Path(os.getenv("ATLAS_ARTIFACT_DIR", str(BASE_DIR / "artifacts")))

AA = set("ACDEFGHIKLMNPQRSTVWY")
