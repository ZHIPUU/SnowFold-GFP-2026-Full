import csv
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

from app.config import AA, SOURCE_ROOTS
from app.db import init_db

ROUND_RE = re.compile(r"round[_-]?(\d+)|r(\d+)", re.I)


def sha(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def norm_seq(s: str | None) -> str | None:
    if not s:
        return None
    s = re.sub(r"\s+", "", str(s)).upper()
    if not s or not re.match(r"^[A-Z*X-]+$", s):
        return None
    return s


def round_key_from_path(path: Path) -> str | None:
    text = str(path).replace("\\", "/")
    matches = ROUND_RE.findall(text)
    if not matches:
        return None
    last = matches[-1]
    n = last[0] or last[1]
    return f"R{int(n)}"


def file_role(path: Path) -> str:
    name = path.name.lower()
    if name.startswith("final_") or "final_6" in name or "top6" in name:
        return "final_top"
    if name.startswith("submission"):
        return "submission"
    if "all_passed" in name:
        return "all_passed"
    if "progress" in name or "phase" in name:
        return "progress"
    if name.endswith(".fa"):
        return "fasta"
    if name.endswith(".pdb"):
        return "pdb"
    if name.endswith(".jsonl"):
        return "jsonl"
    if name.endswith(".md"):
        return "document"
    if name.endswith(".py"):
        return "script"
    return "artifact"


def insert_artifact(conn, path: Path):
    st = path.stat()
    rk = round_key_from_path(path)
    cur = conn.execute(
        """INSERT OR IGNORE INTO artifacts(round_key,path,file_type,role,size_bytes,modified_time)
        VALUES(?,?,?,?,?,datetime(?, 'unixepoch'))""",
        (rk, str(path), path.suffix.lower().lstrip("."), file_role(path), st.st_size, st.st_mtime),
    )
    if cur.lastrowid:
        return cur.lastrowid
    row = conn.execute("SELECT id FROM artifacts WHERE path=?", (str(path),)).fetchone()
    return row["id"] if row else None


def upsert_round(conn, round_key: str | None, docs_path: str | None = None):
    if not round_key:
        return
    n = int(round_key[1:]) if round_key[1:].isdigit() else None
    conn.execute(
        "INSERT OR IGNORE INTO rounds(round_key, round_number, title, docs_path) VALUES(?,?,?,?)",
        (round_key, n, f"Round {n}" if n else round_key, docs_path),
    )


def upsert_sequence(conn, seq: str, round_key: str | None, parent: str | None = None, metrics: dict[str, Any] | None = None):
    sid = sha(seq)
    valid = int(all(c in AA for c in seq))
    starts = int(seq.startswith("M"))
    score = get_float(metrics, ["score", "sort_score", "final_score", "best_score"]) if metrics else None
    ptm = get_float(metrics, ["ptm", "pTM"]) if metrics else None
    plddt = normalize_conf(get_float(metrics, ["plddt", "plddt_mean", "global_plddt", "global_score"])) if metrics else None
    chromo = normalize_conf(get_float(metrics, ["chromo", "plddt_chromo", "chromo_plddt", "chromo_score"])) if metrics else None
    conn.execute(
        """INSERT OR IGNORE INTO sequences(id,sequence,length,starts_with_m,valid_aa,source_round,best_score,best_ptm,best_plddt,best_chromo,best_parent)
        VALUES(?,?,?,?,?,?,?,?,?,?,?)""",
        (sid, seq, len(seq), starts, valid, round_key, score, ptm, plddt, chromo, parent),
    )
    if score is not None:
        conn.execute(
            """UPDATE sequences SET
            best_score=CASE WHEN best_score IS NULL OR ? > best_score THEN ? ELSE best_score END,
            best_ptm=CASE WHEN best_score IS NULL OR ? >= best_score THEN ? ELSE best_ptm END,
            best_plddt=CASE WHEN best_score IS NULL OR ? >= best_score THEN ? ELSE best_plddt END,
            best_chromo=CASE WHEN best_score IS NULL OR ? >= best_score THEN ? ELSE best_chromo END,
            best_parent=CASE WHEN best_score IS NULL OR ? >= best_score THEN ? ELSE best_parent END
            WHERE id=?""",
            (score, score, score, ptm, score, plddt, score, chromo, score, parent, sid),
        )
    return sid


def get_float(d: dict[str, Any], keys: list[str]):
    for k in keys:
        if k in d and d[k] not in (None, ""):
            try:
                return float(d[k])
            except Exception:
                pass
    return None


def normalize_conf(v):
    if v is None:
        return None
    try:
        v = float(v)
    except Exception:
        return None
    # pLDDT/chromophore confidence appears as either 0-1 or 0-100 across rounds.
    # Store consistently as 0-1.
    return v / 100.0 if v > 1.5 else v


def insert_metric(conn, sid, round_key, artifact_id, item, rank=None):
    score = get_float(item, ["score", "sort_score", "final_score", "best_score"])
    ptm = get_float(item, ["ptm", "pTM"])
    plddt = normalize_conf(get_float(item, ["plddt", "plddt_mean", "global_plddt", "global_score"]))
    chromo = normalize_conf(get_float(item, ["chromo", "plddt_chromo", "chromo_plddt", "chromo_score"]))
    recycles = item.get("recycles") or item.get("best_recycles") or item.get("num_recycles")
    try:
        recycles = int(recycles) if recycles is not None else None
    except Exception:
        recycles = None
    passed = item.get("passes") if "passes" in item else item.get("all_pass")
    passed = None if passed is None else int(bool(passed))
    parent = str(item.get("parent", ""))[:200] if isinstance(item, dict) else None
    name = str(item.get("name", item.get("Seq_ID", "")))[:200]
    conn.execute(
        """INSERT INTO metrics(sequence_id,round_key,artifact_id,metric_context,score,ptm,plddt,chromo,recycles,rank,passed,parent,name)
        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        (sid, round_key, artifact_id, file_role(Path(str(artifact_id))), score, ptm, plddt, chromo, recycles, rank, passed, parent, name),
    )
    if parent:
        conn.execute(
            "INSERT INTO lineage_edges(parent_label,child_sequence_id,round_key,edge_type,weight) VALUES(?,?,?,?,?)",
            (parent, sid, round_key, "reported_parent", score),
        )


def parse_json_items(obj: Any):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        # Handle grouped results like {control: [...], r15_top3: [...]}
        out = []
        for k, v in obj.items():
            if isinstance(v, list):
                for it in v:
                    if isinstance(it, dict):
                        it = dict(it); it.setdefault("group", k); out.append(it)
            elif isinstance(v, dict) and any(x in v for x in ["seq", "Sequence", "sequence"]):
                vv = dict(v); vv.setdefault("group", k); out.append(vv)
        return out
    return []


def import_json(conn, path: Path, artifact_id: int):
    rk = round_key_from_path(path)
    try:
        obj = json.load(open(path, encoding="utf-8"))
    except Exception:
        return 0
    count = 0
    for rank, item in enumerate(parse_json_items(obj), 1):
        if not isinstance(item, dict):
            continue
        seq = norm_seq(item.get("seq") or item.get("Sequence") or item.get("sequence"))
        if not seq:
            continue
        sid = upsert_sequence(conn, seq, rk, item.get("parent"), item)
        insert_metric(conn, sid, rk, artifact_id, item, rank)
        count += 1
    return count


def import_csv(conn, path: Path, artifact_id: int):
    rk = round_key_from_path(path)
    count = 0
    try:
        with open(path, newline="", encoding="utf-8-sig") as f:
            for rank, row in enumerate(csv.DictReader(f), 1):
                seq = norm_seq(row.get("Sequence") or row.get("seq") or row.get("sequence"))
                if not seq:
                    continue
                sid = upsert_sequence(conn, seq, rk, row.get("Parent") or row.get("parent"), row)
                insert_metric(conn, sid, rk, artifact_id, row, rank)
                if "submission" in path.name.lower():
                    conn.execute(
                        "INSERT INTO submissions(round_key,sequence_id,seq_id,team_name,artifact_id) VALUES(?,?,?,?,?)",
                        (rk, sid, row.get("Seq_ID"), row.get("Team_Name"), artifact_id),
                    )
                count += 1
    except Exception:
        pass
    return count


def import_fasta(conn, path: Path, artifact_id: int):
    rk = round_key_from_path(path)
    count = 0
    name = ""
    seq = ""
    def flush():
        nonlocal count, name, seq
        s = norm_seq(seq)
        if s:
            parent = path.parent.parent.name if path.parent else None
            item = {"name": name, "parent": parent}
            sid = upsert_sequence(conn, s, rk, parent, item)
            insert_metric(conn, sid, rk, artifact_id, item, None)
            count += 1
    try:
        for line in open(path, encoding="utf-8", errors="ignore"):
            line = line.strip()
            if line.startswith(">"):
                flush(); name = line[1:]; seq = ""
            elif line:
                seq += line
        flush()
    except Exception:
        pass
    return count


def import_md(conn, path: Path, artifact_id: int):
    rk = round_key_from_path(path)
    try:
        body = open(path, encoding="utf-8", errors="ignore").read()
    except Exception:
        return 0
    title = next((ln.lstrip("# ").strip() for ln in body.splitlines() if ln.startswith("#")), path.stem)
    conn.execute("INSERT OR IGNORE INTO documents(round_key,path,title,body) VALUES(?,?,?,?)", (rk, str(path), title, body))
    return 1


def scan_files():
    exts = {".json", ".csv", ".fa", ".fasta", ".pdb", ".jsonl", ".md", ".py", ".log"}
    for root in SOURCE_ROOTS:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.is_file() and p.suffix.lower() in exts:
                # Skip huge third-party tool trees lightly
                sp = str(p).replace("\\", "/")
                if "/node_modules/" in sp or "/.git/" in sp:
                    continue
                yield p


def main():
    conn = init_db()
    imported = 0
    for p in scan_files():
        rk = round_key_from_path(p)
        upsert_round(conn, rk, None)
        aid = insert_artifact(conn, p)
        suffix = p.suffix.lower()
        n = 0
        if suffix == ".json": n = import_json(conn, p, aid)
        elif suffix == ".csv": n = import_csv(conn, p, aid)
        elif suffix in (".fa", ".fasta"): n = import_fasta(conn, p, aid)
        elif suffix == ".md": n = import_md(conn, p, aid)
        if n:
            conn.execute("UPDATE artifacts SET parsed=1 WHERE id=?", (aid,))
            imported += n
    # Update best sequence per round
    for row in conn.execute("SELECT round_key, MAX(score) AS m FROM metrics WHERE round_key IS NOT NULL GROUP BY round_key"):
        best = conn.execute("SELECT sequence_id, score FROM metrics WHERE round_key=? AND score=? LIMIT 1", (row["round_key"], row["m"])).fetchone()
        if best:
            conn.execute("UPDATE rounds SET best_score=?, best_sequence_id=? WHERE round_key=?", (best["score"], best["sequence_id"], row["round_key"]))
    conn.commit()
    print({
        "sequences": conn.execute("SELECT COUNT(*) c FROM sequences").fetchone()["c"],
        "metrics": conn.execute("SELECT COUNT(*) c FROM metrics").fetchone()["c"],
        "artifacts": conn.execute("SELECT COUNT(*) c FROM artifacts").fetchone()["c"],
        "rounds": conn.execute("SELECT COUNT(*) c FROM rounds").fetchone()["c"],
        "imported_records": imported,
    })

if __name__ == "__main__":
    main()
