from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pathlib import Path
import sqlite3
import json

from .db import init_db, connect
from .config import DB_PATH, ARTIFACT_DIR

app = FastAPI(
    title="Protein Design Atlas API",
    description="SQLite-backed atlas for GFP/protein design rounds, sequences, metrics, lineage, artifacts, and plots.",
    version="0.1.0",
)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

@app.on_event("startup")
def startup():
    init_db(DB_PATH)


def rows(sql, params=()):
    conn = connect(DB_PATH)
    cur = conn.execute(sql, params)
    out = [dict(r) for r in cur.fetchall()]
    conn.close()
    return out

@app.get("/api/health")
def health():
    return {"ok": True, "db": str(DB_PATH), "exists": DB_PATH.exists()}

@app.post("/api/import/run")
def run_import():
    from .services.importer import main
    main()
    return stats()

@app.get("/api/stats")
def stats():
    conn = connect(DB_PATH)
    result = {}
    for t in ["rounds", "sequences", "metrics", "artifacts", "documents", "submissions", "lineage_edges"]:
        result[t] = conn.execute(f"SELECT COUNT(*) c FROM {t}").fetchone()["c"]
    result["best"] = dict(conn.execute("SELECT id, best_score, best_ptm, best_plddt, best_chromo, source_round FROM sequences WHERE best_score BETWEEN 0 AND 1 ORDER BY best_score DESC LIMIT 1").fetchone() or {})
    conn.close()
    return result

@app.get("/api/rounds")
def get_rounds():
    return rows("SELECT * FROM rounds ORDER BY COALESCE(round_number, 9999), round_key")

@app.get("/api/sequences")
def get_sequences(round_key: str | None = None, min_score: float | None = None, limit: int = Query(200, le=2000), offset: int = 0):
    where = []
    params = []
    if round_key:
        where.append("source_round=?"); params.append(round_key)
    if min_score is not None:
        where.append("best_score>=?"); params.append(min_score)
    clause = "WHERE " + " AND ".join(where) if where else ""
    params += [limit, offset]
    return rows(f"SELECT * FROM sequences {clause} ORDER BY best_score DESC NULLS LAST LIMIT ? OFFSET ?", params)

@app.get("/api/sequences/{sequence_id}")
def get_sequence(sequence_id: str):
    seq = rows("SELECT * FROM sequences WHERE id=?", (sequence_id,))
    metrics = rows("SELECT * FROM metrics WHERE sequence_id=? ORDER BY score DESC NULLS LAST", (sequence_id,))
    edges = rows("SELECT * FROM lineage_edges WHERE child_sequence_id=? OR parent_sequence_id=?", (sequence_id, sequence_id))
    return {"sequence": seq[0] if seq else None, "metrics": metrics, "edges": edges}

@app.get("/api/metrics/top")
def top_metrics(limit: int = Query(50, le=1000), round_key: str | None = None):
    if round_key:
        return rows("SELECT m.*, s.sequence FROM metrics m JOIN sequences s ON s.id=m.sequence_id WHERE m.round_key=? AND m.score BETWEEN 0 AND 1 ORDER BY m.score DESC LIMIT ?", (round_key, limit))
    return rows("SELECT m.*, s.sequence FROM metrics m JOIN sequences s ON s.id=m.sequence_id WHERE m.score BETWEEN 0 AND 1 ORDER BY m.score DESC LIMIT ?", (limit,))

@app.get("/api/artifacts")
def get_artifacts(round_key: str | None = None, role: str | None = None, limit: int = 300):
    where=[]; params=[]
    if round_key: where.append("round_key=?"); params.append(round_key)
    if role: where.append("role=?"); params.append(role)
    clause="WHERE "+" AND ".join(where) if where else ""
    params.append(limit)
    return rows(f"SELECT * FROM artifacts {clause} ORDER BY modified_time DESC LIMIT ?", params)

@app.get("/api/documents")
def get_documents(round_key: str | None = None):
    if round_key:
        return rows("SELECT id, round_key, path, title, substr(body,1,1000) AS preview FROM documents WHERE round_key=? ORDER BY path", (round_key,))
    return rows("SELECT id, round_key, path, title, substr(body,1,1000) AS preview FROM documents ORDER BY round_key, path")

@app.get("/api/documents/{doc_id}")
def get_document(doc_id: int):
    r = rows("SELECT * FROM documents WHERE id=?", (doc_id,))
    return r[0] if r else {}

@app.get("/api/graph/lineage")
def lineage_graph(min_score: float = 0.90, limit: int = 500):
    nodes = rows("SELECT id, source_round, best_score, best_ptm, best_plddt, best_chromo, length FROM sequences WHERE best_score>=? ORDER BY best_score DESC LIMIT ?", (min_score, limit))
    ids = {n["id"] for n in nodes}
    edges = rows("SELECT * FROM lineage_edges WHERE child_sequence_id IN (%s) LIMIT 2000" % ",".join("?" for _ in ids), tuple(ids)) if ids else []
    cy_nodes = [{"data": {"id": n["id"], "label": f"{n.get('source_round') or ''}\n{(n.get('best_score') or 0):.4f}", **n}} for n in nodes]
    cy_edges = [{"data": {"id": f"e{e['id']}", "source": e.get("parent_sequence_id") or e.get("parent_label") or "unknown", "target": e["child_sequence_id"], **e}} for e in edges]
    # Add parent label placeholder nodes if needed
    existing = {n["data"]["id"] for n in cy_nodes}
    for e in cy_edges:
        src = e["data"]["source"]
        if src not in existing:
            cy_nodes.append({"data": {"id": src, "label": str(src)[:32], "source_round": "parent_label", "best_score": 0}})
            existing.add(src)
    return {"nodes": cy_nodes, "edges": cy_edges}

@app.get("/api/plots/round-trend")
def round_trend():
    return rows("SELECT round_key, round_number, best_score FROM rounds WHERE best_score IS NOT NULL ORDER BY round_number")

@app.get("/api/plots/score-scatter")
def score_scatter(limit: int = 5000):
    return rows("SELECT s.id, s.source_round, s.best_score AS score, s.best_ptm AS ptm, s.best_plddt AS plddt, s.best_chromo AS chromo, s.length FROM sequences s WHERE s.best_score IS NOT NULL ORDER BY s.best_score DESC LIMIT ?", (limit,))

@app.get("/api/artifacts/file")
def file(path: str):
    p = Path(path)
    if not p.exists():
        return {"error": "missing"}
    return FileResponse(p)
