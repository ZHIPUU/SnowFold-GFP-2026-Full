import sqlite3
from pathlib import Path
from .config import DB_PATH

SCHEMA = """
PRAGMA journal_mode=WAL;
PRAGMA foreign_keys=ON;

CREATE TABLE IF NOT EXISTS rounds (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  round_key TEXT UNIQUE,
  round_number INTEGER,
  title TEXT,
  summary TEXT,
  status TEXT,
  best_score REAL,
  best_sequence_id TEXT,
  docs_path TEXT,
  created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS sequences (
  id TEXT PRIMARY KEY,
  sequence TEXT NOT NULL UNIQUE,
  length INTEGER,
  starts_with_m INTEGER,
  valid_aa INTEGER,
  source_round TEXT,
  best_score REAL,
  best_ptm REAL,
  best_plddt REAL,
  best_chromo REAL,
  best_parent TEXT,
  created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS experiments (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  round_key TEXT,
  name TEXT,
  experiment_type TEXT,
  generator TEXT,
  temperatures TEXT,
  fixed_positions TEXT,
  recycles TEXT,
  candidate_count INTEGER,
  passed_count INTEGER,
  notes TEXT,
  FOREIGN KEY(round_key) REFERENCES rounds(round_key)
);

CREATE TABLE IF NOT EXISTS metrics (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  sequence_id TEXT,
  round_key TEXT,
  artifact_id INTEGER,
  metric_context TEXT,
  score REAL,
  ptm REAL,
  plddt REAL,
  chromo REAL,
  recycles INTEGER,
  rank INTEGER,
  passed INTEGER,
  parent TEXT,
  name TEXT,
  FOREIGN KEY(sequence_id) REFERENCES sequences(id),
  FOREIGN KEY(round_key) REFERENCES rounds(round_key)
);

CREATE TABLE IF NOT EXISTS lineage_edges (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  parent_label TEXT,
  parent_sequence_id TEXT,
  child_sequence_id TEXT,
  round_key TEXT,
  edge_type TEXT,
  weight REAL,
  hamming_distance INTEGER,
  FOREIGN KEY(child_sequence_id) REFERENCES sequences(id)
);

CREATE TABLE IF NOT EXISTS artifacts (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  round_key TEXT,
  path TEXT UNIQUE,
  file_type TEXT,
  role TEXT,
  size_bytes INTEGER,
  modified_time TEXT,
  parsed INTEGER DEFAULT 0
);

CREATE TABLE IF NOT EXISTS submissions (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  round_key TEXT,
  sequence_id TEXT,
  seq_id TEXT,
  team_name TEXT,
  artifact_id INTEGER,
  FOREIGN KEY(sequence_id) REFERENCES sequences(id)
);

CREATE TABLE IF NOT EXISTS documents (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  round_key TEXT,
  path TEXT UNIQUE,
  title TEXT,
  body TEXT
);

CREATE INDEX IF NOT EXISTS idx_metrics_score ON metrics(score DESC);
CREATE INDEX IF NOT EXISTS idx_metrics_round ON metrics(round_key);
CREATE INDEX IF NOT EXISTS idx_sequences_best ON sequences(best_score DESC);
CREATE INDEX IF NOT EXISTS idx_edges_child ON lineage_edges(child_sequence_id);
"""

def connect(path: Path | None = None):
    db = Path(path or DB_PATH)
    db.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys=ON")
    return conn

def init_db(path: Path | None = None):
    conn = connect(path)
    conn.executescript(SCHEMA)
    conn.commit()
    return conn
