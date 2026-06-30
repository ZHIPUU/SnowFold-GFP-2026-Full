import sqlite3, csv, json
from pathlib import Path

ROOT = Path(r'D:\workspace\protein-design-atlas')
DB = ROOT / 'data' / 'design.db'
OUT = ROOT / 'artifacts' / 'report_v1_data'
OUT.mkdir(parents=True, exist_ok=True)
con = sqlite3.connect(DB)
con.row_factory = sqlite3.Row

def dump(sql, path):
    rows = con.execute(sql).fetchall()
    if not rows:
        return 0
    with open(OUT / path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows([dict(r) for r in rows])
    return len(rows)

counts = {}
counts['rounds'] = dump("""
SELECT round_key, round_number, best_score FROM rounds
WHERE best_score BETWEEN 0 AND 1
ORDER BY round_number
""", 'rounds.csv')
counts['seq_scores'] = dump("""
SELECT id, source_round, best_score AS score, best_ptm AS ptm, best_plddt AS plddt,
       best_chromo AS chromo, best_parent AS parent, length, sequence
FROM sequences
WHERE best_score BETWEEN 0 AND 1 AND best_ptm IS NOT NULL AND best_chromo IS NOT NULL
ORDER BY best_score DESC
LIMIT 20000
""", 'seq_scores.csv')
counts['top_candidates'] = dump("""
SELECT id, source_round, best_score AS score, best_ptm AS ptm, best_plddt AS plddt,
       best_chromo AS chromo, best_parent AS parent, length, sequence
FROM sequences
WHERE best_score BETWEEN 0 AND 1
ORDER BY best_score DESC
LIMIT 80
""", 'top_candidates.csv')
counts['metrics_top'] = dump("""
SELECT round_key, parent, score, ptm, plddt, chromo, name, sequence
FROM metrics m JOIN sequences s ON s.id=m.sequence_id
WHERE score BETWEEN 0 AND 1
ORDER BY score DESC
LIMIT 1000
""", 'metrics_top1000.csv')

# Hamming matrix for top 40 sequences with equal length
seqs = [dict(r) for r in con.execute("""
SELECT id, source_round, best_score AS score, sequence
FROM sequences WHERE best_score BETWEEN 0 AND 1
ORDER BY best_score DESC LIMIT 40
""").fetchall()]
mat_rows=[]
for i,a in enumerate(seqs):
    for j,b in enumerate(seqs):
        L=min(len(a['sequence']), len(b['sequence']))
        d=sum(x!=y for x,y in zip(a['sequence'][:L], b['sequence'][:L]))
        mat_rows.append({'i':i+1,'j':j+1,'seq_i':f"{a['source_round']}_{i+1}",'seq_j':f"{b['source_round']}_{j+1}",'score_i':a['score'],'score_j':b['score'],'hamming':d,'identity':1-d/max(L,1)})
with open(OUT/'hamming_top40.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f, fieldnames=mat_rows[0].keys()); w.writeheader(); w.writerows(mat_rows)
counts['hamming_top40']=len(mat_rows)

# Parent contribution table
counts['parent_contrib'] = dump("""
SELECT COALESCE(best_parent, 'unknown') AS parent, source_round,
       COUNT(*) AS n,
       MAX(best_score) AS max_score,
       AVG(best_score) AS mean_score,
       AVG(best_ptm) AS mean_ptm,
       AVG(best_chromo) AS mean_chromo
FROM sequences
WHERE best_score >= 0.93 AND best_score <= 1
GROUP BY source_round, parent
ORDER BY max_score DESC
LIMIT 200
""", 'parent_contrib.csv')

# Summary json
summary = dict(con.execute("SELECT COUNT(*) sequences FROM sequences").fetchone())
summary.update(dict(con.execute("SELECT COUNT(*) metrics FROM metrics").fetchone()))
summary.update(dict(con.execute("SELECT COUNT(*) rounds FROM rounds").fetchone()))
summary['best'] = dict(con.execute("SELECT source_round, best_score, best_ptm, best_plddt, best_chromo FROM sequences WHERE best_score BETWEEN 0 AND 1 ORDER BY best_score DESC LIMIT 1").fetchone())
summary['exports'] = counts
with open(OUT/'summary.json','w',encoding='utf-8') as f: json.dump(summary,f,indent=2,ensure_ascii=False)
print(json.dumps(summary, ensure_ascii=False, indent=2))
con.close()
