import sqlite3, json
con=sqlite3.connect(r'D:\workspace\protein-design-atlas\data\design.db')
con.row_factory=sqlite3.Row
print('=== DB stats ===')
for t in ['sequences','metrics','rounds','artifacts','documents','lineage_edges']:
    print(f"  {t}: {con.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0]}")
print('\n=== Top 10 sequences (by best_score) ===')
for r in con.execute("SELECT source_round,best_score,best_ptm,best_plddt,best_chromo,best_parent,length FROM sequences WHERE best_score BETWEEN 0 AND 1 ORDER BY best_score DESC LIMIT 10"):
    print(f"  {r['source_round']:5s} score={r['best_score']:.4f} pTM={r['best_ptm']:.4f} pLDDT={r['best_plddt']:.3f} chromo={r['best_chromo']:.3f} parent={r['best_parent']} len={r['length']}")
print('\n=== Round best scores ===')
for r in con.execute("SELECT round_key,round_number,best_score FROM rounds WHERE best_score BETWEEN 0 AND 1 ORDER BY round_number"):
    print(f"  {r['round_key']:5s} (#{r['round_number']:2d}): {r['best_score']:.4f}")
print('\n=== Score distribution by round (score>0.90) ===')
for r in con.execute("""
SELECT s.source_round, COUNT(*) n, MAX(s.best_score) max_s, AVG(s.best_score) avg_s, AVG(s.best_chromo) avg_c
FROM sequences s WHERE s.best_score>0.90 AND s.best_score<=1
GROUP BY s.source_round ORDER BY s.source_round
"""):
    print(f"  {r['source_round']:5s}: n={r['n']:6d} max={r['max_s']:.4f} avg={r['avg_s']:.4f} avg_chromo={r['avg_c']:.4f}")
print('\n=== Parent contribution (score>0.93) ===')
for r in con.execute("""
SELECT COALESCE(best_parent,'?') parent, source_round, COUNT(*) n, MAX(best_score) max_s, AVG(best_score) avg_s
FROM sequences WHERE best_score>0.93 AND best_score<=1
GROUP BY source_round, best_parent ORDER BY max_s DESC LIMIT 15
"""):
    print(f"  {r['source_round']:5s} {r['parent'][:20]:20s}: n={r['n']:5d} max={r['max_s']:.4f} avg={r['avg_s']:.4f}")
con.close()
