import sqlite3, os
p = r'D:\workspace\protein-design-atlas\data\design.db'
con=sqlite3.connect(p)
cur=con.cursor()
for sql,label in [
    ('UPDATE metrics SET plddt = plddt/100.0 WHERE plddt > 1.5','metrics plddt scaled'),
    ('UPDATE metrics SET chromo = chromo/100.0 WHERE chromo > 1.5','metrics chromo scaled'),
    ('UPDATE sequences SET best_plddt = best_plddt/100.0 WHERE best_plddt > 1.5','seq plddt scaled'),
    ('UPDATE sequences SET best_chromo = best_chromo/100.0 WHERE best_chromo > 1.5','seq chromo scaled'),
]:
    cur.execute(sql)
    print(label, cur.rowcount)
con.commit()
print('max/min chromo', con.execute('select max(best_chromo), min(best_chromo) from sequences where best_chromo is not null').fetchone())
print('max/min plddt', con.execute('select max(best_plddt), min(best_plddt) from sequences where best_plddt is not null').fetchone())
con.close()
