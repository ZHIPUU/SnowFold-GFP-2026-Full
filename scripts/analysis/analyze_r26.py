import json
paths={
'R22':r'D:\workspace\round22\final_6_r22.json',
'R23':r'D:\workspace\round23\final_6_r23.json',
'R26':r'D:\workspace\round26_local\final_6_r26.json',
}
for name,p in paths.items():
    try:
        data=json.load(open(p))
    except Exception as e:
        print(name, 'ERR', e); continue
    print(f'=== {name} ===')
    for i,c in enumerate(data[:6],1):
        print(f"{i}. score={c['score']:.4f} pTM={c['ptm']:.4f} pLDDT={c['plddt']:.3f} chromo={c['chromo']:.3f} parent={c.get('parent','')}")
    print()
r26=json.load(open(paths['R26']))
r22=json.load(open(paths['R22']))
r23=json.load(open(paths['R23']))
print('R26 top1 vs R24 current reference 0.9447:', f"{r26[0]['score']-0.9447:+.4f}")
print('R26 top1 vs R22:', f"{r26[0]['score']-r22[0]['score']:+.4f}")
print('R26/R22 overlap', len({c['seq'] for c in r26} & {c['seq'] for c in r22}))
print('R26/R23 overlap', len({c['seq'] for c in r26} & {c['seq'] for c in r23}))
print('R26 unique passed count', len(json.load(open(r'D:\workspace\round26_local\all_passed.json'))))
