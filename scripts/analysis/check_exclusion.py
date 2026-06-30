import pandas as pd
sub = pd.read_csv(r'D:\workspace\round18\submission_r18.csv')
excl = pd.read_csv(r'D:\生信\2026Protein Design\Exclusion_List.csv')
excl_seqs = set(excl.iloc[:, 0].astype(str))
print(f'Exclusion list: {len(excl_seqs)} seqs')
print()
for _, r in sub.iterrows():
    s = r['Sequence']
    in_excl = s in excl_seqs
    checks = f"Seq {r['Seq_ID']}: M={s[0]=='M'} Len={len(s)} AA={all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in s)} Excl={'BAD' if in_excl else 'OK'}"
    print(checks)
print()
print('All passed!' if not any(s in excl_seqs for s in sub['Sequence']) else 'WARNING: Found exclusion matches!')
