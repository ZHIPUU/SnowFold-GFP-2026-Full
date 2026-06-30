import os
fa_dir = '/root/autodl-tmp/r20/mpnn_out'
for d in os.listdir(fa_dir):
    seqs_dir = os.path.join(fa_dir, d, 'seqs')
    if not os.path.isdir(seqs_dir):
        continue
    for root, dirs, files in os.walk(seqs_dir):
        for f in files:
            if f.endswith('.fa'):
                fp = os.path.join(root, f)
                n = 0
                parent_count = 0
                with open(fp) as fh:
                    cur_name = ''
                    for line in fh:
                        line = line.rstrip('\n')
                        if line.startswith('>'):
                            n += 1
                            cur_name = line
                            if 'designed_chains' in line and 'sample=' not in line.split(',')[0]:
                                parent_count += 1
                print(f'{d}: {n} sequences (1 parent + {n-1} designs)')

# Also test parse_fa manually
def parse_fa(paths):
    seqs = []
    for p in paths:
        with open(p, encoding='utf-8', errors='replace') as f:
            n, s = '', ''
            for l in f:
                l = l.strip()
                if l.startswith('>'):
                    if s: seqs.append({'name': n, 'seq': s})
                    n, s = l[1:], ''
                else:
                    if l: s += l
            if s: seqs.append({'name': n, 'seq': s})
    return seqs

seqs_dir = '/root/autodl-tmp/r20/mpnn_out/T=0.2, sample=25, score=0.6229, global_score=0.6383, seq_recovery=0.8291/seqs'
files = []
for f in os.listdir(seqs_dir):
    if f.endswith('.fa'):
        files.append(os.path.join(seqs_dir, f))
print(f'\nparse_fa input: {len(files)} file(s)')
parsed = parse_fa(files)
print(f'parse_fa output: {len(parsed)} sequences')
print(f'first parsed: {parsed[0]["name"][:80]}')
print(f'last parsed: {parsed[-1]["name"][:80]}')
