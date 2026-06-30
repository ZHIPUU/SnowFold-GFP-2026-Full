fa = '/root/autodl-tmp/r20/mpnn_out/T=0.2, sample=25, score=0.6229, global_score=0.6383, seq_recovery=0.8291/seqs/T=0.2, sample=25, score=0.6229, global_score=0.6383, seq_recovery=0.8291.fa'
n = 0
with open(fa) as f:
    for line in f:
        if line.startswith('>'):
            n += 1
        elif n in [1,2,3,4,5,6,7,8,9,10,20,50,100,500,1001]:
            print(f'seq[{n}]: {line.strip()[:30]}')
            if n == 1001:
                break
