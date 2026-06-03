import sys,subprocess,os,numpy as np,pandas as pd
ns=sys.argv[1]
PY='/home/hohndor/nesy-genome/moose/.venv/bin/python'
WD=f'/home/hohndor/nesy-genome/deepgozero'
os.chdir(WD)
# 1) train fasta + diamond db
tr=pd.read_pickle(f'data/{ns}/train_data.pkl')
trf=f'/tmp/{ns}_train.fa'
with open(trf,'w') as fo:
    for r in tr.itertuples(): fo.write(f'>{r.proteins}\n{r.sequences}\n')
subprocess.run(f'/home/hohndor/nesy-genome/tools/diamond makedb --in {trf} -d /tmp/{ns}_train -p 4',shell=True,capture_output=True)
# 2) diamond test vs train -> max pident per test protein
m8=f'/tmp/{ns}_test_vs_train.m8'
subprocess.run(f'/home/hohndor/nesy-genome/tools/diamond blastp -q data/{ns}/test_data.fa -d /tmp/{ns}_train -o {m8} -p 4 --max-target-seqs 5 --quiet',shell=True)
maxid={}
for line in open(m8):
    c=line.split('\t'); q=c[0]; pid=float(c[2])
    if q==c[1]: continue
    maxid[q]=max(maxid.get(q,0),pid)
def buck(v):
    if v<30:return '<30'
    if v<50:return '30-50'
    if v<70:return '50-70'
    if v<90:return '70-90'
    return '90-100'
def fmax(P,G):
    best=0
    gc=G.sum(1)
    for t in range(1,101):
        thr=t/100; pred=(P>=thr); pc=pred.sum(1); tp=(pred*G).sum(1); has=pc>0
        if has.sum()==0: continue
        prec=(tp[has]/pc[has]).mean(); rec=(tp/np.maximum(gc,1)).mean()
        if prec+rec>0: best=max(best,2*prec*rec/(prec+rec))
    return best
terms=pd.read_pickle(f'data/{ns}/terms.pkl')['gos'].values.flatten(); td={t:i for i,t in enumerate(terms)}
order=['<30','30-50','50-70','70-90','90-100']
print(f'=== {ns} memorization: Fmax by max-identity-to-train ===')
print(f'{"base":16}'+''.join(f'{b:>9}' for b in order)+'   n_per_bucket')
for base in ['mlp','deepgocnn','proteinfer_raw','deepgozero']:
    fn=f'data/{ns}/predictions_{base}.pkl'
    if not os.path.exists(fn): continue
    df=pd.read_pickle(fn)
    prot=df['proteins'].values; P=np.stack(df['preds'].values).astype(np.float32)
    G=np.zeros((len(df),len(terms)),dtype=np.float32)
    for i,r in enumerate(df.itertuples()):
        for go in r.prop_annotations:
            j=td.get(go)
            if j is not None: G[i,j]=1
    bk=np.array([buck(maxid.get(p,0)) for p in prot])
    row=f'{base:16}'; ns_cnt=[]
    for b in order:
        idx=np.where(bk==b)[0]; ns_cnt.append(len(idx))
        row+= f'{fmax(P[idx],G[idx]):9.3f}' if len(idx)>20 else f'{"-":>9}'
    print(row+'   '+str(ns_cnt))
