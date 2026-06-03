import sys,subprocess,os,numpy as np,pandas as pd
ns=sys.argv[1]; DMD='/home/hohndor/nesy-genome/tools/diamond'
os.chdir('/home/hohndor/nesy-genome/deepgozero')
terms=pd.read_pickle(f'data/{ns}/terms.pkl')['gos'].values.flatten(); td={t:i for i,t in enumerate(terms)}; nt=len(terms)
# ---- Naive: train term frequency (propagated) ----
tr=pd.read_pickle(f'data/{ns}/train_data.pkl')
freq=np.zeros(nt)
for r in tr.itertuples():
    for go in r.prop_annotations:
        j=td.get(go)
        if j is not None: freq[j]+=1
freq/=len(tr)
# ---- similarity: diamond test vs train ----
trf=f'/tmp/{ns}c_train.fa'
with open(trf,'w') as fo:
    for r in tr.itertuples(): fo.write(f'>{r.proteins}\n{r.sequences}\n')
subprocess.run(f'{DMD} makedb --in {trf} -d /tmp/{ns}c_train -p 4',shell=True,capture_output=True)
m8=f'/tmp/{ns}c.m8'
subprocess.run(f'{DMD} blastp -q data/{ns}/test_data.fa -d /tmp/{ns}c_train -o {m8} -p 4 --max-target-seqs 5 --quiet',shell=True)
maxid={}
for line in open(m8):
    c=line.split('\t')
    if c[0]==c[1]: continue
    maxid[c[0]]=max(maxid.get(c[0],0),float(c[2]))
def buck(v): return '<30' if v<30 else '30-50' if v<50 else '50-70' if v<70 else '70-90' if v<90 else '90-100'
def fmax(P,G):
    best=0; gc=G.sum(1)
    for t in range(1,101):
        thr=t/100; pred=(P>=thr); pc=pred.sum(1); tp=(pred*G).sum(1); has=pc>0
        if has.sum()==0: continue
        prec=(tp[has]/pc[has]).mean(); rec=(tp/np.maximum(gc,1)).mean()
        if prec+rec>0: best=max(best,2*prec*rec/(prec+rec))
    return best
# reference df for gold + proteins + blast
ref=pd.read_pickle(f'data/{ns}/predictions_deepgozero_blast.pkl')
prot=ref['proteins'].values
G=np.zeros((len(ref),nt),dtype=np.float32)
for i,r in enumerate(ref.itertuples()):
    for go in r.prop_annotations:
        j=td.get(go)
        if j is not None: G[i,j]=1
bk=np.array([buck(maxid.get(p,0)) for p in prot])
order=['<30','30-50','50-70','70-90','90-100']
# predictor score matrices
preds={}
preds['Naive']=np.tile(freq,(len(ref),1)).astype(np.float32)
preds['BLAST']=np.stack(ref['blast_preds'].values).astype(np.float32)
for base in ['mlp','deepgocnn','proteinfer_raw','deepgozero']:
    fn=f'data/{ns}/predictions_{base}.pkl'
    if os.path.exists(fn):
        d=pd.read_pickle(fn); preds[base]=np.stack(d['preds'].values).astype(np.float32)
print(f'=== {ns}: CAFA-style Fmax by max-identity-to-train (n={[int((bk==b).sum()) for b in order]}) ===')
print(f'{"predictor":16}'+''.join(f'{b:>9}' for b in order)+f'{"ALL":>9}')
for name,P in preds.items():
    row=f'{name:16}'
    for b in order:
        idx=np.where(bk==b)[0]
        row+= f'{fmax(P[idx],G[idx]):9.3f}' if len(idx)>20 else f'{"-":>9}'
    row+=f'{fmax(P,G):9.3f}'
    print(row)
