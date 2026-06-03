import sys,os,numpy as np,pandas as pd
ns=sys.argv[1]; os.chdir('/home/hohndor/nesy-genome/deepgozero')
terms=pd.read_pickle(f'data/{ns}/terms.pkl')['gos'].values.flatten(); td={t:i for i,t in enumerate(terms)}; nt=len(terms)
tr=pd.read_pickle(f'data/{ns}/train_data.pkl')
# seen InterPro domains in training + their frequency
seen={}
for r in tr.itertuples():
    for ip in (r.interpros or []): seen[ip]=seen.get(ip,0)+1
# Naive freq
freq=np.zeros(nt)
for r in tr.itertuples():
    for go in r.prop_annotations:
        j=td.get(go); 
        if j is not None: freq[j]+=1
freq/=len(tr)
ref=pd.read_pickle(f'data/{ns}/predictions_deepgozero_blast.pkl')
G=np.zeros((len(ref),nt),dtype=np.float32)
for i,r in enumerate(ref.itertuples()):
    for go in r.prop_annotations:
        j=td.get(go)
        if j is not None: G[i,j]=1
# domain-novelty bucket per test protein
def dbuck(ips):
    ips=ips or []
    if len(ips)==0: return 'no-IPR'
    shared=[ip for ip in ips if ip in seen]
    if len(shared)==0: return 'all-novel'
    frac=len(shared)/len(ips)
    return 'partial' if frac<0.999 else 'all-seen'
bk=np.array([dbuck(r.interpros) for r in ref.itertuples()])
order=['no-IPR','all-novel','partial','all-seen']
def fmax(P,Gm):
    best=0; gc=Gm.sum(1)
    for t in range(1,101):
        thr=t/100; pred=(P>=thr); pc=pred.sum(1); tp=(pred*Gm).sum(1); has=pc>0
        if has.sum()==0: continue
        prec=(tp[has]/pc[has]).mean(); rec=(tp/np.maximum(gc,1)).mean()
        if prec+rec>0: best=max(best,2*prec*rec/(prec+rec))
    return best
preds={'Naive':np.tile(freq,(len(ref),1)).astype(np.float32),
       'BLAST':np.stack(ref['blast_preds'].values).astype(np.float32)}
for base in ['mlp','deepgocnn','proteinfer_raw','deepgozero']:
    fn=f'data/{ns}/predictions_{base}.pkl'
    if os.path.exists(fn): preds[base]=np.stack(pd.read_pickle(fn)['preds'].values).astype(np.float32)
print(f'=== {ns}: Fmax by InterPro-DOMAIN novelty (n={[int((bk==b).sum()) for b in order]}) ===')
print(f'{"predictor":16}'+''.join(f'{b:>10}' for b in order)+f'{"ALL":>9}')
for name,P in preds.items():
    row=f'{name:16}'
    for b in order:
        idx=np.where(bk==b)[0]
        row+= f'{fmax(P[idx],G[idx]):10.3f}' if len(idx)>20 else f'{"-":>10}'
    row+=f'{fmax(P,G):9.3f}'
    print(row)
