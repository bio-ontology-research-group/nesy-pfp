import argparse,numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score
LO=1e-12
def parse_norm_full(path):
    nf1,nf2,nf3,nf4=[],[],[],[]
    for line in open(path):
        s=line.strip().replace('_',':')
        if ' SubClassOf ' not in s: continue
        l,r=s.split(' SubClassOf ')
        if ' and ' in l: a,b=l.split(' and '); nf2.append((a,b,r))
        elif ' some ' in l: R,C=l.split(' some '); nf3.append((R,C,r))
        elif ' some ' in r: R,D=r.split(' some '); nf4.append((l,R,D))
        else: nf1.append((l,r))
    return nf1,nf2,nf3,nf4
def noisy_or(cands):
    acc=np.zeros_like(cands[0])
    for c in cands: acc=acc+np.log1p(-np.clip(c,0,1-1e-9))
    return 1-np.exp(acc)
ap=argparse.ArgumentParser()
ap.add_argument('--ont',default='cc'); ap.add_argument('--base',default='mlp_raw')
ap.add_argument('--minpos',type=int,default=10)
a=ap.parse_args()
nf1,nf2,nf3,nf4=parse_norm_full('data/go.norm')
terms=pd.read_pickle(f'data/{a.ont}/terms.pkl')['gos'].values.flatten()
tdict={t:i for i,t in enumerate(terms)}
df=pd.read_pickle(f'data/{a.ont}/predictions_{a.base}.pkl')
P=np.stack(df['preds'].values).astype(np.float64)
gold=[set(r.prop_annotations) for r in df.itertuples()]
# definitions grounded in in-vocab terms
def2={}; 
for x,y,e in nf2:
    if x in tdict and y in tdict and e in tdict: def2.setdefault(e,[]).append((tdict[x],tdict[y]))
def3={}
for R,C,d in nf3:
    if C in tdict and d in tdict: def3.setdefault(d,[]).append(tdict[C])
defined=set(def2)|set(def3)
print(f'{a.ont} base={a.base}: terms={len(terms)} proteins={P.shape[0]} defined-terms(grounded)={len(defined)}')
rec_aucs=[]; dir_aucs=[]; n=0
for e in defined:
    j=tdict[e]; y=np.array([1 if e in g else 0 for g in gold])
    if y.sum()<a.minpos or y.sum()==len(y): continue
    cands=[]
    for (ia,ib) in def2.get(e,[]): cands.append(P[:,ia]*P[:,ib])
    for ic in def3.get(e,[]): cands.append(P[:,ic])
    if not cands: continue
    recon=noisy_or(cands)
    rec_aucs.append(roc_auc_score(y,recon)); dir_aucs.append(roc_auc_score(y,P[:,j])); n+=1
print(f'evaluable defined terms (>= {a.minpos} pos): {n}')
print(f'macro WMC-reconstruct AUC = {np.mean(rec_aucs):.4f}')
print(f'macro base-direct     AUC = {np.mean(dir_aucs):.4f}')
print(f'reconstruct >= direct on {np.mean(np.array(rec_aucs)>=np.array(dir_aucs)):.1%} of terms; recon AUC>0.7 on {np.mean(np.array(rec_aucs)>0.7):.1%}')
np.savez(f'data/{a.ont}_{a.base}_zeroscale.npz',rec=rec_aucs,dir=dir_aucs)
