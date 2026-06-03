import sys,numpy as np,pandas as pd,torch as th
from torch import nn
from sklearn.metrics import roc_auc_score
th.set_num_threads(4)
ont=sys.argv[1]; minpos=int(sys.argv[2]) if len(sys.argv)>2 else 10
dr='data'
# ---------- WMC closure (go_zero_wmc) over an in-vocab base ----------
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
nf1,nf2,nf3,nf4=parse_norm_full(f'{dr}/go.norm')
def2={}; def3={}; sub1={}
for a,b,e in nf2: def2.setdefault(e,[]).append((a,b))
for R,C,d in nf3: def3.setdefault(d,[]).append(('E',R,C))
for a,b in nf1: sub1.setdefault(b,[]).append(a)
terms=pd.read_pickle(f'{dr}/{ont}/terms.pkl')['gos'].values.flatten(); tdict={t:i for i,t in enumerate(terms)}
base=pd.read_pickle(f'{dr}/{ont}/predictions_deepgozero.pkl')
Pin=np.stack(base['preds'].values).astype(np.float64); nprot=Pin.shape[0]; ZERO=np.zeros(nprot)
def noisy_or(c):
    acc=np.zeros_like(c[0])
    for x in c: acc=acc+np.log1p(-np.clip(x,0,1-1e-9))
    return 1-np.exp(acc)
memo={}
def score(node):
    if node in memo: return memo[node]
    if isinstance(node,str) and node in tdict: memo[node]=Pin[:,tdict[node]]; return memo[node]
    memo[node]=ZERO; cands=[]
    if isinstance(node,tuple): cands.append(score(node[2]))
    for a,b in def2.get(node,[]): cands.append(score(a)*score(b))
    for z in def3.get(node,[]): cands.append(score(z))
    memo[node]=noisy_or(cands) if cands else ZERO
    return memo[node]
def wmc_can(t): return (t in def2) or (t in def3)   # has a grounded definition
# ---------- DeepGOZero predict_zero baseline ----------
class MLPBlock(nn.Module):
    def __init__(self,i,o,layer_norm=True,dropout=0.1,activation=nn.ReLU):
        super().__init__(); self.linear=nn.Linear(i,o); self.activation=activation()
        self.layer_norm=nn.LayerNorm(o) if layer_norm else None
        self.dropout=nn.Dropout(dropout) if dropout else None
    def forward(self,x):
        x=self.activation(self.linear(x))
        if self.layer_norm is not None: x=self.layer_norm(x)
        if self.dropout is not None: x=self.dropout(x)
        return x
class Residual(nn.Module):
    def __init__(self,fn): super().__init__(); self.fn=fn
    def forward(self,x): return x+self.fn(x)
class DGEL(nn.Module):
    def __init__(self,nip,ng,nz,nr,ed=1024,hd=1024):
        super().__init__()
        self.net=nn.Sequential(MLPBlock(nip,hd),Residual(MLPBlock(hd,hd)))
        self.hasFuncIndex=th.LongTensor([nr])
        self.go_embed=nn.Embedding(ng+nz,ed); self.go_norm=nn.BatchNorm1d(ed)
        self.go_rad=nn.Embedding(ng+nz,1); self.rel_embed=nn.Embedding(nr+1,ed)
    def predict_zero(self,f,data):
        x=self.net(f); ge=self.go_embed(data); hf=self.rel_embed(self.hasFuncIndex)
        hg=ge+hf; gr=th.abs(self.go_rad(data).view(1,-1))
        return th.sigmoid(th.matmul(x,hg.T)+gr)
def load_nf_idx(go_file,td):
    rels={}; z={}
    def gi(g):
        if g in td: return td[g]
        if g in z: return z[g]
        z[g]=len(td)+len(z); return z[g]
    def gr(r):
        if r not in rels: rels[r]=len(rels)
        return rels[r]
    for line in open(go_file):
        line=line.strip().replace('_',':')
        if 'SubClassOf' not in line: continue
        L,R=line.split(' SubClassOf ')
        if len(L)==10 and len(R)==10: gi(L);gi(R)
        elif ' and ' in L: a,b=L.split(' and '); gi(a);gi(b);gi(R)
        elif ' some ' in L: rel,c=L.split(' some '); gr(rel);gi(c);gi(R)
        elif ' some ' in R: rel,d=R.split(' some '); gr(rel);gi(d);gi(L)
    return rels,z
zterms=pd.read_pickle(f'{dr}/{ont}/terms_zero_10.pkl')['gos'].values.flatten(); ztd={t:i for i,t in enumerate(zterms)}
iprs=pd.read_pickle(f'{dr}/{ont}/interpros.pkl')['interpros'].values; ipd={v:i for i,v in enumerate(iprs)}
rels,zc=load_nf_idx(f'{dr}/go.norm',ztd)
test=pd.read_pickle(f'{dr}/{ont}/test_data.pkl').reset_index(drop=True)
feats=th.zeros((len(test),len(ipd)),dtype=th.float32)
for i,row in enumerate(test.itertuples()):
    for ip in row.interpros:
        if ip in ipd: feats[i,ipd[ip]]=1
net=DGEL(len(ipd),len(zterms),len(zc),len(rels)); net.load_state_dict(th.load(f'{dr}/{ont}/deepgozero_zero_10.th',map_location='cpu')); net.eval()
gold=[set(r.prop_annotations) for r in test.itertuples()]
# proteins must align between base preds and test_data: both are test set; assume same order
assert len(test)==nprot, f'len test {len(test)} vs base {nprot}'
# eval terms: zero_10 terms with >=minpos positives
wmc_a=[]; dgz_a=[]; both=0; wmc_only_cov=0
for t in zterms:
    y=np.array([1 if t in g else 0 for g in gold])
    if y.sum()<minpos or y.sum()==len(y): continue
    if t in ztd: idx=ztd[t]
    elif t in zc: idx=zc[t]
    else: continue
    with th.no_grad(): sd=net.predict_zero(feats,th.LongTensor([idx])).numpy().ravel()
    dgz=roc_auc_score(y,sd)
    if wmc_can(t):
        memo.clear(); sw=score(t)
        if np.ptp(sw)==0: continue
        wmc=roc_auc_score(y,sw); wmc_a.append(wmc); dgz_a.append(dgz); both+=1
print(f'{ont}: zero_10={len(zterms)} | evaluable(both, >={minpos}pos, WMC-groundable)={both}')
if both:
    wa=np.array(wmc_a); da=np.array(dgz_a)
    print(f'  macro WMC-AUC   = {wa.mean():.4f}')
    print(f'  macro DGZ0-AUC  = {da.mean():.4f}')
    print(f'  WMC > DGZ0 on {(wa>da).mean():.1%} of terms; WMC>=0.7 on {(wa>=0.7).mean():.1%}; DGZ0>=0.7 on {(da>=0.7).mean():.1%}')
    from scipy.stats import wilcoxon
    try: print('  Wilcoxon p =',wilcoxon(wa,da).pvalue)
    except Exception as e: print('  wilcoxon NA',e)
    np.savez(f'data/{ont}_h2h.npz',wmc=wa,dgz=da)
