import sys,random,numpy as np,networkx as nx,time
sys.path.insert(0,'/home/hohndor/nesy-genome/deepgozero')
from jt_wmc import jt_wmc,loopy_bp,brute_wmc
NORM='/home/hohndor/nesy-genome/deepgozero/data/go.norm'
def parse_nf1():
    e=[]
    for line in open(NORM):
        s=line.strip().replace('_',':')
        if ' SubClassOf ' not in s: continue
        l,r=s.split(' SubClassOf ')
        if ' and ' in l or ' some ' in l or ' some ' in r: continue
        e.append((l,r))
    return e
ns=sys.argv[1] if len(sys.argv)>1 else 'cc'
terms={x.strip().replace('_',':') for x in open(f'/home/hohndor/nesy-genome/deepgozero/data/{ns}/terms.txt') if x.strip()}
E=[(c,p) for c,p in parse_nf1() if c in terms and p in terms]
V=sorted({x for e in E for x in e})
random.seed(7); pr={v:random.uniform(0.05,0.95) for v in V}
print(f'{ns}: |V|={len(V)} |E|={len(E)}')
G=nx.Graph(); G.add_nodes_from(V); G.add_edges_from(E)
comps=list(nx.connected_components(G)); print('components:',len(comps),'sizes',sorted([len(c) for c in comps],reverse=True)[:5])
t0=time.time(); jt={}; tws=[]
for comp in comps:
    sub=[e for e in E if e[0] in comp and e[1] in comp]
    m,tw=jt_wmc(sorted(comp),sub,{v:pr[v] for v in comp}); jt.update(m); tws.append(tw)
print('JT exact done in %.1fs, max component treewidth=%d'%(time.time()-t0,max(tws)))
# validate JT vs brute on small components (<=18 nodes)
import itertools
for comp in comps:
    if len(comp)<=16:
        sub=[e for e in E if e[0] in comp and e[1] in comp]
        bf=brute_wmc(sorted(comp),sub,{v:pr[v] for v in comp})
        err=max(abs(bf[v]-jt[v]) for v in comp)
        print('  validate component size %d: max|JT-brute|=%.2e'%(len(comp),err)); break
t0=time.time(); bp=loopy_bp(V,E,pr,iters=100,damp=0.4); print('loopy BP done in %.1fs'%(time.time()-t0))
errs=np.array([abs(jt[v]-bp[v]) for v in V])
indeg={v:0 for v in V}
for c,p in E: indeg[c]+=1   # child->parent; multi-parent = indeg(child as subclass)>1
mp=np.array([1 if sum(1 for c,p in E if c==v)>1 else 0 for v in V])
print('JT vs BP error: mean=%.4f  max=%.4f  p95=%.4f'%(errs.mean(),errs.max(),np.percentile(errs,95)))
print('mean error at multi-parent nodes=%.4f vs single-parent=%.4f'%(errs[mp==1].mean(),errs[mp==0].mean()))
print('frac terms |err|>0.01: %.3f  |err|>0.05: %.3f'%((errs>0.01).mean(),(errs>0.05).mean()))
np.save(f'/home/hohndor/nesy-genome/deepgozero/{ns}_jtbp_err.npy',errs)
