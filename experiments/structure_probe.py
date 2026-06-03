import pickle,sys
from collections import defaultdict
NORM='/home/hohndor/nesy-genome/deepgozero/data/go.norm'
def parse(path):
    nf1=[]; nf2=[]; nf3=[]; nf4=[]
    for line in open(path):
        s=line.strip().replace('_',':')
        if ' SubClassOf ' not in s: continue
        l,r=s.split(' SubClassOf ')
        if ' and ' in l:
            a,b=l.split(' and '); nf2.append((a,b,r))
        elif ' some ' in l:
            R,C=l.split(' some '); nf3.append(('E:'+R+':'+C, r))
        elif ' some ' in r:
            R,D=r.split(' some '); nf4.append((l,'E:'+R+':'+D))
        else:
            nf1.append((l,r))
    return nf1,nf2,nf3,nf4
nf1,nf2,nf3,nf4=parse(NORM)
print('axioms: nf1=%d nf2=%d nf3=%d nf4=%d'%(len(nf1),len(nf2),len(nf3),len(nf4)))

# ---- union-find for components + circuit rank on the PRIMAL graph of a factor set ----
class UF:
    def __init__(s): s.p={}
    def find(s,x):
        s.p.setdefault(x,x)
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def union(s,a,b):
        ra,rb=s.find(a),s.find(b)
        if ra!=rb: s.p[ra]=rb

def primal_stats(factors, name, restrict=None):
    # factors: list of tuples of atoms (each factor = a clique over its atoms)
    edges=set(); V=set()
    for f in factors:
        atoms=[a for a in f if (restrict is None or a in restrict)]
        if len(atoms)<2: 
            for a in atoms: V.add(a)
            continue
        for a in atoms: V.add(a)
        for i in range(len(atoms)):
            for j in range(i+1,len(atoms)):
                e=tuple(sorted((atoms[i],atoms[j]))); edges.add(e)
    uf=UF()
    for v in V: uf.find(v)
    for a,b in edges: uf.union(a,b)
    comps=len({uf.find(v) for v in V}) if V else 0
    Vn=len(V); En=len(edges)
    circuit_rank=En - Vn + comps   # independent cycles; 0 iff forest
    # multi-parent (nf1 indegree) within restrict
    indeg=defaultdict(int)
    for c,p in nf1:
        if restrict is None or (c in restrict and p in restrict): indeg[c]+=1
    mp=sum(1 for k,v in indeg.items() if v>1)
    maxin=max(indeg.values()) if indeg else 0
    print('[%s] |V|=%d |E|=%d comps=%d circuit_rank=%d  forest=%s  multiparent=%d maxindeg=%d'%(
        name,Vn,En,comps,circuit_rank,circuit_rank==0,mp,maxin))
    return V,edges

# nf1 hierarchy graph (what the true-path WMC in go_wmc_eval uses)
nf1_factors=[(c,p) for c,p in nf1]
# full theory primal (all normal forms)
full_factors=[(c,p) for c,p in nf1]+[(a,b,r) for a,b,r in nf2]+[(n,r) for n,r in nf3]+[(l,n) for l,n in nf4]

print('=== GLOBAL ===')
primal_stats(nf1_factors,'nf1 hierarchy (global)')
primal_stats(full_factors,'full EL theory (global)')

for ns in ['cc','mf','bp']:
    terms={l.strip().replace('_',':') for l in open(f'/home/hohndor/nesy-genome/deepgozero/data/{ns}/terms.txt') if l.strip()}
    print(f'=== namespace {ns}: {len(terms)} predicted terms ===')
    primal_stats(nf1_factors,f'{ns} nf1 hierarchy',restrict=terms)
    Vf,Ef=primal_stats(full_factors,f'{ns} full EL theory',restrict=terms)
    # min-degree treewidth UPPER bound on the induced full-theory primal graph
    adj={v:set() for v in Vf}
    for a,b in Ef: adj[a].add(b); adj[b].add(a)
    import copy; A={v:set(n) for v,n in adj.items()}; tw=0
    while A:
        v=min(A,key=lambda x:len(A[x])); nb=A[v]; tw=max(tw,len(nb))
        for u in nb:
            A[u]|=(nb-{u}); A[u].discard(v)
        del A[v]
        for u in list(A): A[u].discard(v)
    print('   [%s] min-degree treewidth upper bound = %d'%(ns,tw))
