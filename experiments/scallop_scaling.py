import sys,time,random,itertools
NORM='/home/hohndor/nesy-genome/deepgozero/data/go.norm'
import scallopy
def nf1(ns):
    terms={x.strip().replace('_',':') for x in open(f'/home/hohndor/nesy-genome/deepgozero/data/{ns}/terms.txt') if x.strip()}
    E=[]
    for line in open(NORM):
        s=line.strip().replace('_',':')
        if ' SubClassOf ' not in s: continue
        l,r=s.split(' SubClassOf ')
        if ' and ' in l or ' some ' in l or ' some ' in r: continue
        if l in terms and r in terms: E.append((l,r))
    return E
def scallop_holds(nodes,edges,pr,k):
    ctx=scallopy.ScallopContext(provenance='topkproofs',k=k)
    ctx.add_relation('raw',(str,)); ctx.add_facts('raw',[(pr[v],(v,)) for v in nodes])
    ctx.add_relation('sub',(str,str)); ctx.add_facts('sub',[(None,(c,p)) for c,p in edges])
    ctx.add_rule('holds(T) = raw(T)')
    ctx.add_rule('holds(P) = sub(C,P), holds(C)')
    ctx.run()
    return {t[1][0]:t[0] for t in ctx.relation('holds')}
def exact_holds(nodes,edges,pr):
    # exact P(holds(t)) = prob t in upward closure of sampled raw set; brute force over raw subsets
    par={v:[] for v in nodes}
    for c,p in edges: par.setdefault(p,[]).append(c)  # p holds if any child c holds
    idx={v:i for i,v in enumerate(nodes)}; n=len(nodes); acc={v:0.0 for v in nodes}
    for bits in itertools.product([0,1],repeat=n):
        raw={v:bits[idx[v]] for v in nodes}
        w=1.0
        for v in nodes: w*= pr[v] if raw[v] else (1-pr[v])
        # upward closure
        holds=dict(raw); changed=True
        while changed:
            changed=False
            for c,p in edges:
                if holds[c] and not holds[p]: holds[p]=1; changed=True
        for v in nodes:
            if holds[v]: acc[v]+=w
    return acc
# --- diamond sanity ---
nodes=['a','b','c','d']; edges=[('d','b'),('d','c'),('b','a'),('c','a')]; pr={v:0.5 for v in nodes}
ex=exact_holds(nodes,edges,pr)
print('=== diamond P(holds): exact vs Scallop topk ===')
for k in [1,3,10]:
    sc=scallop_holds(nodes,edges,pr,k)
    print(f'  k={k}: a exact={ex["a"]:.4f} scallop={sc.get("a",0):.4f}  (a has 2 reconvergent paths from d)')
# --- GO scaling + drift ---
ns=sys.argv[1] if len(sys.argv)>1 else 'cc'
E=nf1(ns); allnodes=sorted({x for e in E for x in e})
import networkx as nx
G=nx.DiGraph(); G.add_edges_from(E)
UG=G.to_undirected()
def bfs(n):
    start=max(UG.degree,key=lambda x:x[1])[0]
    from collections import deque; s={start}; q=deque([start]); seen=[start]
    while q and len(seen)<n:
        u=q.popleft()
        for v in UG.neighbors(u):
            if v not in s: s.add(v); seen.append(v); q.append(v)
            if len(seen)>=n: break
    sub=set(seen[:n]); return list(sub),[(c,p) for c,p in E if c in sub and p in sub]
random.seed(1)
print(f'=== Scallop scaling on {ns} (k=3) ===')
print(f'{"nodes":>7}{"edges":>7}{"time_s":>9}{"maxdrift_vs_exact":>18}')
for n in [30,60,120,250,500,1000,3000,len(allnodes)]:
    if n>len(allnodes): break
    nodes,edges=bfs(n); pr={v:random.uniform(0.1,0.9) for v in nodes}
    t0=time.time()
    try:
        sc=scallop_holds(nodes,edges,pr,3); dt=time.time()-t0
        drift='-'
        if len(nodes)<=16:
            ex=exact_holds(nodes,edges,pr); drift=f'{max(abs(sc.get(v,0)-ex[v]) for v in nodes):.4f}'
        print(f'{len(nodes):>7}{len(edges):>7}{dt:>9.1f}{drift:>18}',flush=True)
    except Exception as e:
        print(f'{len(nodes):>7}{len(edges):>7}   FAIL {str(e)[:40]}',flush=True); break
