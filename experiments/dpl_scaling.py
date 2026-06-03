import sys,os,time,random,subprocess,signal,numpy as np,networkx as nx
sys.path.insert(0,'/home/hohndor/nesy-genome/deepgozero')
from networkx.algorithms.approximation import treewidth_min_fill_in
from jt_wmc import jt_wmc
NORM='/home/hohndor/nesy-genome/deepgozero/data/go.norm'
PL='/home/hohndor/nesy-genome/nesy-venv/bin/problog'
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
def bfs_sub(G,n):
    start=max(G.degree,key=lambda x:x[1])[0]
    seen=[start]; from collections import deque; q=deque([start]); s={start}
    while q and len(seen)<n:
        u=q.popleft()
        for v in G.neighbors(u):
            if v not in s: s.add(v); seen.append(v); q.append(v)
            if len(seen)>=n: break
    return G.subgraph(seen[:n]).copy()
def goid(t): return 't'+t.replace('GO:','').replace(':','_')
def run_problog(sub,priors,timeout=300):
    lines=[]
    for v in sub.nodes(): lines.append(f'{priors[v]:.4f}::{goid(v)}.')
    for c,p in sub.edges():
        # orient child->parent by original direction unknown in undirected; use both? we kept (c,p) dir via edge attr
        pass
    # need directed edges: rebuild from stored
    for (c,p) in sub.graph['edges']:
        lines.append(f'viol :- {goid(c)}, \\+ {goid(p)}.')
    lines.append('evidence(viol,false).')
    qn=list(sub.nodes())[:1]
    for v in qn: lines.append(f'query({goid(v)}).')
    prog='\n'.join(lines)
    f=f'/tmp/dpl_{len(sub)}.pl'; open(f,'w').write(prog)
    t0=time.time()
    try:
        r=subprocess.run(f'ulimit -v 16000000; {PL} {f}',shell=True,capture_output=True,text=True,timeout=timeout)
        dt=time.time()-t0
        ok=(goid(qn[0]) in r.stdout)
        return ('OK' if ok else 'FAIL',dt,r.stdout.strip().split(chr(10))[0][:60]+'|'+r.stderr.strip()[-60:])
    except subprocess.TimeoutExpired:
        return ('TIMEOUT',timeout,'')
ns=sys.argv[1] if len(sys.argv)>1 else 'cc'
E=nf1(ns); G=nx.Graph(); G.add_edges_from(E)
random.seed(1)
print(f'=== DeepProbLog/ProbLog true-path WMC scaling on {ns} ===')
print(f'{"nodes":>7}{"edges":>7}{"tw":>5}{"problog":>10}{"time_s":>9}')
sizes=[30,100,250,400,600,800,1200,len(G)]
for n in sizes:
    sub=bfs_sub(G,n) if n<len(G) else G.copy()
    sub.graph['edges']=[(c,p) for c,p in E if c in sub and p in sub]
    pr={v:random.uniform(0.1,0.9) for v in sub.nodes()}
    tw,_=treewidth_min_fill_in(sub)
    status,dt,msg=run_problog(sub,pr)
    print(f'{sub.number_of_nodes():>7}{sub.number_of_edges():>7}{tw:>5}{status:>10}{dt:>9.1f}  {msg if status!=chr(79)+chr(75) else ""}',flush=True)
    if status=='TIMEOUT' or status=='FAIL': 
        print('  (stop sweep at first failure)'); break
