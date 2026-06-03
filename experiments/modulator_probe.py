import sys,numpy as np,networkx as nx,time
from networkx.algorithms.approximation import treewidth_min_fill_in
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
ALL=parse_nf1()
def tw_of(G):
    tws=[]
    for comp in nx.connected_components(G):
        if len(comp)<3: continue
        tw,_=treewidth_min_fill_in(G.subgraph(comp)); tws.append(tw)
    return max(tws) if tws else 0
for ns in ['mf','bp']:
    terms={x.strip().replace('_',':') for x in open(f'/home/hohndor/nesy-genome/deepgozero/data/{ns}/terms.txt') if x.strip()}
    E=[(c,p) for c,p in ALL if c in terms and p in terms]
    G=nx.Graph(); G.add_edges_from(E)
    N=G.number_of_nodes(); TARGET=14; cut=[]
    print(f'{ns}: N={N} start tw={tw_of(G)}',flush=True)
    t0=time.time()
    while True:
        tw=tw_of(G)
        if tw<=TARGET:
            print(f'  REACHED tw<={TARGET} after cutting {len(cut)} hubs ({len(cut)/N:.2%}); remaining terms exactly-computable: {G.number_of_nodes()}/{N} = {G.number_of_nodes()/N:.2%}  [{time.time()-t0:.0f}s]',flush=True); break
        # remove top-degree nodes (reconvergence hubs) in a batch
        k=max(5,int(0.01*G.number_of_nodes()))
        hubs=[v for v,_ in sorted(G.degree(),key=lambda x:-x[1])[:k]]
        G.remove_nodes_from(hubs); cut+=hubs
        print(f'  cut {len(cut)} (tw was {tw}); now N={G.number_of_nodes()}',flush=True)
        if len(cut)>0.3*N: print('  >30% cut, abort'); break
