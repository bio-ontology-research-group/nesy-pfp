import sys,random,numpy as np,networkx as nx,time
sys.path.insert(0,'/home/hohndor/nesy-genome/deepgozero')
from jt_wmc import jt_wmc,loopy_bp
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
for ns in ['mf','bp']:
    terms={x.strip().replace('_',':') for x in open(f'/home/hohndor/nesy-genome/deepgozero/data/{ns}/terms.txt') if x.strip()}
    E=[(c,p) for c,p in ALL if c in terms and p in terms]
    V=sorted({x for e in E for x in e})
    G=nx.Graph(); G.add_nodes_from(V); G.add_edges_from(E)
    comps=list(nx.connected_components(G))
    TRACT=18; tract_terms=0; tws=[]
    for comp in comps:
        sub=G.subgraph(comp)
        tw,_=treewidth_min_fill_in(sub); tws.append(tw)
        if tw<=TRACT: tract_terms+=len(comp)
    print(f'{ns}: |V|={len(V)} comps={len(comps)} max_comp_tw={max(tws)} | terms in exactly-tractable(tw<={TRACT}) components: {tract_terms}/{len(V)} = {tract_terms/len(V):.2%}')
