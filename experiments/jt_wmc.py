import sys,itertools,random,numpy as np,networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in

# ---- factor = (vars tuple, np.array over 2^len in row-major of those vars) ----
def prior_factor(v,p): return ((v,), np.array([1.0-p,p]))
def edge_factor(c,p):  # child c => parent p : forbid c=1,p=0
    # order vars (min,max); table indexed [val_lo, val_hi]
    a,b=(c,p) if c<p else (p,c)
    T=np.ones((2,2))
    # set the (c=1,p=0) cell to 0
    cell=[0,0]; cell[0 if c==a else 1]=1; cell[1 if p==b else 0]=0
    T[cell[0],cell[1]]=0.0
    return ((a,b),T)

def fmul(f1,f2):
    v1,t1=f1; v2,t2=f2; vs=tuple(sorted(set(v1)|set(v2)))
    idx=lambda vv:[vs.index(x) for x in vv]
    r1=t1.reshape([2]*len(v1)); r2=t2.reshape([2]*len(v2))
    b1=np.ones([2]*len(vs)); b2=np.ones([2]*len(vs))
    b1=np.moveaxis(np.broadcast_to(r1.reshape(t1.shape+(1,)*(len(vs)-len(v1))),[2]*len(vs)),range(len(v1)),idx(v1))
    b2=np.moveaxis(np.broadcast_to(r2.reshape(t2.shape+(1,)*(len(vs)-len(v2))),[2]*len(vs)),range(len(v2)),idx(v2))
    return (vs,b1*b2)
def fmarg(f,keep):
    vs,t=f; keep=tuple(sorted(keep)); ax=tuple(i for i,x in enumerate(vs) if x not in keep)
    r=t.reshape([2]*len(vs)); s=r.sum(axis=ax) if ax else r
    return (keep,np.asarray(s).reshape([2]*len(keep)))

def brute_wmc(nodes,edges,priors):
    nodes=list(nodes); Z=0.0; num={n:0.0 for n in nodes}
    for bits in itertools.product([0,1],repeat=len(nodes)):
        a=dict(zip(nodes,bits))
        if any(a[c]==1 and a[p]==0 for c,p in edges): continue
        w=1.0
        for n in nodes: w*= priors[n] if a[n] else 1-priors[n]
        Z+=w
        for n in nodes:
            if a[n]: num[n]+=w
    return {n:num[n]/Z for n in nodes}

def jt_wmc(nodes,edges,priors):
    G=nx.Graph(); G.add_nodes_from(nodes)
    for c,p in edges: G.add_edge(c,p)
    tw,TD=treewidth_min_fill_in(G)   # TD: tree of frozenset bags
    bags=list(TD.nodes()); 
    # assign factors to a bag containing their scope
    fac={b:[] for b in bags}
    def host(scope):
        for b in bags:
            if set(scope)<=set(b): return b
        return None
    for n in nodes:
        b=host((n,)); fac[b].append(prior_factor(n,priors[n]))
    for c,p in edges:
        b=host((c,p))
        if b is None:  # safety: shouldn't happen for a valid TD of G
            b=host((c,)); 
        fac[b].append(edge_factor(c,p))
    # initial bag potential
    def bagpot(b):
        vs=tuple(sorted(b)); pot=(vs,np.ones([2]*len(vs)))
        for f in fac[b]: pot=fmul(pot,f)
        return pot
    pot={b:bagpot(b) for b in bags}
    T=nx.Graph(TD)
    root=bags[0]; order=list(nx.dfs_postorder_nodes(T,root))
    msg={}
    # collect (leaves->root)
    parent={root:None}
    for u,v in nx.dfs_edges(T,root): parent[v]=u
    for b in order:
        pa=parent[b]
        if pa is None: continue
        sep=tuple(sorted(set(b)&set(pa)))
        m=pot[b]
        for ch in T.neighbors(b):
            if ch!=pa and (ch,b) in msg: m=fmul(m,msg[(ch,b)])
        msg[(b,pa)]=fmarg(m,sep)
    # distribute (root->leaves)
    for b in [root]+[v for u,v in nx.dfs_edges(T,root)]:
        for ch in T.neighbors(b):
            if parent[ch]==b:
                sep=tuple(sorted(set(b)&set(ch)))
                m=pot[b]
                for nb in T.neighbors(b):
                    if nb!=ch and (nb,b) in msg: m=fmul(m,msg[(nb,b)])
                msg[(b,ch)]=fmarg(m,sep)
    # bag beliefs -> node marginals
    marg={}
    for b in bags:
        bel=pot[b]
        for nb in T.neighbors(b):
            if (nb,b) in msg: bel=fmul(bel,msg[(nb,b)])
        vs,t=bel; Z=t.sum()
        for n in b:
            if n in marg: continue
            f=fmarg((vs,t),(n,)); marg[n]=float(f[1][1]/Z)
    return marg,tw

def loopy_bp(nodes,edges,priors,iters=200,damp=0.3):
    # sum-product on the factor graph (edge factors + unary priors)
    facs=[edge_factor(c,p) for c,p in edges]
    # messages var->factor and factor->var as 2-vectors
    nfac=len(facs); 
    v2f={}; f2v={}
    for fi,(vs,_) in enumerate(facs):
        for v in vs: v2f[(v,fi)]=np.array([1.0,1.0]); f2v[(fi,v)]=np.array([1.0,1.0])
    var_facs={n:[] for n in nodes}
    for fi,(vs,_) in enumerate(facs):
        for v in vs: var_facs[v].append(fi)
    for _ in range(iters):
        # var->factor
        for v in nodes:
            base=np.array([1-priors[v],priors[v]])
            for fi in var_facs[v]:
                m=base.copy()
                for fj in var_facs[v]:
                    if fj!=fi: m=m*f2v[(fj,v)]
                m=m/ (m.sum() or 1)
                v2f[(v,fi)]=damp*v2f[(v,fi)]+(1-damp)*m
        # factor->var
        for fi,(vs,T) in enumerate(facs):
            for vi,v in enumerate(vs):
                other=vs[1-vi]
                inm=v2f[(other,fi)]
                M=T if vi==0 else T.T
                out=M@inm
                out=out/(out.sum() or 1)
                f2v[(fi,v)]=damp*f2v[(fi,v)]+(1-damp)*out
    marg={}
    for v in nodes:
        b=np.array([1-priors[v],priors[v]])
        for fi in var_facs[v]: b=b*f2v[(fi,v)]
        b=b/(b.sum() or 1); marg[v]=float(b[1])
    return marg

if __name__=='__main__':
    # UNIT TEST: diamond  d=>b,d=>c,b=>a,c=>a ; all priors 1/2  -> exact P(a)=5/6, BP=7/10
    nodes=['a','b','c','d']; edges=[('d','b'),('d','c'),('b','a'),('c','a')]
    pr={n:0.5 for n in nodes}
    bf=brute_wmc(nodes,edges,pr); jt,tw=jt_wmc(nodes,edges,pr); bp=loopy_bp(nodes,edges,pr)
    print('UNIT diamond: brute a=%.4f jt a=%.4f (tw=%d) bp a=%.4f  [expect 0.8333/0.8333/0.7000]'%(bf['a'],jt['a'],tw,bp['a']))
    assert abs(bf['a']-5/6)<1e-9 and abs(jt['a']-5/6)<1e-6, 'JT/brute mismatch'
    print('JT==brute on diamond OK')
