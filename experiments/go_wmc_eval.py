#!/usr/bin/env python
"""GO-hierarchy WMC benchmark: raw vs max-propagation vs soft-WMC closure.

Symbolic structure = the EL theory `go.norm` (DeepGOZero's normalised GO):
  nf1  A SubClassOf B            (A -> B)            -- subsumption (is_a/part_of/regulates)
  nf2  A and B SubClassOf E      (A & B -> E)        -- GO-Plus logical definitions (conjunction)

Max-propagation only exploits nf1 upward (parent = max over descendants). It
*cannot* infer E from the joint presence of its definitional parts A and B (nf2).
The soft-WMC closure conditions the independent per-term marginals on the EL theory
via a monotone noisy-OR fixpoint (the same log-space scatter/gather scheme used for
metabolic producibility), so nf2 adds recall the heuristic structurally misses.

Run with a torch venv:  python go_wmc_eval.py --ont cc --model mlp_raw
"""
import argparse, numpy as np, pandas as pd, torch as th
from sklearn.metrics import average_precision_score

LO = 1e-12


def parse_norm(path, terms_dict):
    """nf1 = [(x,y): x->y], nf2 = [(x,z,y): x&z->y], restricted to vocab terms."""
    nf1, nf2 = [], []
    with open(path) as f:
        for line in f:
            s = line.strip().replace('_', ':')
            if ' SubClassOf ' not in s:
                continue
            left, right = s.split(' SubClassOf ')
            if 'some' in left or 'some' in right:
                continue  # nf3/nf4 (existential roles) -- not per-class observable here
            if ' and ' in left:                       # A and B SubClassOf E
                a, b = left.split(' and ')
                if a in terms_dict and b in terms_dict and right in terms_dict:
                    nf2.append((terms_dict[a], terms_dict[b], terms_dict[right]))
            elif len(left) == 10 and len(right) == 10:  # A SubClassOf B
                if left in terms_dict and right in terms_dict:
                    nf1.append((terms_dict[left], terms_dict[right]))
    return nf1, nf2


def maxprop(P, nf1, n_terms, iters=64):
    """Deterministic ancestor-max closure over nf1 (the universal heuristic)."""
    q = P.clone()
    src = th.tensor([x for x, y in nf1]); dst = th.tensor([y for x, y in nf1])
    for _ in range(iters):
        newq = q.clone()
        # parent column dst gets max(own, children's current scores)
        newq.scatter_reduce_(1, dst.expand(q.shape[0], -1), q[:, src], 'amax')
        if th.allclose(newq, q, atol=1e-7):
            q = newq; break
        q = newq
    return q


def minprop(P, nf1, n_terms, iters=64):
    """Deterministic precision-direction dual of maxprop: enforce child <= parent by
    LOWERING children (true-path consistency by suppression instead of inflation)."""
    q = P.clone()
    src = th.tensor([x for x, y in nf1]); dst = th.tensor([y for x, y in nf1])
    for _ in range(iters):
        newq = q.clone()
        # child column src gets min(own, each parent's current score)
        newq.scatter_reduce_(1, src.expand(q.shape[0], -1), q[:, dst], 'amin')
        if th.allclose(newq, q, atol=1e-7):
            q = newq; break
        q = newq
    return q


def wmc_soft(P, nf1, nf2, n_terms, iters=128, downward=True):
    """Mean-field WMC marginals under the EL theory, conditioning the independent
    per-term priors on the true-path constraint.

    The exact marginal for a single implication A->B is b/(1-a(1-b)) -- a *bounded,
    normalised* soft-OR that RAISES the parent yet keeps it finite (unlike naive
    noisy-OR, which over-inflates and collapses precision). Generalised to a node Y
    with bodies {X} (nf1 children X->Y, nf2 conjunctions X&Z->Y):

        s_Y = 1 - prod_bodies (1 - body_prob)            # soft-OR support from below
        q_Y = (p_Y + (1-p_Y) s_Y) / Z   with the normaliser folded in as
        q_Y = p_Y / ( p_Y + (1-p_Y) * prod_bodies(1-body_prob) )   when p_Y is the
              independent prior and consistency forbids (body=1, Y=0).

    Downward pass (precision): a child must not exceed any parent -> cap q_X by min
    parent. Upward raises recall, downward enforces the dual; together they are the
    coherent joint the independence assumption breaks.
    """
    e = P.clamp(LO, 1 - LO)
    q = e.clone()
    n1s = th.tensor([x for x, y in nf1]); n1d = th.tensor([y for x, y in nf1])
    if nf2:
        n2a = th.tensor([a for a, b, y in nf2]); n2b = th.tensor([b for a, b, y in nf2])
        n2y = th.tensor([y for a, b, y in nf2])
    nrow = q.shape[0]
    for _ in range(iters):
        # ---- upward: normalised bounded soft-OR into each parent ----
        logprod = th.zeros_like(q)                       # sum log(1-body) into col Y
        logprod.scatter_add_(1, n1d.expand(nrow, -1), th.log1p(-q[:, n1s].clamp(LO, 1 - LO)))
        if nf2:
            joint = (q[:, n2a] * q[:, n2b]).clamp(LO, 1 - LO)
            logprod.scatter_add_(1, n2y.expand(nrow, -1), th.log1p(-joint))
        prod = th.exp(logprod)                           # prod_bodies (1 - body_prob)
        newq = e / (e + (1 - e) * prod).clamp_min(LO)    # normalised marginal
        # ---- downward: child <= min parent (true-path precision) ----
        if downward:
            newq.scatter_reduce_(1, n1s.expand(nrow, -1), newq[:, n1d], 'amin')
        newq = newq.clamp(LO, 1 - LO)
        if th.allclose(newq, q, atol=1e-6):
            q = newq; break
        q = newq
    return q


def gold_matrix(test_df, terms_dict, n_prot, n_terms):
    G = np.zeros((n_prot, n_terms), dtype=np.float32)
    for i, row in enumerate(test_df.itertuples()):
        for go_id in row.prop_annotations:
            j = terms_dict.get(go_id)
            if j is not None:
                G[i, j] = 1
    return G


def fmax(P, G):
    """protein-centric Fmax over a 0..1 threshold sweep (DeepGOZero convention)."""
    best = 0.0; bt = 0.0
    gold_cnt = G.sum(1)
    for t in range(1, 101):
        thr = t / 100.0
        pred = (P >= thr)
        pcnt = pred.sum(1)
        tp = (pred * G).sum(1)
        has = pcnt > 0
        if has.sum() == 0:
            continue
        prec = (tp[has] / pcnt[has]).mean()
        rec = (tp / np.maximum(gold_cnt, 1)).mean()
        if prec + rec > 0:
            f = 2 * prec * rec / (prec + rec)
            if f > best:
                best = f; bt = thr
    return best, bt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-root', default='data')
    ap.add_argument('--ont', default='cc')
    ap.add_argument('--model', default='mlp_raw', help='predictions_<model>.pkl provides raw preds')
    args = ap.parse_args()

    terms = pd.read_pickle(f'{args.data_root}/{args.ont}/terms.pkl')['gos'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    n_terms = len(terms)
    nf1, nf2 = parse_norm(f'{args.data_root}/go.norm', terms_dict)
    print(f'{args.ont}: {n_terms} terms | nf1(A->B)={len(nf1)} nf2(A&B->E)={len(nf2)}')

    df = pd.read_pickle(f'{args.data_root}/{args.ont}/predictions_{args.model}.pkl')
    P = th.tensor(np.stack(df['preds'].values), dtype=th.float32)
    n_prot = P.shape[0]
    G = gold_matrix(df, terms_dict, n_prot, n_terms)
    Gt = th.tensor(G)
    print(f'proteins={n_prot} gold pairs={int(G.sum())} prevalence={G.mean():.5f}')

    variants = {
        'raw':            P,
        'maxprop':        maxprop(P, nf1, n_terms),
        'minprop':        minprop(P, nf1, n_terms),
        'wmc-up(nf1)':    wmc_soft(P, nf1, [], n_terms, downward=False),
        'wmc(nf1)':       wmc_soft(P, nf1, [], n_terms),
        'wmc(nf1+nf2)':   wmc_soft(P, nf1, nf2, n_terms),
    }
    print(f'{"variant":16} {"Fmax":>7} {"thr":>5} {"AUPR":>7}')
    for name, M in variants.items():
        Mn = M.numpy()
        f, t = fmax(Mn, G)
        aupr = average_precision_score(G.ravel(), Mn.ravel())
        print(f'{name:16} {f:7.3f} {t:5.2f} {aupr:7.3f}')


if __name__ == '__main__':
    main()
