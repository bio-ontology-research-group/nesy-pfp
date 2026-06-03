# When the soft true-path fixpoint equals exact WMC — theory (from the Lean development)

*Every claim below is anchored by a kernel-checked theorem in
`metabolic-pfp/theory/lean-wmc-approx/WMCApprox.lean` (Lean 4 core, no `sorry`; axioms ⊆
{`propext`, `Quot.sound`}). Theorem names are given in `monospace`.*

## Setup

GO prediction emits independent per-term marginals `p_i = σ(logit_i)`. The ontology imposes
the **true-path constraint** — for every subsumption edge `child ⊑ parent`, `child ⟹ parent`.
The principled object is the **weighted model count (WMC)**: condition the independent product
distribution on the conjunction `⋀_edges (¬child ∨ parent)` and read per-term marginals
(`wmcMarg`). The universal heuristic (`max`-propagation, DeepGO* default) and our soft
noisy-OR closure both *approximate* this. The question is exactly *when the approximation is
the WMC*.

## 1. The approximation is belief propagation; it is exact on polytrees, wrong at reconvergences

The soft closure is belief propagation on the Horn constraint graph: it combines incoming
messages **as if independent**. For a node with two parents this is the independence form

    P(sink) = p_d·p_b·p_c / ( p_d·p_b·p_c + (1−p_d) )      (`bpSink`).

- **Independent parents (polytree) → exact.** When the two parents are independent roots,
  the WMC marginal equals this form exactly: at `p=½`, both are `1/5`
  (`polytree_sink_exact`).

- **Shared-ancestor parents (diamond) → wrong.** Add a common ancestor `a` of both parents
  (the reconvergence that the structure probe shows is pervasive in GO). The *same* form now
  disagrees with the WMC: exact `1/6` vs BP `1/5` at the sink (`diamond_sink_bp_wrong`), and
  at the reconvergent node itself exact **5/6** vs BP **7/10**
  (`diamond_reconv_inexact`, with the values pinned by `diamond_exact_a_is_5_6` and
  `diamond_bp_a_is_7_10`).

The mechanism is precisely the shared-ancestor correlation: BP equals exact inference on the
*computation-tree unrolling* that duplicates the shared ancestor, so it is exact iff the
unrolling is the identity — iff the constraint graph is a **polytree**. A two-parent node is
not itself the problem; a *reconvergence* is.

## 2. The closure always converges

The synchronous true-path closure `step` (a parent is activated by any active child) is

- **extensive** — it only ever adds (`step_ext`),
- **monotone** in the active set (`step_mono`), and
- **stable at its fixpoint** — once a step adds nothing, every further iteration is fixed
  (`fixed_stable`, which needs *no axioms at all*).

Because the active set grows monotonically inside a finite vocabulary, iteration reaches the
least fixpoint in at most `depth` steps. So convergence is unconditional; only *what it
converges to* is structure-dependent (§1).

## 3. Why this is the design spec for modularization (ties to the structure probe)

The empirical probe (`go_structure_probe.md`) measured exactly the quantity §1 says governs
the error — the reconvergence / treewidth of each namespace:

| namespace | circuit rank | treewidth UB | regime |
|---|---|---|---|
| cc | 654 | **11** | near-exact: soft fixpoint ≈ WMC; modules buy little |
| mf | 1 938 | **44** | borderline: modular exact WMC feasible and worth demonstrating |
| bp | 24 432 | **1 810** | massively reconvergent: monolithic exact WMC impossible; modules essential, residual loops bounded by §1 |

So the theory and the data agree: the approximation is exact on the tree fragment, and its
error lives **exactly at the reconvergences**, which are rare in `cc` (tw 11) and pervasive in
`bp` (tw 1810). This is the cut criterion for locality modules: place module boundaries around
the reconvergent cores, compile each module to an exact SDD, and let the soundly-coordinated
inter-module messages (the `kobayashi-marust` context calculus, lifted from Boolean
consequences to weighted separator marginals) carry the residual — exact where the
module-interface graph is a tree, §1-bounded approximate where it is not.

## Honest scope

What is machine-checked: convergence (general), and exactness-on-polytrees vs
inexactness-at-reconvergence (concrete witnesses, `decide`). What is **not yet** proved in
full generality: the *all-priors* tree-exactness theorem (the witnesses instantiate it; the
parametric proof is the next Lean increment, and is the EL analogue of moose's
`complete_EL`). The quantitative error bound on general DAGs is `#P`-hard in closed form and
is reported empirically against the modular-exact SDD on `cc`/`mf`.
