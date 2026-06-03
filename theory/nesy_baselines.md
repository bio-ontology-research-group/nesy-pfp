# DeepProbLog vs Scallop on GO true-path WMC: the exact-wall and the approx-drift

The two canonical NeSy frameworks, run on the same GO true-path inference, to position the
modular/soft-WMC approach. Both use the actual installed tools (ProbLog 2.2.10 = DeepProbLog's
inference engine; Scallop 0.2.4 `topkproofs`). Scripts: `dpl_scaling.py`, `scallop_scaling.py`.

## Correctness — both reproduce our exact WMC

- **ProbLog** (true-path constraint as evidence, `evidence(violated,false)`): on the diamond,
  `P(go(a)) = 0.8333 = 5/6`, `P(go(d)) = 1/6` — exactly our Lean theorem and junction-tree value.
- **Scallop** (upward derivation `holds(P) :- sub(C,P), holds(C)`): the chain a→x,b→x at p=½
  gives `P(holds(x)) = 0.875 = 1−½³`, exact.

## DeepProbLog / ProbLog — EXACT, but hits a hard scalability wall

True-path WMC on GO subgraphs of increasing size (knowledge compilation → exact WMC):

| ontology | 100 terms | 250 terms | 400 terms |
|---|---|---|---|
| mf | OK 0.3 s | OK 2.0 s | **TIMEOUT (>300 s)** |
| bp | OK 0.3 s | OK 0.8 s | **TIMEOUT (>300 s)** |

A sharp crossover at ~250→400 terms. Full GO is **2,800 (cc) / 6,900 (mf) / 21,000 (bp)**
terms — one to two orders of magnitude beyond ProbLog's practical limit. The exact-NeSy wall is
empirical here and *rigorous* via the structure probe (bp treewidth 272 ⇒ `2²⁷²` circuit).

## Scallop — SCALES, but the top-k approximation drifts at reconvergences

`topkproofs` runs full GO trivially: **cc 1,000 terms in 0.1 s; mf 6,796 in 0.2 s; bp 21,324
(all) in 0.9 s** — 2–3 orders of magnitude past ProbLog. But top-k truncates proofs, and the
error lives **exactly at reconvergent (multi-proof) nodes** — the regime our Lean theorem
(`diamond_reconv_inexact`) characterises:

- **Diamond** `P(holds(a))`, exact 0.9375: Scallop **k=1 → 0.500** (−47%), **k=3 → 0.875**,
  **k=10 → 0.9375** (exact).
- **14-node bp balls (circuit rank 7–12):** max drift vs exact — **k=1: 0.15–0.30**,
  **k=3: ≤0.013**, **k=10: 0.000**.

Drift grows with proof multiplicity (reconvergence) and shrinks with k. Crucially, at full-GO
scale (bp circuit rank **24,432**) a hub term has astronomically many proofs, so any fixed k
systematically undercounts there — **and this error is unverifiable**, because no exact
reference is computable at that scale (that is the whole #P-hardness point). Scallop scales by
giving up the guarantee.

## Where our approach sits (the differentiator)

| method | scales to full GO? | accuracy guarantee |
|---|---|---|
| DeepProbLog (exact KC) | **no** (wall ~400 terms) | exact (when it runs) |
| Scallop (top-k) | yes (<1 s) | **none at scale** (top-k drift, unverifiable) |
| **ours (soft-WMC + modular)** | yes | **exact on the tree fragment (Lean-proven), modular-exact over 92–100 % of GO, soft-fixpoint validated vs exact JT on cc (mean err 6×10⁻⁴)** |

So the comparison is honest in both directions: exact NeSy (DeepProbLog) is correct but cannot
reach GO scale; approximate NeSy (Scallop) reaches scale but trades away the accuracy guarantee,
with error concentrated exactly where our theory predicts (reconvergences). The contribution is
the method that is **both** scalable **and** accuracy-characterised — machine-checked exact on
the tree fragment, modular-exact on the low-treewidth bulk, and bounded-approximate only on the
small reconvergent core, rather than top-k-truncated everywhere with no guarantee.

## Learned approximate inference — A-NeSI and NeSyDM (a different category)

DeepProbLog and Scallop are *inference engines*: run on a fixed WMC, so the GO scaling sweep is
apples-to-apples. **A-NeSI** (van Krieken, NeurIPS 2023 — amortized inference: a neural surrogate
*trained* to predict the WMC output) and **NeSyDM** (van Krieken 2025 — masked discrete diffusion
over the symbolic variables, modelling the joint to break independence) are *end-to-end training
methods*, not inference engines. We do not re-run them on GO for two reasons:

1. **They are not the native comparison.** Putting them on GO means *training* a surrogate /
   diffusion model per ontology. A-NeSI's amortised inference must be *trained on samples of the
   exact WMC* — which on GO are intractable to generate at scale (the very #P-hardness we attack),
   so its amortisation is bootstrapped from an oracle one does not have. Our layer computes the
   exact marginal *directly* on the tree fragment, with no inference-training data.
2. **We already have their empirical behaviour on a structured task** (the WMS-NeSy carry-chain /
   MNIST-addition, `archive/wms-nesy/results/`, 3 seeds), and it is the relevant signal — both
   *scalable approximate* methods **collapse at moderate structure depth**:

   | method | N=1 (1 digit) | N=4 (4-digit carry chain) |
   |---|---|---|
   | A-NeSI (amortised) | 0.968 | **0.225** |
   | NeSyDM (diffusion) | 0.962 | **0.03–0.13** |

   At N=4 the carry constraint couples four digits; both learned approximators break down
   (A-NeSI 0.97→0.23; NeSyDM 0.96→~0.03), exactly the independence-assumption / reasoning-shortcut
   failure this project targets. So the learned-approximate route does not escape the wall — it
   relocates it from *exact-inference cost* (DeepProbLog) to *approximation quality on coupled
   structure* (Scallop's reconvergence drift; A-NeSI/NeSyDM's depth collapse), in every case
   without an exactness guarantee.

**Positioning.** Across the NeSy spectrum — exact knowledge compilation (DeepProbLog), top-k
provenance (Scallop), amortised inference (A-NeSI), joint diffusion (NeSyDM) — each trades
correctness for scale and offers no guarantee on the coupled-structure regime. Our contribution
is the one that is exact and **machine-checked** on the tree fragment, modular-exact over
92–100% of GO, and bounded-approximate only on the small reconvergent core.

Reproduce: `dpl.sbatch` (job 4528), `scallop.sbatch` (job 4525); venvs `nesy-venv` (ProbLog),
`scallop-venv` (Scallop 0.2.4, cp310 wheel from GitHub release 0.2.4). A-NeSI/NeSyDM carry-chain
logs: `archive/wms-nesy/results/{anesi,nesydm}_*seed*.log`; NeSyDM impl `moose/moose/baselines/nesydm.py`.
