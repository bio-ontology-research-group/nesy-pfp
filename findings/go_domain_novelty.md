# What the GO predictors actually memorize: InterPro-domain novelty is the failure axis

Prompted by skepticism of ProteInfer's apparent "robustness": its authors prove (class-
activation mapping) that it keys on domain-like motif regions, and it underperforms BLAST on
the UniRef50 remote-homology split. So the robustness I first reported was measured on the
WRONG axis — *global sequence identity* — whereas a motif learner should be stratified by
*domain novelty*. Here we bucket test proteins by **InterPro-domain overlap with training**
(`go_domain.py`, job 4508): `no-IPR` (no domains), `all-novel` (no domain seen in training),
`partial`, `all-seen`.

## Result — Fmax by InterPro-domain novelty (mf; cc/bp identical shape)

| predictor | no-IPR | all-novel | partial | all-seen | ALL |
|---|---|---|---|---|---|
| Naive (frequency floor) | 0.306 | 0.368 | 0.363 | 0.305 | 0.326 |
| **BLAST** (homology) | 0.191 | **0.025** | 0.590 | 0.692 | 0.623 |
| **MLP** (InterPro input) | 0.309 | **0.369** | 0.659 | 0.694 | 0.657 |
| **DeepGOZero** (InterPro input) | 0.306 | **0.367** | 0.649 | 0.698 | 0.657 |
| DeepGOCNN (seq CNN) | 0.314 | 0.404 | 0.491 | 0.392 | 0.430 |
| **ProteInfer** (seq CNN) | **0.523** | **0.679** | 0.696 | 0.702 | 0.697 |

(bp: BLAST 0.49→0.042 on all-novel; ProteInfer 0.63→0.60→0.51 stays high. cc: same pattern.)

## Findings

1. **BLAST collapses to ≈0 on domain-novel proteins (0.025).** Homology transfer is dead when
   no training protein shares a domain — the clean memorization signature, now on the correct
   (domain) axis rather than global identity.
2. **MLP and DeepGOZero — the lab's own InterPro-feature models — collapse to the Naive
   frequency floor on novel/absent domains** (all-novel 0.37 = Naive 0.37; no-IPR 0.31 = Naive
   0.31). Their *input is the InterPro vector*; a novel domain is an uninformative input, so
   they retain no signal beyond label frequency. **They are InterPro-domain memorizers.**
3. **ProteInfer is the exception and generalises beyond InterPro.** It holds 0.68 on
   all-novel-domain proteins and **0.52 on no-IPR proteins — the only method above the floor
   there.** So it learns **sub-InterPro motifs** that transfer across InterPro-domain
   boundaries; InterPro-novelty does not break it (its mild 0.70→0.52 drop on no-IPR is the
   closest visible failure, pointing at genuinely novel sub-motifs/folds as its true limit).
4. **DeepGOCNN is weak and near-Naive throughout** — a poor base.

## Why this matters for the NMI thesis

The **failure regime is domain novelty**, and it bites exactly the **feature-based predictors
(DeepGOZero, MLP)** that the GO-WMC layer sits on top of: on domain-novel inputs they fall to
the frequency floor (0.37), losing all functional signal. That is the within-vocab analogue of
the zero-shot regime — a concrete, large population (mf: 281 all-novel + 78 no-IPR proteins)
where the base predictor fails and ontology structure is the only remaining handle. **Next
experiment:** does `base + WMC` lift DeepGOZero/MLP on domain-novel proteins (stratified Fmax)?
If the hierarchy/definition closure recovers signal where the InterPro features are blank, that
is a within-vocab failure-regime win — stronger than the (neutral) all-proteins within-vocab
result and complementary to the zero-shot-term result.

Honest caveat: ProteInfer, being a finer motif learner, does *not* collapse on this axis, so a
WMC lift would apply to the feature-based bases, not to the strongest sequence base. Script:
`go_domain.py`.

## Follow-up: does WMC lift the feature-based bases on domain-novel proteins? NO (job 4518)

| variant | no-IPR | all-novel | partial | all-seen |
|---|---|---|---|---|
| DeepGOZero raw / +WMC | 0.306 / 0.309 | 0.367 / **0.365** | 0.649 / 0.651 | 0.698 / 0.697 |
| MLP raw / +WMC | 0.309 / 0.309 | 0.369 / **0.369** | 0.659 / 0.658 | 0.694 / 0.692 |

**WMC is neutral on domain-novel proteins — it does not recover the collapsed base.** The
mechanism is decisive and sharpens the whole thesis: WMC composes/propagates the base's *own*
per-term marginals. On a domain-novel protein the feature-based base has collapsed to the
frequency floor for **every** term, so there is no protein-specific signal on any component for
the closure to compose. **WMC amplifies existing perception signal across the ontology; it
cannot manufacture signal where perception has failed.**

### The precise boundary (a key NMI delimitation)
- **Zero-shot TERM** — base scores the protein's *component* terms but not the held-out target:
  a *compositional* gap → **WMC works** (training-free, matches trained zero-shot, AUC 0.91).
- **Domain-novel PROTEIN** — base scores *no* term for this protein (perception failure):
  no informative parts → **WMC fails**; the fix is better perception (e.g. ProteInfer's
  sub-motifs), not structure.

So the WMC layer's GO value lives strictly in the zero-shot-*term* regime (external,
compositional structure over a base that still perceives the parts), NOT the domain-novel-
*protein* regime. This is the cleanest statement of the project's core principle: WMC helps iff
the structure is external to supervision AND the base retains signal on the structure's
components.
