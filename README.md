# Verified weighted-model-counting for protein function prediction

> **Protein function predictors memorize; a machine-checked weighted-model-counting (WMC)
> layer generalizes where they fail.**

This repository documents the experiments, analyses, and machine-checked proofs behind our
study of the *conditional-independence assumption* in neuro-symbolic protein function
prediction (PFP). It accompanies the manuscript of the same name.

The work has two halves, joined by one principle:

1. **Critique.** Under similarity- and domain-controlled evaluation against the canonical
   CAFA *Naive* and *BLAST* baselines, leading PFP predictors (CLEAN for EC numbers, the
   DeepGO family for GO) succeed by **memorizing homology and curated domain features** and
   **collapse out-of-distribution**.
2. **Remedy.** Compiling the *symbolic structure* of the target space — ontology axioms for
   GO, metabolic producibility for EC — into a differentiable **weighted-model-counting layer
   on top of any base predictor** restores generalization, but only in a precisely delimited
   regime, and it is **provably exact on the tree fragment** and **modular-exact over
   92–100 % of GO**.

**Governing principle.** *WMC helps iff the structure is external to supervision **and** the
base predictor retains signal on the structure's components.*

---

## Headline results

| # | Finding | Evidence |
|---|---------|----------|
| 1 | CLEAN's EC accuracy collapses **0.86 → 0.07** with sequence identity to its own training set, on provably-unseen post-publication enzymes | `findings/clean_temporal_finding.md` |
| 2 | Against CAFA baselines, **BLAST → 0.025** and **MLP/DeepGOZero → the Naive floor (0.37)** on domain-novel proteins; **ProteInfer** is the exception (0.68), learning sub-InterPro motifs | `findings/go_cafa_baselines.md`, `findings/go_domain_novelty.md` |
| 3 | **Training-free** ontology-WMC gives any base predictor zero-shot reach that **matches a *trained* zero-shot extrapolator** (mf macro-AUC **0.9125 vs 0.9118**, p=0.83), while Naive/BLAST are ≈0.5 | `findings/go_nmi_benchmark.md`, `experiments/scale_zero_h2h.py` |
| 4 | The layer helps the **compositional** failure regime (zero-shot terms) but **not** the **perception** failure regime (domain-novel proteins) | `findings/go_domain_novelty.md` |
| 5 | EC producibility-WMC is a similarity-independent **candidate layer**: CLEAN **0.491 → 0.556** recall, but ~72-reaction blast radius at 0.69 precision (intrinsic, growth-coupled) | `findings/clean_prior_wmc_result.md`, `findings/clean_wmc_precision.md` |
| 6 | The soft WMC closure is **belief propagation**: exact on trees, error confined to reconvergences — **machine-checked in Lean 4** (general theorem for all priors) | `theory/approximation_theory.md`, `theory/lean/` |
| 7 | Constraint-graph treewidth **cc 9 / mf 31 / bp 272**; cutting **<8 % reconvergence hubs** makes **92–100 % of GO exactly computable** | `theory/modular_wmc_result.md`, `theory/marust_coordination.md` |
| 8 | NeSy frameworks: **DeepProbLog** (exact) walls at ~400 terms; **Scallop** (top-k) scales but drifts at reconvergences; **A-NeSI/NeSyDM** collapse on coupled structure (N=4) | `theory/nesy_baselines.md` |

---

## Repository structure

```
experiments/   evaluation + analysis scripts (GO and EC)
findings/       per-experiment write-ups with numbers and honest scope
theory/         approximation theory, modular WMC, NeSy baselines
theory/lean/    machine-checked Lean 4 developments (no `sorry`)
```

### Experiments (`experiments/`)
| script | what it does |
|---|---|
| `go_wmc_eval.py` | within-vocab GO: raw vs max-propagation vs soft-WMC closure (Fmax/AUPR) |
| `scale_zero_h2h.py` | zero-shot GO: training-free WMC vs trained DeepGOZero-zero (macro-AUC) |
| `go_zero_scale.py` | masked-definition reconstruction of compositionally-defined terms |
| `go_mem.py` | Fmax stratified by sequence identity to training (memorization) |
| `go_cafa.py` | CAFA Naive + BLAST + learned predictors, Fmax by similarity |
| `go_domain.py` | Fmax stratified by **InterPro-domain novelty** (the true failure axis) |
| `structure_probe.py` | per-namespace circuit rank, multi-parent count, treewidth bound |
| `jt_wmc.py` | junction-tree exact WMC + loopy BP (validated on the Lean diamond) |
| `cc_modular.py` | exact JT-WMC vs soft fixpoint on cc; error localized at reconvergences |
| `tract_probe.py`, `modulator_probe.py` | treewidth ladder + irreducible-core (<8 %) |
| `dpl_scaling.py` | DeepProbLog/ProbLog exact-WMC scaling sweep (the wall) |
| `scallop_scaling.py` | Scallop top-k scaling + drift-vs-reconvergence |

### Machine-checked theory (`theory/lean/`)
- `lean-wmc-approx/WMCApprox.lean` — Lean 4 **core**, no `sorry`, axioms ⊆ {`propext`,
  `Quot.sound`}. Exactness witnesses (polytree exact; diamond wrong; 5/6 vs 7/10) and
  unconditional convergence of the true-path closure.
- `lean-wmc-approx-general/WMCBattery.lean` — breadth battery: single-edge BP = exact WMC
  across an 81-prior grid; polytree exact / diamond wrong across priors (all `decide`).
- `wmc_general/WmcGeneral/Star.lean` — **general theorem** (mathlib): for a parent with `k`
  independent children under the true-path constraint, the exact WMC marginal equals the
  soft-OR update *for all k and all priors* (`star_marginal_eq_softOR`), no `sorry`.

Build the core developments with a Lean 4 toolchain (`leanprover/lean4:v4.30.0`):
```bash
cd theory/lean/lean-wmc-approx       && lean WMCApprox.lean
cd theory/lean/wmc_general           && lake exe cache get && lake build   # needs mathlib
```

---

## Reproduction notes

- **Data.** GO experiments use the DeepGOZero SwissProt release (terms, train/test splits,
  per-model prediction pickles, InterPro annotations, `go.norm` EL axioms). EC experiments use
  post-2023 reviewed UniProt enzymes and organism genome-scale models built with `gapsmith`.
- **NeSy baselines.** ProbLog 2.2.10 (DeepProbLog's inference engine) and Scallop 0.2.4
  (`topkproofs`, cp310 wheel from the GitHub release) both reproduce our exact WMC on the
  validation cases (diamond `P(a)=5/6`; chain `0.875`). A-NeSI / NeSyDM carry-chain numbers are
  from the WMS-NeSy strand.
- Scripts are written for a SLURM cluster; each finding doc lists the exact job and inputs.

## Honest scope

Every claim is bounded, and the bounds are part of the contribution: within-vocabulary GO WMC
is **neutral**; the GO zero-shot result is a **match** to a trained extrapolator, not a beat;
WMC does **not** repair perception failure (domain-novel proteins); the EC layer is
**candidate-generation** at a stated precision; exact monolithic WMC is impossible on bp
(treewidth 272) and is recovered only **modularly**, with a small bounded-approximate core.

## License

Code and documentation released under the MIT License (see `LICENSE`).

## Manuscript

The full manuscript source (LaTeX + bibliography) is in `paper/` and is maintained on Overleaf. It consolidates every experiment and proof documented here: the memorization critique, the training-free WMC layer and its competence boundary, the machine-checked approximation theory, the modular construction, and the full neuro-symbolic-framework comparison.
