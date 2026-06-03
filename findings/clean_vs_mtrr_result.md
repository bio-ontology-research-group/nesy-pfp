# MTRR vs CLEAN at scale: complementary, not dominant (honest negative for the strong claim)

The single-case secondary-metabolism recovery (`clean_partB_secondary_recovery.md`) is
here generalized to a **method** — *Metabolite-Targeted Reaction Recovery* (MTRR) — and
evaluated at scale against CLEAN on held-out reactions, stratified by sequence identity
to CLEAN's own training set.

## Method (MTRR)

For each evaluable reaction in a gap-filled gapsmith model: delete its stoichiometry
column **and** strip its sequence evidence from `Reactions.tbl` (simulating an
uncharacterised enzyme), then run `gapsmith fill --target <metabolite>` and ask whether
the reaction (or, at EC level, *any* reaction carrying one of its true EC numbers) is
re-added from metabolic necessity alone — blind to sequence. Two targets per reaction:
the reaction's own non-cofactor **product** (`prod`) and **biomass** (`cpd11416`).

Scale: 225 sampled reactions across the 30 organism models, balanced into 5 buckets of
sequence identity to CLEAN's `split100` training set (45 reactions/bucket). EC-level
functional recovery (the fair metric: credit if the recovered route carries a correct EC).
SLURM array job 4422 on unimatrix01 (node002/003/004).

## MTRR recovery (EC-level), by identity-to-CLEAN-training

| target | ≤30 | 30–50 | 50–70 | 70–90 | 90–100 | overall |
|---|---|---|---|---|---|---|
| `prod` (product-targeted) | 0.067 | 0.222 | 0.178 | 0.089 | 0.133 | **0.138** |
| `biomass` | 0.022 | 0.000 | 0.022 | 0.067 | 0.022 | 0.029 |

The product target carries essentially all recovery; biomass adds nothing (confirms the
Part B negative — growth-coupled producibility is blind to peripheral enzymes). MTRR is
**flat across similarity** (no decay with dissimilarity) but **weak** (~14%).

## Head-to-head vs CLEAN (216 shared reactions, EC-level, per reaction = any catalysing gene)

CLEAN inference on the 485 catalysing genes (node005 GPU, job 4438).

| identity-to-CLEAN-train | CLEAN | MTRR | n |
|---|---|---|---|
| ≤30% | **0.089** | 0.067 | 45 |
| 30–50% | 0.581 | 0.233 | 43 |
| 50–70% | 0.581 | 0.186 | 43 |
| 70–90% | 0.500 | 0.095 | 42 |
| 90–100% | 0.721 | 0.140 | 43 |
| **overall** | **0.491** | 0.144 | 216 |

Confusion (216 reactions): both 18 · CLEAN-only 88 · MTRR-only 13 · neither 97.

## Findings (honest)

1. **CLEAN collapses with dissimilarity — reconfirmed.** 0.72 → 0.089 from the
   highest to the lowest identity bucket (8× drop), an independent replication of the
   Part A memorisation signature (`clean_temporal_finding.md`) on a fresh reaction set.

2. **MTRR does NOT rescue the failure regime — the strong claim is REFUTED.** At ≤30%
   identity MTRR recovers 3/45 (0.067), statistically tied with CLEAN's floor (4/45).
   McNemar on the ≤50% region is **p = 0.007 *in CLEAN's favour*** (CLEAN-only 24 vs
   MTRR-only 8). Structure-only gap-fill recovery is too sparse to beat CLEAN even where
   CLEAN is weakest. Goal success-criterion (ii) ("ESM+WMC ≥ CLEAN in ≤50% buckets") is
   **not met** by the structure-only limb.

3. **MTRR is genuinely complementary, not a replacement.** It recovers **13 reactions
   CLEAN missed** (8 at ≤50% identity, 2 at ≤30%), lifting union coverage 0.491 → **0.551**.
   These are similarity-independent discoveries the homology/contrastive model cannot make,
   sourced purely from metabolic structure.

## Why (3-condition principle, sharpened again)

Gap-fill/producibility structure fails **condition 3 (coverage)**: it only reaches enzymes
whose product is a non-cofactor metabolite on a *unique* ModelSEED route inside the
producible scope. Most novel enzymes are peripheral/specialised and fall outside that
scope, so MTRR can recover only a flat ~14% minority. This is the same coverage wall that
falsified the BacDive route. Contrast with GO zero-shot WMC, which satisfies all three
conditions (logical definitions compositionally cover the label space) and *does* beat its
baseline — that remains the strand's clean positive, alongside Part A.

## Defensible paper claim

Not "WMC beats CLEAN." Rather: **(a)** CLEAN's accuracy is similarity-driven and collapses
on dissimilar enzymes (memorisation); **(b)** a sequence-blind, metabolic-structure
recovery method is *similarity-independent* and *complementary*, recovering a real set of
reactions CLEAN misses (union lift 0.49 → 0.55), but **(c)** structure-only recovery is
sparse and does not dominate — coverage of the target space is the binding constraint,
exactly as the 3-condition principle predicts.

## Reproduce
```
# EC-level MTRR (cluster): partb/mtrr_array.sbatch -> results/mtrr_res_*.pkl
ssh unimatrix01 'python ~/nesy-genome/partb/agg_ec.py'           # -> mtrr_ec_byrxn.pkl
# CLEAN: partb/clean_mtrr.sbatch on node005 -> mtrr_clean_maxsep.csv
ssh unimatrix01 'python ~/nesy-genome/partb/clean_vs_mtrr.py'    # -> clean_vs_mtrr.pkl + tables
```
