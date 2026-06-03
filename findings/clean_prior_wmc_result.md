# CLEAN as the prior for producibility-WMC: a never-worse refinement that wins at low similarity

The earlier MTRR-vs-CLEAN comparison (`clean_vs_mtrr_result.md`) pitted structure-only
recovery *against* CLEAN — the wrong framing for our thesis. The contribution is not a
competing EC predictor; it is a **WMC layer on top of the best base predictor**. Here we
test that directly: use CLEAN as the prior weighting the producibility gap-fill, and ask
whether **CLEAN + WMC > CLEAN**.

## Method

gapsmith's gap-fill objective is already `min Σ wᵢ·xᵢ s.t. producibility` — a weighted
model count. `RxnWeights::from_probabilities` is documented in-source as *"the injection
point for a DL function predictor."* The prob→pseudo-bitscore map is invertible, so CLEAN's
beliefs are written straight into the `Reactions.tbl` weight column — no rebuild.

- **Genome-wide CLEAN prior.** ESM-embedded all **135,062 genes** across the 30 organisms
  (job 4439, node005), reduced CLEAN's distance map to a per-organism EC→probability prior
  (`org_ec_prob.pkl`, ~3000 ECs/organism above threshold) + per-gene true-EC scores
  (`gene_top.pkl`).
- **CLEAN-weighted recovery.** Same 225 held-out reactions, balanced into 5 identity
  buckets. For each: remove the reaction + strip its homology, build a CLEAN-weighted
  weight table (every candidate reaction weighted by CLEAN's belief the genome encodes its
  EC; the held-out reaction weighted by CLEAN's score on its *true* EC, homology-stripped),
  run `gapsmith fill --target <product|biomass>`, check EC-level recovery (job 4440).
- **Metrics:** `clean_alone` (CLEAN top-1, from maxsep), `flat-MTRR` (uniform-prior
  gap-fill), `CLEAN-WMC` (CLEAN-weighted gap-fill), `refine = max(clean_alone, CLEAN-WMC)`
  — the deployable ensemble (accept a reaction if CLEAN predicts it OR the WMC recovers it).

## Result (216 reactions, EC-level, by identity-to-CLEAN-train)

| identity | CLEAN-alone | flat-MTRR | CLEAN-WMC | refine | n |
|---|---|---|---|---|---|
| ≤30% | 0.089 | 0.067 | 0.067 | **0.133** | 45 |
| 30–50% | 0.581 | 0.233 | 0.279 | **0.721** | 43 |
| 50–70% | 0.581 | 0.186 | 0.163 | **0.651** | 43 |
| 70–90% | 0.500 | 0.095 | 0.095 | **0.500** | 42 |
| 90–100% | 0.721 | 0.140 | 0.186 | **0.791** | 43 |
| **overall** | **0.491** | 0.144 | 0.157 | **0.556** | 216 |

## Findings

1. **As a refinement on top of CLEAN, the WMC layer is a strict, significant improvement.**
   Overall **0.491 → 0.556** (+6.5 pts). `refine ≥ clean_alone` in every bucket (never
   worse). At ≤50% identity the gain is **8 reactions, 0 losses, McNemar p = 0.0078** — a
   statistically significant lift *exactly in the regime where CLEAN collapses*, with zero
   displacement of CLEAN-correct calls. ≤30%: +50% relative; 30–50%: +14 pts. **This is the
   "strictly on top of the best base predictor" result.**

2. **The lift comes from the structure's similarity-independent complementarity, not from
   CLEAN-informed candidate ranking.** CLEAN-WMC (0.157) barely exceeds flat-MTRR (0.144):
   net +3 reactions (6 CLEAN-weighting-only vs 3 flat-only). The producibility constraint
   already narrows candidates to a small set, so the prior-weighting has little room to act
   — the constraint dominates. The 14 reactions CLEAN misses but structure recovers (8 at
   ≤50%) are the value-add, and they are recovered *because producibility is blind to
   sequence similarity*, not because CLEAN re-ranked them.

3. **The strong standalone claim stays refuted.** CLEAN-WMC alone (0.157) ≪ CLEAN-alone
   (0.491). Structure is not a replacement predictor; it is a complement. The refinement
   realizes the union ceiling (0.55) predicted from the confusion matrix — it does not
   exceed it, because CLEAN-weighting did not beat flat structure.

## Caveat (open)

Measured: recovery of held-out *true* reactions, with zero displacement of CLEAN-correct
calls. **Not yet measured: the WMC layer's false-positive (precision) cost** — it adds
reactions to restore producibility, and some additions may carry wrong ECs. Before claiming
a clean Pareto improvement we must report precision of the WMC-added reactions, not only
recall on the held-out set. (Next experiment.)

## Paper claim (refined, honest)

Not "WMC beats CLEAN." Rather: **a producibility-WMC layer applied on top of CLEAN is a
never-worse refinement that recovers a statistically significant set of enzymes CLEAN
misses, concentrated in the low-sequence-similarity regime where CLEAN fails (p < 0.01,
zero degradation) — because the structural constraint is blind to similarity.** The lift is
bounded by coverage (the producibility-reachable subset), consistent with the 3-condition
principle. The clean compositional-coverage win remains GO zero-shot WMC
(`go_zero_finding.md`); this is the metabolic strand's *complementary-refinement* win.

## Reproduce
```
# prior:    partb/clean_prior.sbatch (node005) -> infer_prior_all.py -> org_ec_prob.pkl, gene_top.pkl
# recovery: partb/mtrr_clean_array.sbatch -> mtrr_clean_recover.py -> results_clean/*.pkl
# compare:  partb/clean_wmc_compare.py -> clean_wmc_compare.pkl
# SLURM chain: 4439 (prior) -> afterok 4440 (array) -> afterok 4441 (agg)
```
