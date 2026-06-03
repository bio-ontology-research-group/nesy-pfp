# CAFA baselines (Naive + BLAST): what the GO predictors actually generalize

A reviewer's check (and the right one): without CAFA's **Naive** (predict each term at its
training frequency) and **BLAST** (diamond annotation transfer) baselines, a predictor that
merely memorizes label frequency + the hierarchy looks good. We added both, Fmax stratified by
sequence identity to training (`go_cafa.py`, job 4496).

## Result — mf (Fmax by max-identity-to-train)

| predictor | <30 | 30–50 | 50–70 | ALL |
|---|---|---|---|---|
| **Naive (frequency)** | 0.343 | 0.318 | 0.406 | 0.326 |
| **BLAST (homology)** | 0.472 | 0.706 | 0.752 | 0.623 |
| MLP | 0.597 | 0.697 | 0.733 | 0.657 |
| DeepGOCNN | 0.431 | 0.427 | 0.536 | 0.430 |
| **ProteInfer** | **0.704** | 0.687 | 0.747 | 0.697 |
| DeepGOZero | 0.593 | 0.698 | 0.741 | 0.657 |

(cc shows the same shape: Naive 0.60 flat, BLAST collapses 0.43→0.71, ProteInfer 0.52 — weak on cc.)

## What the baselines reveal

1. **BLAST is the clean memorization signature:** Fmax collapses with dissimilarity
   (mf 0.47 → 0.75 from <30 to 50–70). Homology transfer is excellent when a similar training
   protein exists and useless when not — the GO analogue of the CLEAN collapse, now shown
   against the *canonical* baseline.
2. **Naive is the flat frequency floor** (mf 0.33). Any predictor near it is "just frequency."
3. **Lab models (MLP, DeepGOZero) beat Naive and show a mild similarity gradient** (0.60→0.73) —
   partial homology dependence, but real signal above frequency.
4. **ProteInfer is a genuine generalizer — the hypothesis that it "just memorizes frequency/
   hierarchy" is REFUTED.** At <30% identity it scores 0.704: far above Naive (0.326), above
   BLAST (0.472), above MLP (0.597). And it is *flat* across similarity. So its strength is
   neither frequency nor homology — it is the only base robust in the low-similarity regime.

## Why this sharpens the WMC thesis (and the NMI claim)

WMC helps where the base **fails**. Two consequences:
- **At low sequence similarity, the best base (ProteInfer) does NOT fail** — it is already
  robust — so WMC has little to add there. (This is why within-vocab GO WMC is neutral.)
- **At zero-shot (held-out terms), EVERY predictor fails** — including ProteInfer — because the
  term is absent from the output entirely; Naive and BLAST are both ≈0.5 (no frequency, no
  homolog). There, training-free WMC reaches **0.80–0.91 AUC**, matching the *trained*
  DeepGOZero-zero (mf 0.9125 vs 0.9118; bp 0.797 vs 0.796, p=2e-4; cc 0.827 vs 0.832).

So the honest, baseline-controlled NMI statement is: **the genuine failure regime is zero-shot,
not low-similarity** (a strong base like ProteInfer already handles low similarity); and in that
regime a *training-free* ontology-WMC layer matches a *trained* zero-shot extrapolator, while the
trivial baselines (Naive, BLAST) cannot score the terms at all. The low-similarity *collapse*
story belongs to BLAST/homology and the lab CNN, not to the best learned base.

Scripts: `go_cafa.py` (Naive+BLAST+learned, Fmax by bucket), `scale_zero_h2h.py` (zero-shot).
