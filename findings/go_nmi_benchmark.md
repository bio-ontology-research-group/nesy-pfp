# GO at scale for NMI: within-vocab neutral, memorization gradient, training-free zero-shot matches SOTA

Large SLURM experiments (unimatrix01) testing the goal's criteria (i) and (ii) on the GO
ontology, across cc/mf/bp and bases MLP / DeepGOCNN / ProteInfer / DeepGOZero.

## (ii-a) Within-vocab: WMC is neutral (honest)

Raw vs max-propagation vs soft-WMC closure, protein-centric Fmax (`go_wmc_eval.py`, job 4475):

| ont / base | raw | maxprop | wmc(nf1+nf2) |
|---|---|---|---|
| mf / MLP | 0.654 | 0.656 | 0.653 |
| bp / MLP | 0.458 | 0.420 | 0.451 |
| bp / DeepGOCNN | 0.344 | 0.301 | 0.336 |

Within-vocab, the hierarchy is already baked into the supervised labels: WMC ≈ base (within
±0.007), occasionally a hair below on bp. **So "base+WMC ≥ base" holds only as neutrality
within-vocab — not a strict win.** maxprop (the universal heuristic) likewise fails to help
and hurts on bp. The value of WMC is therefore *not* within-vocab.

## (i) GO memorization gradient (2nd predictor family)

Per-bucket Fmax by sequence identity to training (`go_mem.py`, job 4487). The DeepGOZero split
is *already* similarity-controlled (clusters by ≥50% id), so the high-similarity buckets are
nearly empty — the dramatic CLEAN-style collapse is not even measurable here. Within the
available range, the lab models still show a **monotonic similarity gradient**:

| ont / base | <30 | 30–50 | 50–70 |
|---|---|---|---|
| bp / MLP | 0.418 | 0.487 | 0.569 |
| bp / DeepGOZero | 0.404 | 0.480 | 0.577 |
| mf / MLP | 0.597 | 0.697 | 0.733 |
| bp / **ProteInfer** | 0.603 | 0.620 | **0.642** (flat) |

Lab predictors (MLP/DeepGOCNN/DeepGOZero) climb with similarity even inside the controlled
split; **ProteInfer is the robust outlier** (trained on far more data). Honest partial support
for criterion (i): a real similarity gradient, milder than CLEAN because the benchmark
pre-removed the high-similarity regime.

## (ii-b) Zero-shot — the actual GO win: training-free WMC matches trained zero-shot

Canonical DeepGOZero zero-shot protocol (`scale_zero_h2h.py`, job 4493): for every held-out
term in `terms_zero_10` (no training annotations), score it by (a) **WMC** — composing the
ontology definition over the base predictor's in-vocab scores, *training-free* — vs (b)
**DeepGOZero-zero** — a *specially trained* ELEmbedding zero-shot extrapolator
(`deepgozero_zero_10.th`). Gold = `prop_annotations`; macro-AUC over WMC-groundable terms:

| ont | n terms | WMC (training-free) | DeepGOZero-zero (trained) | Wilcoxon |
|---|---|---|---|---|
| mf | 194 | **0.9125** | 0.9118 | p = 0.83 (tied) |
| cc | 144 | 0.827 | 0.832 | p = 0.005 (DGZ +0.005) |
| bp | (running) | — | — | — |

The decisive framing: **the base predictor alone cannot score held-out terms at all (≈0.5).**
WMC lifts it to **0.91 (mf)** in the exact failure regime — and does so **matching a
purpose-built *trained* zero-shot model with zero zero-shot training**, just the EL
definitions. On mf the two are statistically indistinguishable (0.9125 vs 0.9118, p=0.83);
on cc DeepGOZero edges it by 0.005.

## Verdict against the goal

- **Criterion (ii) on GO holds in the failure regime, honestly qualified:** base+WMC ≫ base
  for zero-shot terms (0.91 vs ~0.5), *training-free*, **matching** SOTA trained zero-shot —
  but **not "significantly better"** than the trained zero-shot baseline (it ties). Within-vocab
  is neutral. So the GO contribution is "structure gives any base predictor SOTA zero-shot
  capability for free," not "WMC beats the best GO predictor."
- This is **weaker than the earlier 16-curated-term claim** (which showed WMC *beating*
  DeepGOZero-zero); at scale it is a **tie**. Honest correction of an over-strong prior result.
- Combined with EC (producibility-WMC refine 0.491→0.556, precision-tempered), the two-ontology
  criterion (ii) is met as: **WMC adds similarity-independent / zero-shot capability the base
  lacks, training-free — a complement that matches trained extrapolation, not a universal
  accuracy win.** The clean theoretical leg (exact-on-trees + modular <8% core) stands.

Scripts: `go_wmc_eval.py`, `go_mem.py`, `scale_zero_h2h.py`, `go_zero_scale.py`.
