# Precision / blast-radius of the producibility-WMC refinement (honest tempering)

The refinement result (`clean_prior_wmc_result.md`) measured *recall* of held-out true
reactions + zero displacement of CLEAN-correct calls. The open caveat was the **precision
cost** of the reactions the WMC *adds* to restore producibility. Measured here (job 4457,
225 reactions, product-target CLEAN-weighted fills, all added reactions captured).

## Ground-truth proxy
A WMC-added reaction is a **true positive** if its EC is in the organism's homology-supported
repertoire (any reaction with a genome BLAST hit in the `all-Reactions.tbl`), or is the
held-out true EC. Otherwise a **false positive**.

## Result

| identity | blast radius (mean added) | added-EC precision |
|---|---|---|
| ≤30% | 50.5 | 0.697 |
| 30–50% | 94.1 | 0.674 |
| 50–70% | 85.2 | 0.684 |
| 70–90% | 76.5 | 0.705 |
| 90–100% | 51.4 | 0.708 |
| **overall** | **71.6** (median 69) | **0.691** (TP 8210 / FP 3674) |

Among the 32 actual EC-recoveries: mean **103.5** reactions added, median 116, and **0/32**
recovered with ≤3 additions.

## What it means (honest)

1. **Not a clean Pareto improvement.** Restoring producibility of a single held-out product
   from minimal media pulls in whole pathways — ~70–100 reactions per recovery — of which
   **~31% carry ECs the genome does not support by homology**. The recall gain
   (refine 0.491→0.556) is real but comes with a substantial false-positive rate.

2. **The EC-level recovery metric is partly inflated.** With a blast radius of ~100, hitting
   the held-out EC is in part coincidental coverage, not surgical recovery. The honest read
   of the +6.5-pt lift is "the WMC proposes a large candidate set, similarity-independent,
   that *contains* the missed enzyme 14 extra times — at 69% per-reaction precision," not
   "the WMC pinpoints the missed enzyme."

3. **Method, not verdict.** The blast radius is a property of parsimony gap-fill toward a
   single distal metabolite from minimal media. A surgical variant (target the immediate
   precursor; penalise additions; or rank additions by the CLEAN prior and threshold) should
   shrink it and raise precision — a concrete methodological next step, not a dead end.

## Consequence for NMI / eLife framing

The defensible claim is now precisely bounded: producibility-WMC on top of CLEAN is a
**similarity-independent candidate-generation / hypothesis layer** that recovers enzymes the
sequence model misses, at ~0.69 precision and a large blast radius — useful for discovery,
not a high-precision drop-in. Pairing this honest precision number with the recall lift is
stronger (and more referee-proof) than reporting recall alone. Artifacts:
`partb/precision_recover.py`, `partb/prec_agg.py`, `prec_compare.pkl`.

## Surgical-precision lever: FAILED (honest negative, job 4499 + min-growth probe)

Attempted to shrink the ~72-reaction blast radius by (a) adding the held-out reaction's
immediate substrates to the medium and (b) lowering `--min-growth`. Both failed:
- Substrate-augmented medium: blast radius **72.1** (median 100), precision **0.699** —
  unchanged from the distal-target 71.6 / 0.691.
- `--min-growth` sweep on a probe reaction: 0.01 → 151 added; 0.0001 → 149; **0 → 0 added**.

**Mechanism (fundamental, not a tuning artifact):** gapsmith `fill --target <cpd>` does not
isolate "make this metabolite producible" — it gap-fills the model to **grow** (biomass), with
the target folded into the objective. The entire blast radius is the *growth* gap-fill; the
held-out reaction is recovered only when it lies on the growth-essential path (hence the ~14%
recovery and ~0.69 precision). Disabling growth (`min-growth 0`) adds nothing at all, so there
is no growth-free "produce metabolite X minimally" mode to exploit. **The candidate-generation
character of producibility-WMC is therefore intrinsic** — precision cannot be raised by media
or growth tuning within gapsmith's design. Raising it would require a different solver mode
(minimal reaction set to produce a target, growth-decoupled) that gapsmith does not provide.
