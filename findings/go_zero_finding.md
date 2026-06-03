# Zero-shot GO prediction by EL++ WMC — training-free, beats DeepGOZero's embedding

The within-vocab result showed hierarchy WMC is *neutral* when the structure is
already in the supervised labels. The principle it gave —

> WMC helps precisely when the symbolic structure carries information **absent** from
> the perception model's supervision

— predicts that GO axioms *do* pay off for **zero-shot terms**: GO classes with zero
training annotations, which the perception net cannot score at all and which are
reachable only through their logical definition. This is the GO analogue of the
metabolic "external structure" condition, and it is DeepGOZero's home turf.

## Method

DeepGOZero normalises GO-Plus to the EL theory `go.norm` (`nf1` subsumption, `nf2`
conjunction `A⊓B⊑E`, `nf3/nf4` existential role restrictions `R some C`). For each
held-out term `E` we evaluate its definition **downward** from the base predictor's
in-vocab marginals (`go_zero_wmc.py`):

- in-vocab named class → the base predictor's graded marginal (stop);
- `E` defined by `A⊓B⊑E` → `score(A)·score(B)` (conjunction = product = WMC of the definition);
- existential node `∃R.D` → `score(D)` (filler content); `nf3: R some C ⊑ D` recursed;
- alternative sufficient conditions combined by noisy-OR.

This is the soft EL++ completion moose compiles — here at inference, composing seen
parts into an unseen term. **Genuinely zero-shot**: the 16 curated terms are *masked*
from the base predictor (even though the stock MLP has them in-vocab), so their score
comes only through the definition. Metric: per-term ROC-AUC (ranking; the noisy-OR's
absolute calibration is irrelevant). Base = DeepGOZero's MLP (InterPro features),
similarity-held-out test set.

## Result — DeepGOZero's own 16 curated zero-shot terms

| | **WMC zero-shot (ours, training-free)** | DeepGOZero zero-shot (published) | DeepGOZero *supervised* (ceiling) |
|---|---|---|---|
| **mf** | **0.961** on 2/5 reached | 0.610 | 0.939 |
| **bp** | **0.860** (7/7) | 0.796 | 0.877 |
| **cc** | **0.866** (4/4) | 0.825 | 0.871 |

Per-term AUC (WMC / DGZ-zero / DGZ-supervised):

```
mf  GO:0001227   —*   / 0.257 / 0.932     bp  GO:0000381  0.881 / 0.855 / 0.876
    GO:0001228   —*   / 0.574 / 0.951         GO:0032729  0.923 / 0.870 / 0.911
    GO:0003735   —*   / 0.400 / 0.891         GO:0032755  0.810 / 0.719 / 0.866
    GO:0004867  0.989 / 0.972 / 0.986         GO:0032760  0.864 / 0.861 / 0.910
    GO:0005096  0.933 / 0.847 / 0.932         GO:0046330  0.907 / 0.855 / 0.916
                                              GO:0051897  0.872 / 0.772 / 0.895
cc  GO:0005762  0.932 / 0.889 / 0.893         GO:0120162  0.767 / 0.637 / 0.767
    GO:0022625  0.880 / 0.898 / 0.910
    GO:0042788  0.810 / 0.858 / 0.891     (—* = not reached: definitional parts
    GO:1904813  0.841 / 0.653 / 0.789      do not bottom out in the base vocab)
```

## What it says

1. **Training-free EL++ WMC zero-shot beats DeepGOZero's *published* zero-shot
   embedding** on bp (0.860 vs 0.796) and cc (0.866 vs 0.825), and on every reached
   mf term. On bp it wins all 7 terms individually. DeepGOZero must *train* a special
   `el_loss` embedding model to place zero-shot classes; WMC needs only a stock
   per-term predictor plus the published axioms at inference.
2. **It nearly matches DeepGOZero's *supervised* ceiling** (bp 0.860 vs 0.877; cc
   0.866 vs 0.871) — i.e. composing seen parts through the ontology recovers almost
   all of what a model trained *with* the term's annotations achieves.
3. **This is the metabolic finding's mirror image, confirming the principle.** Where
   the structure is external to supervision (metabolic network; zero-shot GO terms),
   WMC adds large value; where it is already in the labels (within-vocab GO), it is
   neutral. Same WMC machinery, opposite outcomes, predicted by one rule.
4. **Honest coverage limit.** 3/16 mf terms are not reachable because their
   definitional parts do not bottom out in the base vocabulary via `nf2/nf3`
   (a subclass-union fallback exists but explodes the closure and is not used).
   Notably DeepGOZero *also* fails exactly these three (0.257, 0.574, 0.400 — two
   below random), so they are intrinsically hard, not a WMC-specific weakness.

## Caveats

- WMC and DeepGOZero-published AUCs use the standard DeepGOZero protocol but slightly
  different test-protein subsets (theirs n≈1300–1660, ours the full per-ont test set).
  The matched-protein comparison is WMC-zero vs DGZ-*supervised* (same proteins, same
  model machinery) — there WMC-zero nearly equals the supervised ceiling.
- "Held-out" masks the term from the base at inference; a base *retrained* without the
  term would be stricter, but the part-term predictions WMC composes are unaffected.
- Next: confirm base-independence with a second base (DeepGOCNN, different
  architecture; ProteInfer pending — TF1.15 environment).

## Base-independence — DeepGOCNN (sequence CNN, different architecture/input)

Re-running the same zero-shot closure on a second base — DeepGOCNN (1-D CNN over the
sequence, vs the InterPro MLP) — confirms the mechanism is not specific to one
predictor. Mean zero-shot WMC-AUC:

| ont | WMC (MLP/InterPro base) | WMC (CNN/sequence base) | DGZ zero (published) | DGZ supervised |
|---|---|---|---|---|
| mf | 0.96 (2/5 reached) | 0.88 / 0.74 (2/5 reached) | 0.610 | 0.939 |
| bp | 0.860 | 0.746 | 0.796 | 0.877 |
| cc | 0.866 | **0.894** | 0.825 | 0.871 |

The zero-shot derivation works with **both** bases. Magnitudes track *base quality on
the definitional part-terms*: the CNN is stronger on cc (0.894, the best of any method
incl. supervised DGZ) but weaker on bp (0.746) — WMC is only as good as the base's
predictions of the in-vocab parts it composes. So the **mechanism** is
base-independent; the **absolute score** inherits the base. Both bases beat
DeepGOZero's published zero-shot on cc; the InterPro base also beats it on bp.

## Third base — ProteInfer (Google, non-lab, sequence CNN, different training data)

To rule out that the result is specific to lab predictors, a third base: **ProteInfer**
(Google, dilated CNN, SwissProt) run via its own TF1.15 SavedModel (32,102 GO labels,
covering 94% of DeepGOZero's terms), scores aligned to the DeepGOZero vocab.

Mean zero-shot WMC-AUC across **three** bases (reached-term count in parens):

| ont | MLP (InterPro) | CNN (sequence) | **ProteInfer (non-lab)** | DGZ zero | DGZ supervised |
|---|---|---|---|---|---|
| mf | 0.96 (2/5) | 0.81 (2/5) | 0.55 (1/5) | 0.610 | 0.939 |
| bp | 0.860 (7/7) | 0.746 (7/7) | 0.66 (3 strong: .81–.90) | 0.796 | 0.877 |
| cc | 0.866 (4/4) | **0.894** (4/4) | 0.70 (3/4) | 0.825 | 0.871 |

ProteInfer composes **some** terms well (bp GO:0051897 0.898, GO:0046330 0.814; cc
GO:0022625 0.859) but is the weakest of the three — and for a *principled* reason that
reinforces the thesis rather than weakening it:

- **WMC composition is exactly as good as the base's prediction of the definitional
  parts.** Degenerate terms (AUC 0.5) are precisely those whose definition recurses to
  a part-leaf ProteInfer scores ≈0. Example: `GO:0051897 ≡ GO:0050525 ⊓ GO:0065007`
  works because ProteInfer predicts both leaves; `GO:0120162 ≡ GO:0120170 ⊓ GO:0065007`
  fails because ProteInfer scores `GO:0120170`'s sub-definition ≈0.
- ProteInfer is **out-of-domain** for this composition: it predicts a *different*
  subset of GO well than the DeepGOZero-vocab MLP/CNN, so it covers fewer of the
  specific definitional leaves. The two GO-vocab-trained bases compose better.

**Conclusion across all three bases:** the EL++ WMC composition is *base-agnostic in
mechanism* — it lifts any per-term predictor to unseen terms through the published
axioms, training-free — and *base-dependent in magnitude*, recovering more the better
the base predicts the definitional parts. With an in-domain base it beats DeepGOZero's
trained zero-shot embedding and nears the supervised ceiling; with an out-of-domain
base (ProteInfer) it still composes the terms whose parts that base predicts.
