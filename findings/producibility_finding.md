# Producibility-as-selection — finding (the precision mechanism)

Follow-up to `coupled_finding.md`. The coverage coupling was recall-leaning
and *amplified* false positives. This experiment tests the METEOR-style
**selection** objective in our differentiable-WMC regime, posed as ASP-style
**producibility** rather than MILP flux balance.

## Setup (exact, brute-force WMC — `producibility_prototype.py`)

Atoms `RXN|r` (reaction present) and `MET|m` (metabolite producible).
A metabolite is producible iff reachable from the growth medium (seeds) by
forward firing of *present* reactions — the monotone least-fixpoint /
scope semantics (Schaub–Thiele expansion, meneco-style gap-filling), which
is Horn / deterministic-closure and thus compilable to a tractable circuit.
The biomass target metabolites are observed producible. Each reaction
carries a Bernoulli weight `σ(logit(perception_r) − β)`; `β` is a sparsity
("off by default") prior, the analogue of METEOR's parsimony / cost term.
The loss is `−log WMC`; we read off the per-reaction posterior marginal and
the numerical gradient `d(−logWMC)/d(perception_r)`.

## Results

**1. Necessitated reactions are pinned ON and trained up (recall/gap-fill).**
On `A→B→C→Z` with target `Z`, reactions `r1,r2,r3` get posterior 1.000 and
strong push-up gradient (−5.0 at β=2). This is the gap-filling direction.

**2. No false positive is ever amplified** (the coverage failure does not
recur). A topologically dead FP (`rDead`, substrate never producible) and a
reachable-but-unneeded FP (`rUnneeded`, METEOR's "silent case") both get
gradient **exactly 0** and posterior shrunk by the prior (0.90→0.55 at β=2).
This holds under an 8-decoy width stress with **no underflow** — coverage's
all-false product (`coupled_finding.md`) is gone.

| reaction class | posterior (β=2) | d(−logWMC)/dp | vs. coverage |
|---|---|---|---|
| necessitated (r1–r3) | 1.000 | −5.0 (push up) | up (same) |
| dead-end FP | 0.55 | 0.0 | coverage pushed **up** |
| unneeded-reachable FP | 0.55 | 0.0 | coverage pushed **up** |

**3. Redundant catalysts: the weaker candidate is suppressed
(explaining-away / soft at-most-one).** Two enzymes `r2a` (perception 0.9,
true) and `r2b` (0.6, FP) both catalyse the one reaction the target needs.
At β=1.5 the weaker `r2b` posterior **shrinks 0.60→0.334** while `r2a` holds
at 0.889 and gets **7× the push-up gradient** (−2.46 vs −0.35). The model
prefers to satisfy the requirement through the stronger catalyst and lets
the redundant weaker one fall — the per-reaction false-positive suppression
that coverage could not produce.

## Honest scope (matches METEOR)

The precision gain is **targeted, not universal** — exactly METEOR's
character (its own ablation, `metabolic-model/clean_main.tex:542`, reports a
modest ROC-AUC 0.697→0.713 even at 90% genome removal). The objective bites
on (a) reactions whose removal breaks biomass (essential) and (b) redundant
catalysts for a required reaction. A reachable FP that is simply never
*needed* is left to the uniform sparsity prior, not discriminated — METEOR's
"silent case" (`clean_main.tex:706`). Full mass-balance (the MILP) catches
some of those; producibility is the Boolean relaxation that lands inside our
differentiable, modular WMC.

## Why ASP, not MILP (the framing)

METEOR's MILP couples binary selection to **continuous steady-state flux**
(`Sv=0`, `v_biomass≥γ`) — real-valued linear algebra that cannot be a WMC
circuit. ASP producibility keeps only the **qualitative reachability**, a
monotone Horn fixpoint = moose's home regime, and genome-scale producibility
is one large fixpoint that motivates the **modular** WMC (subsystem modules,
boundary atoms = shared currency metabolites). The contribution is
**METEOR made differentiable**: instead of a post-hoc MILP filter on fixed
probabilities, the producibility constraint becomes a WMC layer that trains
the perception net end-to-end (what DPL's independence assumption could not
do). Use `clingo` as an *offline oracle* (ground the program, enumerate
minimal completions, validate foundedness on cyclic networks via stable-model
semantics) — the trainable object stays the WMC circuit.

## Real E. coli (BiGG e_coli_core) + ESM-2 perception

`moose/nesy/producibility.py` + `experiments/producibility_ecoli.py`.
e_coli_core: 75 internal reactions, 72 metabolites, 16 biomass precursors,
genes carry direct UniProt ids (GPR -> ESM-2). The differentiable
`SelectionWMC` layer enumerates the feasible reaction configurations once
(exact, ≤ ~20 candidate reactions) and forms a vectorised log-WMC.

**Bootstrap fix (a real scope-analysis issue).** PTS glucose uptake
`GLCpts: glc__D_e + pep_c -> g6p_c + pyr_c` is autocatalytic in `pep_c`;
network expansion cannot bootstrap it, so with only the medium seeded *no*
biomass precursor is reachable. Seeding `pep_c` as a primer (standard in
scope analysis / meneco) unlocks all 16 — this is exactly the foundedness
concern the ASP literature handles, surfacing on real data.

**The FP-suppression framing breaks on real data — and that is informative.**
Every e_coli_core reaction is a genuine E. coli reaction, so "inject decoys"
is ill-posed, and (matching the prototype + METEOR's silent case) likelihood
WMC only *prior-shrinks* unneeded reactions: with ESM the perception
saturates to 1.0 and the β prior barely moves the read-out. Precision via WMC
is weak/targeted, full stop.

**The robust, defensible effect is completion (gap-filling).** Hide the
annotation of the necessary reactions (train perception toward 0 on them) and
ask whether the producibility-selection *recovers* them. Target `3pg_c`
(necessary: ENO, PGM), β=3, mask 50 %, 3 seeds:

| arm | perception on hidden **necessary** | perception on hidden **unnecessary** |
|---|---|---|
| base (perception only) | 0.000 | 0.000 |
| + producibility WMC | **0.994** | 0.00–0.06 |

The selection loss trains the **perception network** to recover the
necessary-but-unannotated reactions (0.000 → 0.994, robust across seeds)
**without hallucinating** the unnecessary ones (≤0.06). This is differentiable
gap-filling — METEOR's idea made trainable, and it composes with BP taxon
constraints (`taxon_finding.md`).

**Honest scope.** The constraint can only pin what biomass *necessitates*; a
metabolite with redundant production routes (e.g. `accoa_c`) has no strictly
necessary reaction, so no recovery signal — again targeted, matching METEOR.
Collective necessity over the *full* 16-precursor biomass (where many
reactions become jointly required) needs all 75 reactions, i.e. the modular
SDD compile — the next phase and the moose contribution proper.

## Founded SDD encoding + tractability frontier

`build_bounded_theory` (in `producibility.py`) gives a **founded** bounded-step
encoding compilable to an SDD via moose's `compile_theory`: `M|m@k` is
producible-by-step-k iff producible by k-1 or made by a present reaction whose
inputs are all producible-by-(k-1). The level index makes the fixpoint acyclic,
so unlike a flat Clark completion it is correct on cyclic networks — a B↔C
cycle does **not** self-support (`test_producibility_network.py`). Reaction
atoms carry the perception weight; metabolite atoms are functionally determined
and sum out. **Validated:** the SDD model count equals the exact least-fixpoint
feasibility count (e.g. 3pg_c subsystem: 128 = 128).

Tractability scan (e_coli_core, accoa_c subnetworks, cofactor+pep_c seeded):

| reactions | step-depth K | SDD size | compile | exact-match |
|---|---|---|---|---|
| 9  | 3 | 2,284  | 0.07s | ✓ |
| 15 | 3 | 2,938  | 0.19s | ✓ |
| 20 | 3 | 3,094  | 0.26s | ✓ |
| 25–30 | 3 | ~3,000 | 0.27s | (>20, exact infeasible) |
| 31 | **6** | **378,126** | **38s** | — |
| full 75 | — | — | timeout | — |

The cost driver is **step-depth K**, not reaction count: the metabolite layer
is copied K times and the biconditional completion couples the copies, so the
SDD explodes once the producing network gets deep (K≥6). A monolithic
full-biomass compile times out. This is the data that motivates the modular /
SCC-stratified compile (the moose contribution): keep each module shallow and
couple subsystems only through shared boundary metabolites.

## Full-biomass producibility compiles (no modular decomposition needed at
## e_coli_core scale) — and recovers jointly-necessary reactions

The K-layering blow-up is from **reversibility**: all-reversible gives huge
SCCs ([27, 15, 2]) and the SDD times out. But only a *minimal* reversible set
is actually needed to keep all 16 biomass precursors reachable — a greedy
search finds just **3** reactions (`GLUDy`, `PGK`, `RPI`); the rest can be
forward-only. With that min-rev reduction the full theory (74 reactions, 16
targets, true synchronous K=8, 749 atoms) **compiles to an SDD of 13,740 nodes
in ~140–200 s**, model count 1.08×10¹⁷ (validated non-degenerate). So at
e_coli_core scale the monolithic founded compile suffices; modular
decomposition is reserved for genome-scale (iML1515). (Metabolite-major
variable reordering was *worse* than the default — default ordering wins.)

`SddSelectionWMC` then gives the differentiable selection loss over the cached
circuit (reactions weighted, metabolites neutral, targets clamped). Full-biomass
completion (hide 45 of 74 annotations incl. all necessary; ESM-2 perception;
β=2; `producibility_fullbiomass.py`):

| arm | recover **necessary** (17) | recover redundant (28) |
|---|---|---|
| base (perception only) | 0.002 | 0.061 |
| + producibility WMC | **0.998** | 0.475 |

**17/74 reactions are jointly necessary** for the 16 precursors (vs 2 in the
single-target run) — real collective necessity. The selection WMC recovers
**all 17 (0.002 → 0.998)** where perception alone recovers nothing, and
partially recovers the useful-redundant reactions (→0.475). Note every
e_coli_core reaction is a *genuine* E. coli reaction, so raising redundant ones
is correct completion, **not** hallucination — e_coli_core has no true
negatives, so precision (not asserting *absent* reactions) can only be measured
with genuinely-absent decoys. That is the wide-candidate experiment next.

## Wide-candidate precision test (true negatives) — recall vs precision split

Candidates = e_coli_core (min-rev) + 20 decoy reactions from the **yeast**
model iMM904 that E. coli does not have (genuine true-negatives, each sharing
≥1 metabolite with the core). All 20 are dead-ends in the core network (they
need yeast-specific substrates) — i.e. free in the WMC. Masked-core annotations
are supervised toward 0 (so recall is visible); decoys are unannotated and
over-predicted. β=2, `producibility_precision.py`:

| arm | recover necessary | absent decoys (mean) | AUROC(core>decoy) |
|---|---|---|---|
| base (perception only) | 0.002 | 0.881 | 0.365 |
| + producibility | **0.998** | 0.881 | 0.554 |
| + producibility + taxon (`¬R`) | **0.998** | **0.252** | **0.946** |

The split is clean and is the thesis of the whole metabolic strand:

* **Producibility is a recall / completion constraint.** It recovers all 17
  jointly-necessary reactions (0.002 → 0.998) but leaves the absent decoys at
  0.881 — a dead-end reaction is free in the WMC, so the selection gradient on
  it is ~0; producibility *cannot* suppress false positives.
* **Taxon is the precision constraint.** The hard `¬R` exclusion drives the
  absent decoys' perception down (0.881 → 0.252) *through the selection loss*
  (the unit clause makes models with the reaction present worthless, so the
  gradient pushes its probability down — not merely a read-out mask).
* **Together** they give recall *and* precision: AUROC(true core vs absent)
  0.365 → 0.946.

Caveat: the decoys here are all dead-ends, the easy case for taxon. Connectable
decoys (which producibility would *raise*, being recall-leaning) are the harder
case and the reason taxon — not producibility — must carry precision.

## Next

- **Modular / SCC-stratified founded compile** of the full biomass theory:
  Clark completion (no step index) on the DAG part — most of the
  cofactor-seeded network is shallow — and bounded-step unrolling only *within*
  strongly-connected components (the genuine cycles), partitioned by subsystem
  with shared metabolites as boundary atoms + Path B sketching for wide
  boundaries. clingo as the offline foundedness oracle. This is the concrete
  fix for the K-layering blow-up above.
- Wide candidate set including genuinely-absent reactions (iML1515 \ core,
  eukaryote-only) so BP taxon constraints (8,725-term lever) and producibility
  selection act together on real false positives.
