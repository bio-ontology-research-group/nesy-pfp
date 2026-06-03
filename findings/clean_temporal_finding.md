# CLEAN's accuracy is similarity-to-training-driven — temporal, controlled evidence

The de-memorization finding showed CLEAN's headline owed to exact-sequence overlap.
The objection: that's an overlap argument, and CLEAN's own benchmarks (New-392) are
themselves close to training. This closes it with the **gold-standard test**: enzymes
added to SwissProt **after CLEAN was published**, with **similarity to CLEAN's training
controlled by us**.

## Setup

- **Novel enzymes (we collect them):** UniProt reviewed entries with an EC number,
  `date_created ≥ 2023-01-01` (CLEAN trained on ~2022 SwissProt) → **2,744 enzymes,
  2,050 with protein-level (experimental) evidence**. **0 are in CLEAN's training set
  by accession** — all provably unseen.
- **Similarity control (we compute it):** diamond blastp each novel enzyme against
  CLEAN's actual training file (`split100.csv`, 227,362 enzymes) → max % identity to
  any training protein. Bucket by that identity. The novel set is **75% at ≤50%
  identity** — well-populated dissimilar regime.
- **CLEAN run:** its own `split100` pretrained model + GMM, ESM-1b embeddings, maxsep
  inference (node005 GPU via SLURM). Accuracy = predicted EC matches a true UniProt EC.

## Result — CLEAN EC accuracy vs identity-to-its-own-training (experimental enzymes)

| identity to CLEAN training | n | exact-EC | EC level-3 |
|---|---|---|---|
| 90–100% | 69 | **0.855** | 0.942 |
| 50–90% | 428 | 0.472 | 0.841 |
| 30–50% | 751 | 0.252 | 0.699 |
| ≤30% | 226 | 0.097 | 0.487 |
| no detectable hit | 555 | **0.070** | 0.378 |

(ALL-novel n=2,744 gives the same shape: 0.861 → 0.072 exact-EC.)

## What it says

- **CLEAN's performance is overwhelmingly a function of similarity to its training**, a
  **monotonic ~12× collapse** (0.86 → 0.07 exact-EC) from near-identical to no-homolog
  novel enzymes. On genuinely novel, dissimilar enzymes (≤30% identity, 781
  experimental cases) it predicts the correct EC ~7–10% of the time.
- This is **temporal + similarity-controlled**, on enzymes CLEAN could not have seen,
  with the control computed by us against CLEAN's real training file — not CLEAN's own
  (contaminated) New-392 benchmark. It converts "CLEAN's edge is memorization" from an
  overlap argument into a controlled generalization-failure curve.
- **Goal criterion (i) met:** CLEAN accuracy drops monotonically high→low identity.

## Next (Part B — the neuro-symbolic claim)

Show that a perception model refined through **organism-specific metabolic-producibility
WMC** retains reaction-prediction accuracy in exactly the low-similarity buckets where
CLEAN fails. Genome-scale models built per organism with **gapsmith** (`doall`);
bacterial targets with novel low-similarity enzymes: *E. coli, M. tuberculosis,
P. aeruginosa, Streptomyces, A. baumannii*. All compute via SLURM on node005.

Scripts (on cluster `~/nesy-genome`): `novel_bucket.sbatch` (UniProt→diamond buckets),
`clean_novel.sbatch` (CLEAN inference), `clean_bucket_eval.py` (accuracy by bucket).
CLEAN needs a one-line pandas-compat patch (`smallest_10_dist_df[i]` → `dist_lst[i]`)
under the moose venv.
