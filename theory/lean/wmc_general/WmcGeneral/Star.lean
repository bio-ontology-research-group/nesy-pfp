import Mathlib

set_option linter.style.header false
set_option linter.style.longLine false

/-!
  WmcGeneral/Star.lean — the GENERAL (all-priors) tree-exactness theorem.

  Amends the Lean-core development (`lean-wmc-approx`, `lean-wmc-approx-general`)
  with the universal statement those files witness by `decide`: for a STAR node
  (a parent with `k` independent children under the true-path constraint
  `child ⟹ parent`), the exact weighted-model-count marginal of the parent equals
  the soft-OR / belief-propagation update  `qp / (qp + (1 - qp) * ∏ i, (1 - q i))`,
  for ALL `k` and ALL priors. This is the inductive primitive of "BP = exact WMC on
  the tree fragment"; reconvergence (shared ancestors) is what breaks it.
-/

open Finset BigOperators

namespace WmcGeneral

variable {k : ℕ} (qp : ℝ) (q : Fin k → ℝ)

/-- Weight of a world `(p, c)`: parent state `p`, children states `c`. -/
def w (p : Bool) (c : Fin k → Bool) : ℝ :=
  (if p then qp else 1 - qp) * ∏ i, (if c i then q i else 1 - q i)

/-- Unnormalised mass of valid worlds with the parent TRUE (all children valid). -/
noncomputable def num : ℝ := ∑ c : Fin k → Bool, w qp q true c

/-- Partition function: parent-true mass plus the unique valid parent-false world. -/
noncomputable def Zp : ℝ := num qp q + w qp q false (fun _ => false)

/-- Sum over all children assignments of the product factorises (independence). -/
lemma sum_children :
    (∑ c : Fin k → Bool, ∏ i, (if c i then q i else 1 - q i)) = 1 := by
  have h : (∏ i, ∑ b : Bool, (if b then q i else 1 - q i))
         = ∑ c : Fin k → Bool, ∏ i, (if c i then q i else 1 - q i) := by
    rw [Finset.prod_univ_sum]; rw [Fintype.piFinset_univ]
  rw [← h]
  apply Finset.prod_eq_one
  intro i _
  simp only [Fintype.sum_bool, Bool.false_eq_true, if_true, if_false]
  ring

theorem num_eq : num qp q = qp := by
  unfold num w
  simp only [if_true]
  rw [← Finset.mul_sum]
  rw [sum_children]
  ring

theorem Zp_eq : Zp qp q = qp + (1 - qp) * ∏ i, (1 - q i) := by
  unfold Zp w
  rw [num_eq]
  simp only [Bool.false_eq_true, if_false]

/-- **General tree-exactness (star primitive).** -/
theorem star_marginal_eq_softOR :
    num qp q / Zp qp q = qp / (qp + (1 - qp) * ∏ i, (1 - q i)) := by
  rw [num_eq, Zp_eq]

end WmcGeneral
