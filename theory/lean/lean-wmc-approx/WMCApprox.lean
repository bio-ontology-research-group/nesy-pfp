/-
  WMCApprox.lean — when does the soft "true-path" fixpoint equal exact WMC?

  We model the GO true-path constraint (child ⟹ parent) as a Horn theory and
  compare the EXACT weighted model count (WMC) marginal against the belief-
  propagation / noisy-OR approximation the GO code uses.  BP = exact inference
  on the *computation-tree unrolling*; so BP = WMC iff the constraint graph is
  already a tree, and BP ≠ WMC exactly at a reconvergence (diamond).

  All arithmetic is over a kernel-reducible rational `Q = (num,den)`; every
  witness closes by `decide`, so the file's axioms stay within the Lean 4
  standard set {propext, Classical.choice, Quot.sound} (no `native_decide`).
-/

-- ===== kernel-reducible rationals =====
abbrev Q := Int × Int
def Q.add (a b : Q) : Q := (a.1 * b.2 + b.1 * a.2, a.2 * b.2)
def Q.mul (a b : Q) : Q := (a.1 * b.1, a.2 * b.2)
def Q.compl (a : Q) : Q := (a.2 - a.1, a.2)             -- 1 - a   (den unchanged)
def Q.eq (a b : Q) : Bool := decide (a.1 * b.2 = b.1 * a.2)
infixl:65 " +q " => Q.add
infixl:70 " *q " => Q.mul

-- ===== worlds and the true-path constraint =====
-- variables are positions 0..n-1; a world is a Bool list of length n.
def get (w : List Bool) (i : Nat) : Bool := (w[i]?).getD false

-- all bitstrings of length n
def worlds : Nat → List (List Bool)
  | 0 => [[]]
  | n+1 => (worlds n).flatMap (fun w => [false :: w, true :: w])

-- edges are (child, parent): child true ⟹ parent true
def valid (edges : List (Nat × Nat)) (w : List Bool) : Bool :=
  edges.all (fun e => (! get w e.1) || get w e.2)

-- weight of a world = ∏ over vars of (prior i if true else 1 - prior i)
def wWeight (priors : List Q) (w : List Bool) : Q :=
  (List.range priors.length).foldl
    (fun acc i =>
      acc *q (if get w i then priors.getD i (0,1) else Q.compl (priors.getD i (0,1)))) (1,1)

-- sum of a Q-list
def qsum (l : List Q) : Q := l.foldl (· +q ·) (0,1)

-- EXACT WMC marginal P(var t = true | constraint)
def wmcMarg (t : Nat) (n : Nat) (edges : List (Nat × Nat)) (priors : List Q) : Q :=
  let val := (worlds n).filter (valid edges)
  let num := qsum ((val.filter (fun w => get w t)).map (wWeight priors))
  let den := qsum (val.map (wWeight priors))
  (num.1 * den.2, num.2 * den.1)            -- num / den, as a Q


-- ===== BP = exact WMC on the computation-tree UNROLLING =====
-- DIAMOND  (child,parent):  d⟹b, d⟹c, b⟹a, c⟹a   [a=0,b=1,c=2,d=3]
--   a is a shared ancestor of d via two paths (b and c) → reconvergence.
def diamondEdges : List (Nat × Nat) := [(3,1),(3,2),(1,0),(2,0)]
def diamondPriors : List Q := [(1,2),(1,2),(1,2),(1,2)]   -- all 1/2

-- UNROLLED diamond: duplicate shared ancestor a into a1,a2 (BP's tree)
--   [a1=0,a2=1,b=2,c=3,d=4]:  d⟹b,d⟹c, b⟹a1, c⟹a2
def unrollEdges : List (Nat × Nat) := [(4,2),(4,3),(2,0),(3,1)]
def unrollPriors : List Q := [(1,2),(1,2),(1,2),(1,2),(1,2)]

-- POLYTREE: two-parent sink whose parents are INDEPENDENT roots
--   [b=0,c=1,d=2]:  d⟹b, d⟹c   (no shared ancestor)
def polyEdges : List (Nat × Nat) := [(2,0),(2,1)]
def polyPriors : List Q := [(1,2),(1,2),(1,2)]
-- its (trivial) unrolling is itself.

-- ===== the BP "independence" form for a two-parent sink (true-path) =====
-- P(sink) = p_d·p_b·p_c / (p_d·p_b·p_c + (1-p_d))   — valid iff parents independent.
def bpSink (pd pb pc : Q) : Q :=
  let n := pd *q pb *q pc
  let den := n +q (Q.compl pd)
  (n.1 * den.2, n.2 * den.1)

def half : Q := (1,2)

-- (1) POLYTREE: independent parents → BP form is EXACT  (1/5 = 1/5)
theorem polytree_sink_exact :
    Q.eq (wmcMarg 2 3 polyEdges polyPriors) (bpSink half half half) = true := by decide

-- (2) DIAMOND: parents share an ancestor → same BP form is WRONG  (exact 1/6 ≠ BP 1/5)
theorem diamond_sink_bp_wrong :
    Q.eq (wmcMarg 3 4 diamondEdges diamondPriors) (bpSink half half half) = false := by decide

-- (3) DIAMOND reconvergent node a: exact 5/6 ≠ BP (computation-tree unrolling) 7/10
theorem diamond_reconv_inexact :
    Q.eq (wmcMarg 0 4 diamondEdges diamondPriors)
         (wmcMarg 0 5 unrollEdges unrollPriors) = false := by decide

-- exact concrete values, for the prose
theorem diamond_exact_a_is_5_6 :
    Q.eq (wmcMarg 0 4 diamondEdges diamondPriors) (5,6) = true := by decide
theorem diamond_bp_a_is_7_10 :
    Q.eq (wmcMarg 0 5 unrollEdges unrollPriors) (7,10) = true := by decide

-- ===== general convergence of the monotone true-path closure =====
-- one synchronous step: a parent is activated by any active child.
def step (E : List (Nat × Nat)) (s : Nat → Bool) : Nat → Bool :=
  fun v => s v || E.any (fun e => decide (e.2 = v) && s e.1)

-- the closure is EXTENSIVE: it only ever adds.
theorem step_ext (E : List (Nat × Nat)) (s : Nat → Bool) (v : Nat) :
    s v = true → step E s v = true := by
  intro h; simp [step, h]

-- the closure is MONOTONE in the active set.
theorem step_mono (E : List (Nat × Nat)) (s s' : Nat → Bool)
    (hss : ∀ v, s v = true → s' v = true) :
    ∀ v, step E s v = true → step E s' v = true := by
  intro v h
  simp only [step, Bool.or_eq_true, List.any_eq_true] at h ⊢
  rcases h with h | ⟨e, he, hb⟩
  · exact Or.inl (hss v h)
  · refine Or.inr ⟨e, he, ?_⟩
    simp only [Bool.and_eq_true] at hb ⊢
    exact ⟨hb.1, hss e.1 hb.2⟩

-- iterate the closure
def iter (f : α → α) : Nat → α → α
  | 0,    x => x
  | n+1,  x => iter f n (f x)

-- CONVERGENCE: once a step adds nothing, all further iterations are fixed.
theorem fixed_stable (E : List (Nat × Nat)) (s : Nat → Bool)
    (hfix : step E s = s) : ∀ k, iter (step E) k s = s := by
  intro k
  induction k with
  | zero => rfl
  | succ k ih => simp only [iter, hfix]; exact ih

-- ===== axiom audit =====
#print axioms polytree_sink_exact
#print axioms diamond_sink_bp_wrong
#print axioms diamond_reconv_inexact
#print axioms diamond_exact_a_is_5_6
#print axioms step_mono
#print axioms fixed_stable
