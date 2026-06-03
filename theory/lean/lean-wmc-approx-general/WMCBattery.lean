/-
  WMCBattery.lean  (AMENDS lean-wmc-approx; does not modify it)
  Breadth evidence for "BP = exact WMC on the tree fragment, wrong at reconvergences":
  the same equalities/inequalities hold across MANY priors, each kernel-checked by `decide`.
  Standalone (Lean 4 core); axioms stay ⊆ {propext}.
-/
abbrev Q := Int × Int
def Q.add (a b : Q) : Q := (a.1 * b.2 + b.1 * a.2, a.2 * b.2)
def Q.mul (a b : Q) : Q := (a.1 * b.1, a.2 * b.2)
def Q.compl (a : Q) : Q := (a.2 - a.1, a.2)
def Q.eq (a b : Q) : Bool := decide (a.1 * b.2 = b.1 * a.2)
infixl:65 " +q " => Q.add
infixl:70 " *q " => Q.mul

def get (w : List Bool) (i : Nat) : Bool := (w[i]?).getD false
def worlds : Nat → List (List Bool)
  | 0 => [[]]
  | n+1 => (worlds n).flatMap (fun w => [false :: w, true :: w])
def valid (E : List (Nat × Nat)) (w : List Bool) : Bool :=
  E.all (fun e => (! get w e.1) || get w e.2)
def wWeight (pr : List Q) (w : List Bool) : Q :=
  (List.range pr.length).foldl
    (fun acc i => acc *q (if get w i then pr.getD i (0,1) else Q.compl (pr.getD i (0,1)))) (1,1)
def qsum (l : List Q) : Q := l.foldl (· +q ·) (0,1)
def wmcMarg (t n : Nat) (E : List (Nat × Nat)) (pr : List Q) : Q :=
  let v := (worlds n).filter (valid E)
  let num := qsum ((v.filter (fun w => get w t)).map (wWeight pr))
  let den := qsum (v.map (wWeight pr))
  (num.1 * den.2, num.2 * den.1)

-- single-edge BP update  b / (b + (1-b)(1-a))   [parent=index1, child=index0]
def bpEdge (a b : Q) : Q :=
  let den := b +q (Q.compl b *q Q.compl a)
  (b.1 * den.2, b.2 * den.1)

-- a spread of priors (num/den pairs in (0,1))
def PR : List Q := [(1,2),(1,3),(2,3),(1,4),(3,4),(2,5),(7,10),(1,5),(9,10)]

-- ===== SINGLE EDGE: BP update = exact WMC marginal, for ALL these priors =====
-- edge child(0) ⟹ parent(1); check parent marginal over a 9×9 grid of priors.
def edgeOK : Bool :=
  (PR.flatMap (fun a => PR.map (fun b =>
    Q.eq (wmcMarg 1 2 [(0,1)] [a,b]) (bpEdge a b)))).all (· = true)
theorem single_edge_exact_all_priors : edgeOK = true := by decide

-- ===== POLYTREE two-parent sink: BP independence form exact across priors =====
-- sink d(2) ⟹ b(0), d(2) ⟹ c(1); independent parents.
def bpSink (pd pb pc : Q) : Q :=
  let n := pd *q pb *q pc
  let den := n +q Q.compl pd
  (n.1 * den.2, n.2 * den.1)
def polyOK : Bool :=
  (PR.map (fun p => Q.eq (wmcMarg 2 3 [(2,0),(2,1)] [p,p,p]) (bpSink p p p))).all (· = true)
theorem polytree_exact_all_priors : polyOK = true := by decide

-- ===== DIAMOND: shared ancestor ⇒ BP independence form WRONG across priors =====
-- d(3)⟹b(1), d(3)⟹c(2), b(1)⟹a(0), c(2)⟹a(0); sink d shares ancestor a.
def diamondBad : Bool :=
  (PR.map (fun p => Q.eq (wmcMarg 3 4 [(3,1),(3,2),(1,0),(2,0)] [p,p,p,p]) (bpSink p p p))).all (· = false)
theorem diamond_wrong_all_priors : diamondBad = true := by decide

#print axioms single_edge_exact_all_priors
#print axioms polytree_exact_all_priors
#print axioms diamond_wrong_all_priors
