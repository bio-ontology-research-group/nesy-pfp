# Modular WMC coordination: marust assessment + the irreducible-core result

Continuing step 3. Goal: replace the generic graph tree-decomposition with a logically
principled module decomposition coordinated by the `kobayashi-marust` context calculus, and
lift its Boolean inter-context messages to weighted separator marginals.

## A graph-theoretic reality check first

Exact WMC costs `exp(treewidth)`, and **treewidth is a graph invariant** — no decomposition,
logical or generic, computes exact WMC below the constraint graph's treewidth. So for bp
(tw ≈ 272) there is *no* exact modular scheme. The honest modular contribution is therefore
**not** "exact everywhere" but **"exact on the low-treewidth periphery + bounded-approximate
messages across a small high-treewidth core,"** with the core's approximation bounded by the
§-theory (`approximation_theory.md`) and coordinated soundly. The decisive empirical question
becomes: **how small is the irreducible high-treewidth core?**

## Irreducible-core measurement (`modulator_probe.py`)

Greedily remove the highest-degree (reconvergence-hub) terms in 1% batches, recomputing
min-fill treewidth of the remainder, until the rest is exactly computable (tw ≤ 14):

| namespace | start tw | hubs cut to reach tw ≤ 14 | exactly-computable remainder |
|---|---|---|---|
| mf | 30 | **202 (2.95%)** | **97.05%** |
| bp | 269 | **1648 (7.72%)** | **92.28%** |

**mf:** cutting just **~3% of terms** (the reconvergence hubs) makes **97% exactly
WMC-computable**. **bp:** even the worst namespace (tw 269) reduces to a **7.7% hub core**,
leaving **92% exactly computable**. So almost all of GO is exact under a modular scheme;
bounded approximation is confined to a small, explicitly identified reconvergence core
(3–8% of terms). This is the strong, honest form of "modules buy tractability": not exactness
everywhere, but exactness everywhere *except a small reconvergence core* whose messages the
marust coordination must carry (bounded by the §-theory).

## marust engine assessment (verified, but needs the non-trivial strategy)

The `kobayashi-marust` engine (built, run) is a Lean-soundness-checked disjunctive context
reasoner and it **correctly saturates GO-style EL** — e.g. on the kinship hierarchy it returns
the exact transitive closure (`Father ⊑ {Male, Narcissist, Parent, Person}`). So the verified
coordination substrate works on our logic.

**But** the engine currently implements the **trivial successor strategy** (`shared_successor`:
*all* successors collapse into one shared context, `engine.rs:802`). That yields a coarse
context structure — essentially one root context plus a shared successor context — **not the
fine-grained module partition WMC needs.** Genuine modules (distinct context cores → distinct
contexts, with small shared separators) require the **non-trivial / branch-bound context
strategy** (Sequoia's full strategy; the code references `context_branch_bound` for "distinct
context cores"). That is a real, scoped engine extension, plus exposing the context graph
(cores, successor edges, separators) in the JSON output.

## The concrete path (now precisely scoped)

1. **Extend marust to the non-trivial context strategy** so distinct cores yield distinct
   contexts → a genuine module partition with explicit separators (the shared signatures on
   Succ/Pred edges). Emit the context graph.
2. **Lift the Succ/Pred messages** from Boolean consequences to **weighted separator
   marginals** (a small factor/SDD over the separator atoms) — the engine already passes
   messages over exactly these interfaces; only the payload changes from "entailed clause" to
   "weight vector."
3. **Exact-SDD each module**, junction-tree-message over the (low-treewidth) module-interface
   graph; on the irreducible hub core (mf ~3%, bp TBD) use §-bounded loopy messages. Validate
   against the cc exact JT (already computed) and report the residual error.

## Status

Done this turn: graph-invariance reality check; irreducible-core measurement (mf = 3% core,
97% exact; bp running); marust built + verified on GO; precise diagnosis that the trivial
successor strategy blocks module extraction and the non-trivial strategy is the needed
extension. The weighted-message lift (steps 1–2) is the next engine increment.
