# Review: gl.fixed.diff (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.fixed.diff
- Datasets: testset.gl (4-population subset; single-population split for
  the zero-difference case)
- Family mode: analysis (fixed-difference statistics; numerical
  verification against an independent recomputation)
- Checks skipped: test=TRUE simulation distribution not validated
  beyond symmetry/coercion (gl.fdsim is its own review target); Google
  Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — the tloc-with-test warning prints at verbose 0;
  a dead, unreachable tloc.hold block sits inside the pairwise loop
  (tloc has already been coerced to 0 before it) with its own ungated
  cat; matrix dimnames are reassigned on every pass of the double loop;
  header drift.
- **Spec: FAIL** — the documented `mono.rm` parameter ([default TRUE],
  "monomorphic loci ... removed before beginning computations") is
  never referenced in the function body, and the flag logic that stands
  in its place is inverted: gl.filter.monomorphs runs only when
  `loc.metrics.flags$monomorphs` is TRUE (i.e. when the data are
  already monomorph-free) and monomorphs are retained with a warning
  whenever the flag is FALSE. Verified: a 4-population subset with 208
  monomorphic loci (255 → 47 after filtering) runs on all 255, $nloc
  reports up to 242, and mono.rm = FALSE produces results identical to
  the default. The core arithmetic is otherwise sound: the fd matrix
  matches an independent recomputation from the genotype matrix
  exactly, and the p1/p2 positional locus indexing is safe (ftable
  blocks are identically ordered and match locNames).

## Findings

### F1 — mono.rm never used; monomorph logic inverted (HIGH, confidence: certain)

Rule: spec axis (documented parameter silently unused — the
gl.sim.crosses class) + DAT. Location: monomorph block.

Consequences of the current behaviour: monomorphs are effectively never
removed (when present, the flag is FALSE, which routes to
retain-and-warn), so `$nloc` counts and the `$pcfd` denominators are
inflated — on the test subset pcfd is understated roughly fivefold —
and `$gl` is returned with monomorphic loci included. Proposed change:
honor mono.rm — when TRUE (default) filter monomorphs unless the flag
already certifies none; when FALSE retain them with the gated warning.
Escalation-gate class: `$pcfd`, `$nloc`, and `$gl` change for datasets
containing monomorphs under the default; `mono.rm = FALSE` reproduces
the current numbers. Downstream: gl.dist.pop(method = "fixed-diff")
consumes $pcfd and its distances change accordingly; gl.fdsim and
dartR.popgen::gl.collapse consume raw $fd counts, which monomorphic
loci cannot contribute to — unaffected.

### F2 — Warning leak and dead tloc.hold block (MEDIUM, confidence: certain)

Rule: VRB/STY. Location: tloc/test check and pairwise loop.

The "false positives can only be simulated for tloc=0" warning prints
at verbose 0 (verified). Inside the pairwise loop a second
`if (tloc != 0)` block (ungated cat, tloc.hold assigned and never
restored) is unreachable because tloc was already coerced. Proposed
change: gate the warning at verbose >= 2; delete the dead block.

### F3 — verbose>=4 return-listing drift; loop-invariant dimnames (LOW, confidence: certain)

Rule: DOC/STY. Location: reporting block and pairwise loop.

The listing names the last element `$prob` (the list element is
`$pval`), describes `$sdfpos` with the copy-pasted expected-count text
instead of the standard deviation, and calls `$gl` the "input" object
(it is the filtered output). Matrix dimnames and diagonal zeros are
reassigned on every iteration of the double loop — loop-invariant;
moved outside once. No numerical change.

### F4 — Header and style conformance (LOW, docs-only + idioms, confidence: certain)

Rules: DOC1, DOC2, STY. `@return` after `@export`; the details text
promises a per-comparison warning for sample sizes below 5 that does
not exist (the actual behaviour is a global minimum-sample warning at
n < 10, verbose >= 3) — details aligned to behaviour; pb doc notes the
bar displays at verbose >= 2; `!(tloc == 0) & test == TRUE` →
`tloc != 0 && test`; `c(rep(0, nloci))` idiom.

## Report notes (other functions / not fixed here)

- utils.is.fixed documents "@return TRUE ... or FALSE" but returns
  1/0/NA integers — for its own review.
- gl.dist.pop(method = "fixed-diff") divides $pcfd by 100 for a
  distance matrix; with F1 its distances change numerically (they were
  computed over monomorph-inflated denominators). Recorded there when
  gl.dist.pop is reviewed.

## Coverage

Characterization baseline: `tests/testthat/test-gl.fixed.diff.R` — 24
assertions: the six pairwise fd counts verified against an independent
recomputation, symmetry, silence/visibility/input-untouched, monomorph
retention + inert mono.rm (baseline), tloc-test coercion + verbose-0
warning line (baseline), fd-class input round-trip, fatal guards (tloc
range, <2 populations), zero-difference split-population case. All 24
pass on the pre-fix code. The test=TRUE simulation is exercised with
reps = 2 for structure only.

## Approval

All four findings approved via the approval boxes (2026-08-30); F1
approved as recommended (default removes monomorphs; mono.rm = FALSE
preserves the previous numbers).

## Outcome

All four findings applied. Characterization suite: 25/25 pass; the
flipped assertions map to F1 (default nLoc 47, mono.rm = FALSE
reproduces 255/242, raw $fd identical either way) and F2 (no output at
verbose 0). The dimnames/diagonal hoist (F3) produced no unexplained
diffs — with test = FALSE the exp/sd/pval matrices now carry dimnames
(previously bare NA arrays), a labelling improvement within F3's
scope; their NA fill is unchanged. End-to-end at verbose 3, and 4/5
display blocks, clean — the corrected listing text and the "Removing
monomorphic loci" report show as intended; test = TRUE run symmetric.
Caller impact as disclosed in the approval box: gl.fdsim and
dartR.popgen::gl.collapse consume raw $fd (unchanged); gl.dist.pop
(method = "fixed-diff") now receives the corrected $pcfd. NEWS entry
added.

```json
{
  "function": "gl.fixed.diff",
  "package": "dartR.base",
  "family_mode": "analysis",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec", "DAT"], "loc": "R/gl.fixed.diff.r monomorph block", "status": "applied"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB", "STY"], "loc": "R/gl.fixed.diff.r tloc checks", "status": "applied"},
    {"id": "F3", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/gl.fixed.diff.r reporting/loop", "status": "applied"},
    {"id": "F4", "severity": "LOW", "rules": ["DOC1", "DOC2", "STY"], "loc": "R/gl.fixed.diff.r header", "status": "applied"}
  ],
  "datasets": ["testset.gl"],
  "baseline_test": "tests/testthat/test-gl.fixed.diff.R",
  "pr": 278
}
```
