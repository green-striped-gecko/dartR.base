# Review: gl.report.allelerich (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.report.allelerich; Datasets: testset.gl
  (4-population subset, all-NA loci filtered)
- Family mode: report/analysis (rarefied allelic richness; returns
  invisible list — input untouched verified, no history append,
  correct)
- Author: Ching Ching Lau
- Checks skipped: bootstrap CI distribution not validated beyond
  structure (see report note on the statistic mismatch); Google Group
  not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — no `verbose == 0 → plot.display <- FALSE`
  guard, so the plot renders at verbose 0 (verified by page output);
  the lazy signature default `plot.colors.pop = gl.colors("dis")`
  prints the 3-line gl.colors banner at every verbosity including 0
  (verified — the default-argument leak class recorded in the
  gl.colors review, PR #276, where changing the gl.colors default was
  declined; the ratified fix is the caller-side `verbose = 0`
  hand-patch); the package check calls requireNamespace on a length-3
  vector, which silently checks only dplyr (verified), and uses
  `cat(error()) + return(-1)`; boot and Rcpp (both Suggests, both
  required for nboots > 0) are never checked; dead code (an unused p1
  plot is built on every call, a commented parallel block, duplicated
  global-variable declarations); header drift.
- **Spec: FAIL (edges only)** — the core rarefaction is verified
  correct: per-population means match an independent recomputation of
  the El Mousadik & Petit formula exactly (global minimum n = 2 on the
  test subset). The failures are at the edges: an unrecognized
  `error.bar` value crashes downstream with "object 'max_val' not
  found" (verified), and `nboots > 0` silently overrides the user's
  `error.bar` choice to "CI" with no message (verified: SD requested,
  SD absent from the result).

## Findings

### F1 — Verbose-0 plot render and gl.colors banner leak (MEDIUM, confidence: certain)

Rule: VRB/PLT (ratified guard + the gl.colors lazy-default class).
Location: signature and print block.

Proposed change: add `if (verbose == 0) plot.display <- FALSE` after
verbosity is set; change the signature default to
`gl.colors("dis", verbose = 0)`.

### F2 — error.bar unvalidated; silent CI override (MEDIUM, confidence: certain)

Rule: spec/VRB. Location: error-bar blocks.

Unknown `error.bar` values crash with 'max_val' not found; with
nboots > 0 the user's choice is silently replaced by "CI". Proposed
change: validate error.bar upfront (unknown → warn at verbose >= 2 and
coerce to "SD", the house method-coercion pattern); report the
nboots-forces-CI switch at verbose >= 2.

### F3 — Package checks (MEDIUM, confidence: certain)

Rule: DEP/STY. Location: preliminaries and bootstrap block.

requireNamespace over a vector checks only the first package
(verified); `cat(error()) + return(-1)`; boot and Rcpp unchecked
before the bootstrap path. Proposed change: per-package loop with
`stop(error(...))`; check boot and Rcpp when nboots > 0.

### F4 — Dead code and minor validation (LOW, confidence: certain)

Rule: STY. The p1 plot is built on every call and discarded (the
p1/p2 patchwork line is commented out); commented parallel block;
duplicated global-variable NULL/NA declarations merged;
`cat(error()) + stop()` split in the nboots/CI check; boot.method
validated ("ind"/"loc" — a typo currently bootstraps silently as
"ind").

### F5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC2, DOC7. Verbose param canon wording; Custodian line added
to the author field; @return item names aligned with the actual list
element names.

## Report notes (other functions / not fixed here)

- utils.allelic.richness (the bootstrap statistic) rarefies to a
  per-population, per-replicate minimum sample size, whereas the point
  estimates rarefy to the single global minimum — the bootstrap CIs
  therefore target a slightly different statistic than the reported
  means. A methodological question for the authors; recorded, not
  changed.
- utils.allelic.richness compiles its C++ kernel inside the statistic
  via Rcpp::cppFunction (cached after the first call; the first
  bootstrap call carries seconds of compile time) — for its own
  review.
- p2's `coord_cartesian(ylim = c(0.8, max_val))` hardcodes the 0.8
  floor, which would clip bars for datasets with mean richness below
  0.8.

## Coverage

Characterization baseline:
`tests/testthat/test-gl.report.allelerich.R` — 14 assertions: the four
population means verified against an independent recomputation, result
names/popsize/default SD column, 3-line banner at verbose 0
(baseline), page rendered at verbose 0 (baseline), max_val crash
(baseline), bootstrap CI columns + silent SD override (baseline).
All 14 pass on the pre-fix code.

## Approval

All five findings approved via the approval boxes (2026-08-31).

## Outcome

All five findings applied. Characterization suite: 14/14 pass; flipped
assertions map to F1 (no banner lines, no page at verbose 0) and F2
(unknown error.bar coerces to SD instead of crashing); the bootstrap
CI structure and the nboots-forces-CI behaviour are unchanged (the
switch is now announced at verbose >= 2, verified). Richness values
unchanged (still match the independent recomputation). End-to-end at
verbose 3 clean; bootstrap path clean at verbose 2 with the switch
message displayed. NEWS entry added.

```json
{"function": "gl.report.allelerich", "package": "dartR.base", "family_mode": "report",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["VRB", "PLT"], "loc": "R/gl.report.allelerich.r signature/print", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "VRB"], "loc": "R/gl.report.allelerich.r error.bar", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["DEP", "STY"], "loc": "R/gl.report.allelerich.r pkg checks", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.allelerich.r", "status": "applied"},
  {"id": "F5", "severity": "LOW", "rules": ["DOC2", "DOC7"], "loc": "R/gl.report.allelerich.r header", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.report.allelerich.R", "pr": null}
```
