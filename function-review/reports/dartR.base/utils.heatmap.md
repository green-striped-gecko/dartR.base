# Review: utils.heatmap (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.heatmap.
- Datasets: mtcars (the function is genlight-agnostic), testset.gl
  population distances via the caller gl.plot.heatmap;
  gplots::heatmap.2 as the reference implementation (deparse diff +
  behavioural comparison).
- Family mode: plotting utility (vendored fork).
- Author credit in file: Andy Liaw (original); Gentleman, Maechler,
  Huber, Warnes (revisions).
- Checks skipped: Google Group not searched (not available: no
  browser session).

## Verdicts

- **Standards: FAIL** — the file is a vendored FORK of
  gplots::heatmap.2 but nothing in the documentation says so or why;
  a full deparse diff shows the fork is deliberate and coherent:
  dendrogram plotting swapped to dendextend to support COLORED leaf
  labels (colRow/colCol), an auto-computed dendrogram margin,
  ColSideColors/RowSideColors defaulting to NULL instead of missing,
  a local invalid2() replacing gplots' unexported helper, and typo
  fixes ("dendogram"). All fork features verified working.
- **Spec: PASS (one crash)** — row/column reordering is IDENTICAL to
  gplots::heatmap.2 on the same input (verified), the scale/carpet
  return contract holds, side colors validate informatively, and the
  caller gl.plot.heatmap runs clean. The one defect: the added
  auto-margin takes max(nchar(colnames(x))) — a matrix WITHOUT
  dimnames sends -Inf into par(mar=) and crashes ("invalid value
  specified for graphical parameter"), where the gplots original
  handles unnamed matrices.

## Findings

### H1 — Unnamed matrices crash on the auto-margin (MEDIUM, confidence: certain)

Proposed change: fall back to the base margin when colnames are NULL
(margin_dendrogram <- 0.5), preserving the auto-sizing for named
input. All family callers pass named matrices, so behaviour there is
unchanged.

### H2 — Fork provenance documented; visibility (LOW, confidence: certain)

Proposed change: a @details paragraph stating the function is a
modified copy of gplots::heatmap.2 and listing the modifications
(colored dendrogram labels via dendextend, auto margins, NULL
side-color defaults); @keywords internal per the utility-function
policy (stays exported).

## Report notes (other functions / not fixed here)

- gplots remains in DESCRIPTION Imports but gplots:: is used NOWHERE
  in the package since this fork replaced heatmap.2 — the dependency
  could be dropped (a DESCRIPTION decision for the custodian, not
  made here).
- gl.plot.heatmap (the only caller) runs clean; noted for its own
  review.

## Coverage

`tests/testthat/test-utils.heatmap.R` — 10 assertions: reordering
identical to gplots (anchor), unnamed-matrix crash (baseline), fork
features (colored labels, side colors, validation), scale contract.
All pass pre-fix.

## Approval

Both findings approved via the approval boxes (2026-09-01).

## Outcome

Both findings applied. Suite: 10/10 pass; the unnamed-matrix flip maps
to H1; reordering anchor vs gplots unchanged; gl.plot.heatmap caller
smoke completes (its own 2-line verbose-0 leak noted for its review).
PR recorded below.

```json
{"function": "utils.heatmap", "package": "dartR.base", "family_mode": "plot",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "PASS",
 "findings": [
  {"id": "H1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.heatmap.r margin", "status": "applied"},
  {"id": "H2", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.heatmap.r provenance", "status": "applied"}],
 "datasets": ["mtcars", "testset.gl"],
 "baseline_test": "tests/testthat/test-utils.heatmap.R",
 "pr": null}
```
