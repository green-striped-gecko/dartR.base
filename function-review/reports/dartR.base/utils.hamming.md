# Review: utils.hamming (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.hamming. Reviewed in the
  distance wave (five files, one approval round, per-function PRs),
  under the standing member directive that utility functions are not
  for end users.
- Datasets: constructed 10x30 SNP fixture with missing data (plus its
  allele-relabelled mirror), testset.gs, raw-vector sequence fixtures,
  brute-force references.
- Family mode: analysis (distance/resampling kernel).
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: PASS (with notes)** — compact and documented.
- **Spec: FAIL — and the function is an ORPHAN.** No callers remain
  anywhere in the family (gl.filter.hamming and gl.report.hamming use
  utils.hamming.engine). The comparison starts AT the last base of
  the recognition sequence (substr from r; the docs say immediately
  after, i.e. r+1), and the proportion divides by the full minimum
  string length rather than the number of compared bases: one
  difference in four compared bases reports 0.125 (verified), where
  the documented computation gives 0.25.

## Findings

### D8 — Orphan disposition (MEDIUM) [member decision]

Recommended: REMOVE the file (unexported, uncalled, superseded by
the verified utils.hamming.engine; its two defects make it a trap
for future callers). Alternative: fix the off-by-one and the
denominator and keep it as a scalar reference implementation.

## Coverage

test-utils.hamming.R — 2 assertions snapshotting current behaviour
(dropped with the file if removal is approved).

## Approval

All findings approved via the approval boxes (2026-09-01); utils.hamming
disposition: REMOVE.

## Outcome

REMOVED per the approved disposition (file, Rd and its baseline test deleted). No callers anywhere in the family; the verified utils.hamming.engine is the live implementation. PR recorded below.

```json
{"function": "utils.hamming", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "PASS", "verdict_spec": "FAIL",
 "findings": [{"id": "D8", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.hamming.r disposition", "status": "applied"}],
 "datasets": ["constructed", "testset.gs"],
 "baseline_test": "tests/testthat/test-utils.hamming.R",
 "pr": null}
```
