# Review: gl.report.factorloadings (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.factorloadings
- Datasets: glPca object from gl.pcoa(testset.gl) — 111 loadings
- Family mode: report (input is a glPca object, not a genlight)
- Checks skipped: plot rendering verified only for object construction
  (headless run); Google Group not searched (not available: no browser
  session — GitHub issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — the report line and top-N table print
  unconditionally at verbose = 0 (VRB2); `@family matched reports`
  (plural) creates an orphan concept outside the "matched report" family
  (verified in the Rd); header breaks DOC1/DOC2/DOC7; the documented
  `...` parameter (claimed to reach ggsave) is never passed anywhere.
- **Spec: FAIL** — the `@return` doc says "The unchanged genlight object"
  but the input is a glPca and the return is an invisible data.frame of
  the axis loadings; `n.display` beyond the number of loci prints garbage
  NA rows (189 of them at n.display = 300), and `n.display = 0` prints
  one row anyway (the `1:0` slip). The top-N selection itself matches
  independent recomputation; plot.file works without display; the return
  is invisible.

## Findings

### F1 — Results print unconditionally at verbose = 0 (MEDIUM, confidence: certain)

Rule: VRB2 (ratified gate: verbose >= 1). Location:
`R/gl.report.factorloadings.r:112,120`.

The "Reporting factor loadings…" line and the top-N table are ungated
(verified: 17 lines at verbose = 0). Proposed change: gate at
`verbose >= 1`.

### F2 — n.display over/under-index (MEDIUM, confidence: certain)

Rule: spec axis (the gl.report.callrate ind.to.list class). Location:
`R/gl.report.factorloadings.r:120`.

`tmp[1:n.display, ]` over-indexes when n.display exceeds the number of
loadings (verified: 189 NA rows at n.display = 300 with 111 loadings) and
`1:0` prints one row when n.display = 0 (verified). Proposed change:
`head(tmp, n.display)` — prints exactly min(n.display, nrow) rows and
nothing for 0.

### F3 — Cryptic error for an out-of-range axis (LOW, confidence: certain)

Rule: spec axis (input validation). Location:
`R/gl.report.factorloadings.r:105`.

axis = 99 → "subscript out of bounds" (verified). Proposed change:
validate `axis <= ncol(pca$loadings)` with a clear fatal error.

### F4 — Header: wrong @return, plural @family, dead `...` (LOW, docs + plumbing, confidence: certain)

Rules: DOC1, DOC2, DOC5, DOC7. Location:
`R/gl.report.factorloadings.r:1-56,158-163`.

- `@return "The unchanged genlight object"` — the input is a glPca and
  the actual return is an invisible data.frame of the axis loadings.
  Corrected to describe the real return.
- `@family matched reports` (plural) — orphan concept; corrected to
  "matched report".
- `...` is documented as "passed to ggsave" but never forwarded;
  utils.plot.save accepts `...`, so it is now passed through to the save
  call, making the documented behaviour real.
- DOC1 order (`@return` after `@export`), DOC2 wording, DOC7 Custodian
  only.

Proposed change: rewrite the header; forward `...` to utils.plot.save.

### F5 — Datatype check style (LOW, confidence: certain)

Rule: STY. Location: `R/gl.report.factorloadings.r:94-100`.

`cat(error(...)); stop()` prints then raises an empty error; `class(pca)
!= "glPca"` breaks if a subclass ever appears. Proposed change:
`if (!inherits(pca, "glPca")) stop(error(...))`. The `tmp` variable is
also reused for both the top-N table and the save return — renamed in
passing.

## Report notes (other functions / not fixed here)

- gl.filter.factorloadings (next in the queue) carries the standing
  `loclist <- tmp$locus` stray-assignment note from earlier in the
  campaign.
- The glPca loadings carry only polymorphic loci retained by gl.pcoa
  (111 of testset.gl's 255) — expected, not a defect.

## Coverage

Characterization baseline:
`tests/testthat/test-gl.report.factorloadings.R` — 12 assertions:
invisible data.frame return sized to the loadings, top-N vs independent
recomputation, non-glPca rejection, plot.file save without display,
17-line verbose-0 output (baseline), NA-row garbage at n.display = 300 +
one-row print at n.display = 0 (baseline), cryptic axis error (baseline).
All 12 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all five findings
**approved**, including the escalation-gate consequence (F2: clamped
n.display output).

## Outcome

Applied F1–F5. Verification: all 13 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (verbose =
0 fully silent [F1]; no NA rows at n.display = 300, nothing printed at
n.display = 0 [F2]; clear axis error [F3]; `...` forwarded to
utils.plot.save [F4]; inherits-based type check [F5]). Moving out of the
plural family updates cross-references mechanically: the matched-report
siblings gain a factorloadings line, and gl.mahal.assign's "matched
reports" list loses one. **Report note:** the orphan plural family
"matched reports" still holds gl.report.bases, gl.report.fstat,
gl.report.monomorphs and gl.mahal.assign on dev — bases and monomorphs
are corrected in their unmerged PRs (#244/#249); fstat and mahal.assign
await their own reviews. Caller grep across dartR.base + 7 siblings:
only the matched filter's docs — no code caller. dartr2shiny: not
present in the workspace. NEWS entry added.

```json
{
  "function": "gl.report.factorloadings",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["VRB2"], "loc": "R/gl.report.factorloadings.r:112,120", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.report.factorloadings.r:120", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.report.factorloadings.r:105", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC5", "DOC7"], "loc": "R/gl.report.factorloadings.r:1-56,158-163", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.report.factorloadings.r:94-100", "status": "proposed"}
  ],
  "datasets": ["glPca from gl.pcoa(testset.gl)"],
  "baseline_test": "tests/testthat/test-gl.report.factorloadings.R",
  "pr": null
}
```
