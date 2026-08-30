# Review: gl.make.recode.pop (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.make.recode.pop
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT), pop-less variant
- Family mode: modify-adjacent utility (writes a proforma file)
- Checks skipped: FBM path not exercised (not available: no FBM
  fixture); Google Group not searched (not available: no browser session
  — GitHub issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — the indented `@family` tag leaks into the
  rendered help title ("…population names @family data manipulation",
  no \concept — verified in the Rd, the same defect class as PRs
  #256/#258/#260); visible NULL return; SET WORKING DIRECTORY before
  FLAG SCRIPT START; `@param outpath` copy-pasted from a plot function.
- **Spec: FAIL (docs contract only)** — `@return` promises "A vector
  containing the new population names" but the function returns NULL,
  visibly. Functionally exact: 30-row proforma matching the population
  levels, identity round-trip through gl.recode.pop, pop-less objects
  error clearly, SilicoDArT works, verbose = 0 silent, input untouched.

## Findings

### F1 — Visible NULL return contradicting the @return doc (LOW, confidence: certain)

Rule: spec axis (docs contract); VRB5. Location:
`R/gl.make.recode.pop.r:107` and header. Identical to
gl.make.recode.ind F1 (PR #269). Proposed change:
`return(invisible(NULL))` + corrected @return.

### F2 — @family leaks into the help title; header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.make.recode.pop.r:1-49`.

- Line 5 `#'  @family data manipulation` is indented after the title
  continuation — rendered title corrupted, no \concept (verified).
- `@param outpath` says "Directory to save the plot RDS files".
- @details says "recoding individuals" where populations are meant.
- DOC1 order; DOC2 wording; Custodian only; dead commented line; no
  @seealso to gl.recode.pop.

Proposed change: rewrite the header to the ratified template.

### F3 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.make.recode.pop.r:57-67`. Proposed change:
FLAG SCRIPT START directly after SET VERBOSITY.

## Report notes (other functions / not fixed here)

- With this, both proforma-writers match (PR #269 fixed the ind
  sibling).

## Coverage

Characterization baseline: `tests/testthat/test-gl.make.recode.pop.R` —
10 assertions: exact proforma vs population levels, verbose-0 silence,
identity round-trip through gl.recode.pop, pop-less error, SilicoDArT +
input untouched, visible NULL (baseline). All 10 pass on the pre-fix
code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all three findings
**approved**.

## Outcome

Applied F1–F3. Verification: all 10 characterization assertions pass
post-fix; the only diff from baseline is the approved invisible NULL
return. The @family concept now renders and gl.make.recode.pop joins the
"data manipulation" family index — the family's ~20 sibling Rd files
gain the cross-reference line (mechanical roxygen consequence, the same
ripple as PRs #256/#258). Proforma content and the identity round-trip
through gl.recode.pop verified exact. Caller grep across dartR.base + 7
siblings: tutorial bare calls at default verbosity only. dartr2shiny:
not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.make.recode.pop",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "LOW", "rules": ["spec", "VRB5"], "loc": "R/gl.make.recode.pop.r:107", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.make.recode.pop.r:1-49", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.make.recode.pop.r:57-67", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs", "pop-less variant"],
  "baseline_test": "tests/testthat/test-gl.make.recode.pop.R",
  "pr": null
}
```
