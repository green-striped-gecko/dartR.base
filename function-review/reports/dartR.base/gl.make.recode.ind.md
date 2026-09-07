# Review: gl.make.recode.ind (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.make.recode.ind
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT)
- Family mode: modify-adjacent utility (writes a proforma file; the
  object itself is untouched)
- Checks skipped: FBM path not exercised (not available: no FBM
  fixture); Google Group not searched (not available: no browser session
  — GitHub issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — visible NULL return (a bare call prints "NULL"),
  SET WORKING DIRECTORY before FLAG SCRIPT START, header breaks
  DOC1/DOC2/DOC7 with a copy-pasted `@param outpath` describing "plot RDS
  files".
- **Spec: FAIL (docs contract only)** — `@return` promises "A vector
  containing the new individual names" but the function returns NULL. The
  functional behaviour is fully correct: the proforma is written exactly
  (two identical columns matching indNames on both datatypes), round-trips
  through gl.recode.ind as an identity recode, verbose = 0 is silent, and
  the input object is untouched with no history append (appropriate — the
  object is not modified).

## Findings

### F1 — Visible NULL return contradicting the @return doc (LOW, confidence: certain)

Rule: spec axis (docs contract); VRB5. Location:
`R/gl.make.recode.ind.r:100` and header.

`return(NULL)` prints "NULL" on a bare call (verified) while the docs
promise a vector of names. Proposed change: `return(invisible(NULL))`
with the doc corrected to state that the proforma is written to file and
NULL is returned invisibly.

### F2 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.make.recode.ind.r:1-49`.

- `@param outpath` says "Directory to save the plot RDS files" — this
  function writes a recode csv, not plots (copy-paste from a plot
  function).
- `@details` after `@param`; `@return` after `@export` (DOC1); DOC2
  wording variant; `@author` Custodian only (DOC7).
- Dead commented-out progress line; no `@seealso` to gl.recode.ind
  (the docs mention it in prose only).

Proposed change: rewrite the header to the ratified template.

### F3 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.make.recode.ind.r:59-67`. SET WORKING
DIRECTORY precedes FLAG SCRIPT START. Proposed change: move FLAG SCRIPT
START to directly follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- gl.make.recode.pop is the sibling proforma-writer and likely carries
  the same return/doc mismatch — for its own review.
- write.table's default quoting round-trips cleanly through read.csv in
  gl.recode.ind (verified); left as-is.

## Coverage

Characterization baseline: `tests/testthat/test-gl.make.recode.ind.R` —
10 assertions: exact two-column proforma matching indNames, silence at
verbose = 0, round-trip identity recode through gl.recode.ind, SilicoDArT
path, input untouched, visible NULL return (baseline). All 10 pass on
the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all three findings
**approved**.

## Outcome

Applied F1–F3. Verification: all 10 characterization assertions pass
post-fix; the only diff from baseline is the approved invisible NULL
return (no more "NULL" auto-print). Proforma content and the identity
round-trip through gl.recode.ind verified exact on both datatypes.
Caller grep across dartR.base + 7 siblings: no references at all.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.make.recode.ind",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "LOW", "rules": ["spec", "VRB5"], "loc": "R/gl.make.recode.ind.r:100", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.make.recode.ind.r:1-49", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.make.recode.ind.r:59-67", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.make.recode.ind.R",
  "pr": 269
}
```
