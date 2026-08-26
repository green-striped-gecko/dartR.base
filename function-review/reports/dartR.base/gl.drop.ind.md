# Review: gl.drop.ind (dartR.base)

- Family mode: modify
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.4.0
- Package commit: a8b6a3c
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.drop.ind.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — one HIGH-severity data-integrity bug (locus
metric flags only reset at `verbose >= 2`) plus the same three cosmetic/doc
issues found in the sibling `gl.drop.pop` review.

**Spec: Ready** — individual removal, `mono.rm`, `recalc`, and dual
SNP/SilicoDArT support all behave as documented; the one behavioural defect
found is a Standards-axis bug (verbosity controlling data logic), not a
mismatch between docs and intent.

## Findings

**F1 [HIGH, confidence: high] — locus-metric flags only reset when `verbose >= 2` (DAT4)**
`R/gl.drop.ind.r:133-138` —
```r
} else {
    if (verbose >= 2) {
        cat(warn("  Locus metrics not recalculated\n"))
        x <- utils.reset.flags(x, verbose = 0)
    }
}
```
`x <- utils.reset.flags(x, verbose = 0)` is nested *inside* the
`if (verbose >= 2)` block instead of being a sibling statement after it
(compare the sibling `gl.drop.pop.r`, where the equivalent call sits outside
the verbosity check). Confirmed empirically: after dropping 2 individuals
with `recalc = FALSE` (the default), calling at `verbose = 0` leaves 10
locus metrics (`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`,
`CallRate`, `maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) flagged `TRUE`
(current/trustworthy) even though the underlying statistics are now stale;
the identical call at `verbose = 2` correctly resets all of them to `FALSE`.
Failure scenario: any quiet or scripted call (`verbose = 0` or `1`, a very
common default in pipelines) silently leaves every consumer of
`loc.metrics.flags` — other dartR functions or user code checking "is this
metric current?" — trusting statistics that no longer reflect the retained
individuals, with no warning that anything is wrong. Verbosity should only
ever control console output, never data-processing logic.
Proposed change: move `x <- utils.reset.flags(x, verbose = 0)` out of the
`if (verbose >= 2)` block so it always runs on the `recalc = FALSE` path,
matching `gl.drop.pop`'s structure. Note: this reintroduces the same minor,
already-accepted `gl.drop.pop` F5 pattern (the `monomorphs` flag reads
FALSE even right after `mono.rm = TRUE` just confirmed none remain) — not
treated as a separate finding here, since it's the correct, intentional
consequence of `utils.reset.flags()`'s blanket-invalidation design, and
was already flagged for awareness (not fixed) on the sibling function.

**F2 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.drop.ind.r:166` — same issue as `gl.drop.pop` F1. Confirmed: an
unassigned call prints the full 22-line object summary to console.
Proposed change: wrap the final `return(x)` in `invisible()`.

**F3 [MEDIUM, confidence: high] — `@return` out of house order (DOC1)**
`R/gl.drop.ind.r:44` — `@return` is the last roxygen tag, after `@export`,
same pattern as `gl.drop.pop` F2. No functional effect.
Proposed change: move `@return` to directly follow the `@param` block.

**F4 [LOW, confidence: high] — uncoloured console output (VRB2)**
`R/gl.drop.ind.r:99, 144-153` — the verbose>=3 individual-list line and the
10-line "Summary of recoded dataset" block use plain `cat()`/`cat(paste(...))`
rather than `report()`. Confirmed silent at `verbose = 0`; uncoloured at
`verbose = 3`. Same pattern as `gl.drop.pop` F3.
Proposed change: wrap each `cat()` call in `report()`.

**F5 [LOW, confidence: medium] — `verbose` param text deviates from DOC2's canonical wording (proposed rule)**
`R/gl.drop.ind.r:20-22` — identical wording (and identical deviation) to
`gl.drop.pop` F4. Same note applies: this phrasing is widespread across
dartR.base, worth an audit of DOC2's canonical text rather than a
per-function fix.

## Cleared during verification (checked, not a defect)

- **`x$ind.names` vs `indNames(x)`** (line 105 uses `$`, line 81 uses the
  accessor) — confirmed identical results empirically
  (`identical(indNames(x), x$ind.names)` is `TRUE`); not a bug.
- **The recurring `loc.metrics.flags$monomorphs == FALSE` crash** does not
  apply here either: line 109 unconditionally assigns a definite `FALSE`
  literal before the line 119 check.
- Ploidy, individual counts, and history-append all confirmed correct for
  both SNP and SilicoDArT.
- `recalc = TRUE` path is unaffected by F1 — `gl.recalc.metrics()` is called
  unconditionally regardless of verbosity, and correctly leaves all flags
  `TRUE` (accurately, since it just recalculated them).
- `for (case in ind.list)` mutating `ind.list` mid-iteration (lines 80-85)
  is safe in R, same reasoning as the sibling review.
- No `accept =` restriction on `utils.check.datatype()` is correct and
  intentional, confirmed by both datatypes passing in the characterization
  test.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): not run — no
  API-affecting change proposed

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| F1 | approved | Arthur | Escalation gate: consequence (loc.metrics.flags change for verbose=0/1 + recalc=FALSE) explicitly approved |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |

F5 carried no proposed fix (flagged for awareness only) and was not offered
as an approvable change.

## Outcome

- **F1 (flag reset at all verbosity levels)**: applied at
  `R/gl.drop.ind.r:133-139` — moved `x <- utils.reset.flags(x, verbose = 0)`
  out of the `if (verbose >= 2)` block. Verified: `verbose = 0` and
  `verbose = 2` now produce identical (all-FALSE) locus-metric flags after
  dropping individuals with `recalc = FALSE`. Escalation gate: grepped
  `dartR.captive`, `dartR.popgen`, `dartR.sim`, `dartR.spatial`,
  `dartR.sexlinked`, `dartRstartup`, `dartRverse` for `gl.drop.ind(` callers
  — found 8 call sites (dartR.captive: `utils.assignment.r`/`_2`/`_3`/`_4`,
  `gl.filter.parent.offspring.r` x2; dartR.popgen:
  `gl.assign.pa.r`, `gl.assign.on.genotype.r`; dartRstartup: a tutorial
  `.Rmd` chunk, `eval=FALSE`) -- none read `loc.metrics.flags` or any
  affected flag name after calling it, so nothing depends on the stale-flag
  bug. No local `dartr2shiny` checkout exists to check. `NEWS.md` created
  (didn't exist before) with an entry for this fix.
- **F2 (invisible)**: applied at `R/gl.drop.ind.r:168`. Verified: unassigned
  call now prints 0 lines (was 22).
- **F3 (@return order)**: applied -- moved to directly follow `@param`.
  `devtools::document()` re-run; `man/gl.drop.ind.Rd` unchanged, as expected.
- **F4 (colour the output)**: applied to the verbose>=3 individual-list line
  and the 10-line summary block. Verified via a `verbose = 3,
  mono.rm = TRUE, recalc = TRUE` end-to-end run -- completes cleanly.
- Characterization test (`tests/testthat/test-gl.drop.ind.R`): 14/14 pass.
  The one test asserting the pre-fix buggy behaviour was rewritten to
  assert the fix instead (an approved behaviour change, not an unexplained
  diff).
- PR: not yet opened (pending pre-push confirmation).

```json
{
  "function": "gl.drop.ind",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.4.0",
  "commit": "a8b6a3c",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "pending"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "pending"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "pending"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "pending"},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "pending"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "awaiting-approval",
  "pr": null
}
```
