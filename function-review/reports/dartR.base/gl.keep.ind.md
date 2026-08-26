# Review: gl.keep.ind (dartR.base)

- Family mode: modify
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.6.0
- Package commit: 9344ef4
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.keep.ind.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the same HIGH-severity flag-reset bug found in
`gl.drop.ind` (F1 there), present here too, plus the same cosmetic/doc
issues found across the whole `gl.drop.*`/`gl.keep.*` family.

**Spec: Ready** — individual retention, `mono.rm`, `recalc`, and dual
SNP/SilicoDArT support all behave as documented; the defect found is a
Standards-axis bug (verbosity controlling data logic), not a spec mismatch.

## Findings

**F1 [HIGH, confidence: high] — locus-metric flags only reset when `verbose >= 2` (DAT4)**
`R/gl.keep.ind.r:133-138` —
```r
} else {
    if (verbose >= 2) {
        cat(report("  Locus metrics not recalculated\n"))
        x <- utils.reset.flags(x, verbose = 0)
    }
}
```
Identical defect class to `gl.drop.ind` F1 (same file family, apparently
the same copy-paste origin): `x <- utils.reset.flags(x, verbose = 0)` sits
inside the `if (verbose >= 2)` block instead of as a sibling statement.
Confirmed empirically: after retaining 20 individuals with
`recalc = FALSE` (the default), `verbose = 0` leaves 10 locus metrics
(`AvgPIC`, `OneRatioRef`, `OneRatioSnp`, `PICRef`, `PICSnp`, `CallRate`,
`maf`, `FreqHets`, `FreqHomRef`, `FreqHomSnp`) flagged `TRUE` (stale);
`verbose = 2` correctly resets all of them to `FALSE`. Same failure
scenario as `gl.drop.ind` F1: any quiet/scripted call silently leaves
consumers of `loc.metrics.flags` trusting statistics that no longer
reflect the retained individuals.
Proposed change: move `x <- utils.reset.flags(x, verbose = 0)` out of the
`if (verbose >= 2)` block, matching the fix already applied to
`gl.drop.ind` (PR #238) and matching `gl.drop.pop`/`gl.keep.pop`'s
already-correct structure.

**F2 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.keep.ind.r:165` — same issue as every sibling in this family.
Confirmed: an unassigned call prints the full 22-line object summary.
Proposed change: wrap the final `return(x)` in `invisible()`.

**F3 [MEDIUM, confidence: high] — `@return` out of house order (DOC1)**
`R/gl.keep.ind.r:39` — `@return` is the last roxygen tag, after `@export`,
same pattern as every sibling.
Proposed change: move `@return` to directly follow the `@param` block.

**F4 [LOW, confidence: high] — uncoloured console output (VRB2)**
`R/gl.keep.ind.r:143-152` — the 9-line "Summary of recoded dataset" block
uses plain `cat()`/`cat(paste(...))` rather than `report()`. Confirmed
silent at `verbose = 0`.
Proposed change: wrap each `cat()` call in `report()`.

**F5 [LOW, confidence: high] — `@author` states only Custodian, no Author(s) (DOC7)**
`R/gl.keep.ind.r:24` — same gap as `gl.drop.ind` F6 / `gl.keep.pop` F4.
Proposed change per DOC7's remediation default: add
`Author(s): Arthur Georges.`

**F6 [LOW, confidence: medium] — `verbose` param text deviates from DOC2's canonical wording (proposed rule)**
`R/gl.keep.ind.r:20-22` — reads "2, progress log" here (a third variant
from the two seen in siblings). Same standing note: widespread across
dartR.base, worth an audit of DOC2's canonical text rather than a
per-function fix.

## Cleared during verification (checked, not a defect)

- **VRB5 (new this pass)**: confirmed `verbose = 0` is fully text-silent
  (0 lines via `capture.output()`). No plot bundle on this function, so the
  graphics half of VRB5 doesn't apply.
- **DAT7**: no `accept=` override, but the function's docs and logic
  genuinely support both SNP and SilicoDArT (name-based individual
  retention, no dosage-specific math) — correct as unrestricted.
- **Tightened FS8**: single return path, history appended on the correct
  `x` immediately before it.
- **Logic polarity**: `which(x$ind.names %in% ind.list)` is the correct
  positive-selection ("keep") polarity, mirroring `gl.keep.pop`'s
  `which(pop(x) %in% pop.list)` and inverse of `gl.drop.ind`'s
  `which(!x$ind.names %in% ind.list)`.
- Ploidy, individual counts, and history-append all confirmed correct for
  both SNP and SilicoDArT.
- Used the *corrected* monomorphic-detection heuristic from the
  `gl.keep.pop` review (uniform dosage 0 or 2 only, not uniform dosage 1)
  in this baseline's `mono.rm=TRUE` test — no false alarm this time.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs, including a full
  internal-logic trace — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): required for F1
  (escalation gate) — not yet run, pending approval

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| F1 | approved | Arthur | Escalation gate: consequence explicitly approved |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |

F6 carried no proposed fix (flagged for awareness only) and was not
offered as an approvable change.

## Outcome

- **F1 (flag reset at all verbosity levels)**: applied at
  `R/gl.keep.ind.r:133-138` — moved `x <- utils.reset.flags(x, verbose = 0)`
  out of the `if (verbose >= 2)` block. Verified: `verbose = 0` and
  `verbose = 2` now produce identical (all-FALSE) locus-metric flags.
  Escalation gate: grepped `dartR.captive`, `dartR.popgen`, `dartR.sim`,
  `dartR.spatial`, `dartR.sexlinked`, `dartRstartup`, `dartRverse` for
  `gl.keep.ind(` callers — found 7 real call sites (dartR.captive:
  `utils.assignment.r`/`_2`/`_3`/`_4`; dartR.popgen: `gl.nhybrids.r`,
  `gl.assign.pa.r`, `gl.assign.on.genotype.r`), all at `verbose = 0`, none
  reading `loc.metrics.flags` or any affected flag name afterward. No local
  `dartr2shiny` checkout exists. `NEWS.md` entry added.
- **F2 (invisible)**: applied at `R/gl.keep.ind.r:166`. Verified: unassigned
  call now prints 0 lines (was 22).
- **F3 (@return order)**: applied — moved to directly follow `@param`.
- **F4 (colour the output)**: applied to the 9-line summary block.
- **F5 (Author(s) line)**: applied per DOC7.
- `devtools::document()` re-run; `man/gl.keep.ind.Rd`'s `\author{}` section
  updated to match; no other section changed.
- Characterization test (`tests/testthat/test-gl.keep.ind.R`): 13/13 pass.
  The F1 test was rewritten to assert the fix instead of the pre-fix bug
  (an approved behaviour change, not an unexplained diff).
- Verified via a `verbose = 3, mono.rm = TRUE, recalc = TRUE` end-to-end
  run — completes cleanly with coloured output.
- PR: not yet opened (pending pre-push confirmation).

```json
{
  "function": "gl.keep.ind",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.6.0",
  "commit": "9344ef4",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT4", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "not_offered"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": null
}
```
