# Review: gl.keep.pop (dartR.base)

- Family mode: modify
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.5.0
- Package commit: 9344ef4
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.keep.pop.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the same three cosmetic/doc issues found on
`gl.drop.pop`/`gl.drop.ind` (missing `invisible()`, `@return` out of house
order, uncoloured verbose>=3 output), plus a DOC7 gap. Specifically checked
for the `gl.drop.ind` F1 pattern (`utils.reset.flags()` accidentally nested
inside a verbosity guard) — this function does **not** have it; the flag
reset runs unconditionally on the `recalc = FALSE` path.

**Spec: Ready** — population retention, `mono.rm`, `recalc`, `as.pop`, and
dual SNP/SilicoDArT support all behave as documented, verified against a
characterization test covering all of them.

## Findings

**F1 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.keep.pop.r:225` — same issue as the sibling `gl.drop.*` functions.
Confirmed: an unassigned call prints the full 22-line object summary to
console.
Proposed change: wrap the final `return(x)` in `invisible()`.

**F2 [MEDIUM, confidence: high] — `@return` out of house order (DOC1)**
`R/gl.keep.pop.r:47` — `@return` is the last roxygen tag, after `@export`,
same pattern as both `gl.drop.*` siblings.
Proposed change: move `@return` to directly follow the `@param` block.

**F3 [LOW, confidence: high] — uncoloured console output (VRB2)**
`R/gl.keep.pop.r:178-201` — the 20-line "Summary of recoded dataset" block
(both the `as.pop`-set and `as.pop`-null branches) uses plain
`cat()`/`cat(paste(...))` rather than `report()`. Confirmed silent at
`verbose = 0`.
Proposed change: wrap each `cat()` call in `report()`.

**F4 [LOW, confidence: high] — `@author` states only Custodian, no Author(s) (DOC7)**
`R/gl.keep.pop.r:29` — reads `Custodian: Arthur Georges -- Post to...`
with no separate `Author(s):` line, same gap as `gl.drop.ind` F6.
Proposed change per DOC7's remediation default: add
`Author(s): Arthur Georges.`

**F5 [LOW, confidence: medium] — `verbose` param text deviates from DOC2's canonical wording (proposed rule)**
`R/gl.keep.pop.r:25-27` — same wording (and same deviation) as every
sibling in this family. Same standing note: widespread across dartR.base,
worth an audit of DOC2's canonical text rather than a per-function fix.

## Cleared during verification (checked, not a defect)

- **The `gl.drop.ind` F1 pattern** (`utils.reset.flags()` nested inside
  `if (verbose >= 2)`) does **not** apply here — confirmed empirically that
  `verbose = 0` and `verbose = 2` produce identical (all-FALSE)
  locus-metric flags after a `recalc = FALSE` call. The code is written
  correctly (`x <- utils.reset.flags(...)` sits outside the verbosity
  guard); the misleading, unevenly-indented brace layout made this worth
  checking explicitly rather than trusting the visual indent.
- **DAT7**: no `accept=` override on `utils.check.datatype()`, but the
  function's docs and logic genuinely support both SNP and SilicoDArT
  (population-based subsetting, no dosage-specific math) — correct as
  unrestricted.
- **Tightened FS8**: single return path, history appended on the correct
  `x` immediately before it.
- **My own false-positive during verification**: the first `mono.rm=TRUE`
  test run showed a locus surviving monomorph filtering and I initially
  suspected a "filter keeps what it should remove" defect (the exact bug
  class `gl.filter.overshoot` had, per this skill's own fixtures). Traced
  it to my test's own flawed heuristic (`length(unique(col)) <= 1`), which
  wrongly classifies a locus where *every* individual is heterozygous
  (dosage = 1 for all) as monomorphic. That's actually maximally
  polymorphic -- both alleles are present in every individual. Only
  uniform dosage 0 or uniform dosage 2 is truly monomorphic. Fixed the
  test's heuristic, not the function; `gl.filter.monomorphs()` (called
  directly, independent of `gl.keep.pop`) was correct all along.
- Ploidy, individual/population counts, and history-append all confirmed
  correct for both SNP and SilicoDArT.
- No `accept =` restriction on `utils.check.datatype()` is correct and
  intentional, confirmed by both datatypes passing in the characterization
  test.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs, including a full
  internal-logic trace (population selection polarity, mono.rm/recalc
  interaction, as.pop reassignment/restoration) — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): not run — no
  API-affecting change proposed

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |

F5 carried no proposed fix (flagged for awareness only) and was not
offered as an approvable change.

## Outcome

- **F1 (invisible)**: applied at `R/gl.keep.pop.r:227`. Verified: unassigned
  call now prints 0 lines (was 22).
- **F2 (@return order)**: applied -- moved to directly follow `@param`.
- **F3 (colour the output)**: applied to both branches of the 20-line
  summary block.
- **F4 (Author(s) line)**: applied -- `@author` now reads
  `Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to...`.
- `devtools::document()` re-run; `man/gl.keep.pop.Rd`'s `\author{}` section
  updated to match (F4); no other section changed.
- Characterization test (`tests/testthat/test-gl.keep.pop.R`): 16/16 pass,
  no behavioural diff.
- Verified via a `verbose = 3` end-to-end run on both the `as.pop`-null
  (`mono.rm=TRUE, recalc=TRUE`) and `as.pop`-set branches -- both complete
  cleanly with coloured output.
- PR: not yet opened (pending pre-push confirmation).

```json
{
  "function": "gl.keep.pop",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.5.0",
  "commit": "9344ef4",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved"},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC2", "status": "not_offered"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": null
}
```
