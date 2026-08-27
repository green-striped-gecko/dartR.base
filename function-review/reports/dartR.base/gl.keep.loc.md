# Review: gl.keep.loc (dartR.base)

- Family mode: modify
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.6.0
- Package commit: 9344ef4
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.keep.loc.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — two crash-on-documented-usage defects (the
no-arguments path references an unset variable; `first` without `last`
crashes despite the docs promising a default), a range-clamp check testing
the wrong variable, plus the family's usual cosmetic/doc issues.

**Spec: Needs work** — the `@param last` doc promises "[if not specified,
last locus in the dataset]" but the code never implements that default
(confirmed crash), and out-of-range `last` values are silently accepted
with correct-by-accident results instead of the documented clamping
behaviour the sibling warnings imply.

## Findings

**F1 [HIGH, confidence: high] — `first` without `last` crashes; documented default never implemented (DOC5 / Spec)**
`R/gl.keep.loc.r:139-144` — the docs promise `last` "[if not specified,
last locus in the dataset]", but no code ever assigns that default. With
`first = 5` and `last` omitted, `first:last` evaluates `5:NULL` and crashes
with "argument is of length zero". Confirmed empirically.
Failure scenario: any user relying on the documented default — e.g.
`gl.keep.loc(gl, first = 100)` to keep loci 100 through the end — gets an
opaque crash instead.
Proposed change: add `if (is.null(last)) last <- nLoc(x)` at the top of the
range handling (before the bounds checks, which also read `last`).

**F2 [HIGH, confidence: high] — no-arguments call crashes on unset `flag` (FS5)**
`R/gl.keep.loc.r:79-87` — when neither `loc.list` nor `first` is supplied,
the final `else` prints a warning but never assigns `flag`; the very next
statement reads `flag` and crashes with "object 'flag' not found".
Confirmed empirically.
Failure scenario: `gl.keep.loc(gl)` produces a confusing internal-variable
error instead of the clear parameter message the warning intended.
Proposed change: make the no-arguments condition fatal with a clear
message — `stop(error("Fatal Error: Need to specify either a range of loci
to keep, or specific loci to keep\n"))` — consistent with FS5 fail-fast and
with the siblings' "no X listed to keep!" fatal errors.

**F3 [MEDIUM, confidence: high] — range clamp tests `first` but clamps `last` (Spec)**
`R/gl.keep.loc.r:109-118` — the block that warns "Upper limit ... cannot be
greater than the number of loci" is conditioned on `first > nLoc(x)` but
assigns `last <- nLoc(x)`. Two confirmed consequences: (a) when only `last`
is out of bounds (e.g. `first = 1, last = nLoc+500`), no warning fires and
no clamp happens — the result is correct only by accident (R's
out-of-bounds indexing yields `NA`s, which `%in%` then ignores); (b) when
*both* are out of bounds (`first = nLoc+50, last = nLoc+200`), the clamp
fires on the wrong variable, producing `first = 305, last = 255`, the
swap-warning then reverses them, and the user gets a single arbitrary locus
with no error.
Proposed change: condition the clamp on `last > nLoc(x)` (clamping `last`),
and add the analogous `first > nLoc(x)` case as a fatal error since no
sensible interpretation exists.

**F4 [MEDIUM, confidence: high] — missing `invisible()` on return (FS10)**
`R/gl.keep.loc.r:178` — same as every sibling. Confirmed: unassigned call
prints the full 22-line object summary.
Proposed change: wrap the final `return(x2)` in `invisible()`.

**F5 [MEDIUM, confidence: high] — `@return` out of house order (DOC1)**
`R/gl.keep.loc.r:39` — same as every sibling.
Proposed change: move `@return` to directly follow the `@param` block.

**F6 [LOW, confidence: high] — uncoloured console output (VRB2) + stray `)` in two warnings**
`R/gl.keep.loc.r:161-165` — the 4-line summary block uses plain `cat()`.
Also lines 105 and 114: the warning strings end `"...set to 1\n)"` /
`"...\n)"` — a literal stray `)` printed after the newline.
Proposed change: wrap summary `cat()` calls in `report()`; move the stray
`)` out of both string literals.

**F7 [LOW, confidence: high] — `@author` states only Custodian, no Author(s) (DOC7)**
`R/gl.keep.loc.r:24` — same gap as siblings.
Proposed change: add `Author(s): Arthur Georges.`

**F8 [LOW, confidence: medium] — DAT4: no flag invalidation after locus removal (proposed-fix-free note)**
Removing loci does not invalidate per-locus metrics the way removing
individuals does (each surviving locus's metrics remain individually
correct), so the siblings' `utils.reset.flags()` call is legitimately
absent. However `loc.metrics.flags$monomorphs` arguably overstates
confidence after an arbitrary locus subset. Flagged for awareness only —
consistent with how the analogous F5 was handled on `gl.drop.pop`.

## Cleared during verification (checked, not a defect)

- **DAT2/DAT3**: `loc.metrics` re-subset from the original object with the
  same logical mask used for the genotypes — rows stay 1:1 with loci.
  Confirmed via the characterization test (locNames and metrics row counts
  agree).
- **DAT7**: unrestricted `accept=` is correct — name-based locus selection,
  no dosage-specific math, both datatypes pass.
- **Tightened FS8**: history appended on `x2` — the actual returned object
  — including on the empty-list path where `x2 <- x`. Correct.
- **VRB5**: `verbose = 0` fully text-silent (0 lines). No plot bundle.
- Ploidy and individual counts confirmed preserved for both datatypes.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on testset.gl / testset.gs, including
  edge-case probes of every range-parameter combination — run
- FBM path (DAT6): SKIPPED — no FBM-backed fixture readily available
- Sibling `dartR.*` / dartr2shiny caller grep (API1-3): not yet run —
  F1/F2/F3 change error/edge behaviour (crash -> default, crash -> clear
  error, silent-wrong -> clamped/error); required before merge if approved

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | Escalation gate: behaviour-change consequence explicitly approved |
| 2 | approved | Arthur | Escalation gate: behaviour-change consequence explicitly approved |
| 3 | approved | Arthur | Escalation gate: behaviour-change consequence explicitly approved |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |

F8 carried no proposed fix (flagged for awareness only) and was not
offered as an approvable change.

## Outcome

- **F1 (last default)**: applied — `if (is.null(last)) last <- nLoc(x)` at
  the top of the range block. Verified: `first = 250` with `last` omitted
  correctly retains loci 250-255 (6 loci); `first = 5` retains 251.
- **F2 (clear fatal error)**: applied — the no-arguments path now
  `stop(error(...))`s with the parameter message. Verified via
  characterization test. The guard fires only when both `loc.list` and
  `first` are NULL; the zero-length-`loc.list` path (which dartR.popgen's
  `gl.ld.haplotype` and dartR.sim's `gl.sim.WF.run` depend on) still
  warns and returns the object unchanged — covered by a new dedicated test.
- **F3 (range clamp)**: applied — clamp now tests and clamps `last`;
  out-of-range `first` is a fatal error. Verified: `last = nLoc+500` warns
  and clamps to nLoc; `first = nLoc+50` errors clearly.
- **F4 (invisible)**: applied. Verified: unassigned call prints 0 lines
  (was 22).
- **F5 (@return order)**: applied.
- **F6 (colour + stray parens)**: applied — summary block wrapped in
  `report()`; the stray `)` removed from both warning strings.
- **F7 (Author(s) line)**: applied per DOC7.
- Escalation gate for F1-F3: grepped all 7 sibling checkouts — 8 external
  call sites (dartR.popgen x7 across `gl.ld.haplotype.r`/`gl.select.panel.R`,
  dartR.sim x1 in `gl.sim.WF.run.r`), all `loc.list`-based, none pass
  `first`/`last` at all, none affected. No local dartr2shiny checkout
  (re-verified). `NEWS.md` entry added.
- `devtools::document()` re-run; `man/gl.keep.loc.Rd` regenerated.
- Characterization test (`tests/testthat/test-gl.keep.loc.R`): 17/17 pass.
  The three CHARACTERIZATION probes documenting the pre-fix crashes were
  rewritten to assert the fixed behaviour (approved behaviour changes, not
  unexplained diffs).
- `verbose = 3` end-to-end runs on all three selection modes (list, range
  with defaulted `last`, combined) complete cleanly with coloured output.
- PR: [green-striped-gecko/dartR.base#245](https://github.com/green-striped-gecko/dartR.base/pull/245)
  - open, based on `dev_arthur` (PR #238), same manifest-continuity
  reason as PRs #239/#240/#244.

```json
{
  "function": "gl.keep.loc",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.6.0",
  "commit": "9344ef4",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved"},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "approved"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved"},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved"},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "DAT4", "status": "not_offered"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available"],
  "status": "pr-open",
  "pr": 245
}
```
