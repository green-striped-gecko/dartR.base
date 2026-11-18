# Review: gl.collapse (dartR.popgen)

- Family mode: modify (returns an `fd` object with a population-merged
  genlight and recomputed matrices)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 9ea1286 (origin/dev, reviewed state)
- Datasets: testset.gl (30 populations) through
  `gl.fixed.diff(tloc = 0.05)`, the roxygen example; 400 random symmetric
  fixed-difference matrices (4-9 populations) in synthetic `fd` objects,
  compared with connected components computed independently
- Baseline: tests/testthat/test-gl.collapse.R (6 tests, snapshot captured
  pre-review; defects marked `BASELINE (F<n>)`)

**Standards: Needs work** — structure follows the house order, but
arguments are not type-checked, output at verbosity 3 bypasses the crayon
helpers, level 4 is used, and the start flag uses the outdated `build`
argument.
**Spec: Needs work** — whenever the function completes, the groups are
right. But it fails in 61 of 400 random cases: some chains of similar
populations and the "all one group" case. The matrices it returns are
computed at its own `tloc` (default 0), not at the `tloc` the input was
built with.

What works well: in the 339 random cases that complete, the groups equal
the connected components of the "fd <= tpop" graph exactly, and the
testset.gl example collapses 30 populations to 3 as expected.

## Findings

**F1 [HIGH, confidence: high] — populations joined through a chain are
not always merged, and the function then stops (DOC5)**
`R/gl.collapse.r:137-167, 177-197` — the groups are built by one pass
over the pairs of candidate groups. When a chain of similar populations
is spread over groups that the single pass does not join, the same
population is left in two groups. The second `gl.merge.pop()` call then
stops: "Population(s) B, M, S, I, L, Q, Z not present in the dataset".
Failure scenario: 24 of the 363 random matrices that should give two or
more groups stop with this error, for example nine populations
L, S, M, Q, A, Z, I, R, B that form three groups at `tpop = 0`.
Proposed change: build the groups as connected components, repeating
the joining until no group changes (change 1). In the cases that work
today, the groups are the same.

**F2 [HIGH, confidence: high] — the returned matrices use gl.collapse's
`tloc`, not the one used to build `fd` (DOC5; API1, proposed rule)**
`R/gl.collapse.r:50-54, 200-203` — populations are grouped on
`fd$fd`, computed with whatever `tloc` the user gave `gl.fixed.diff()`.
The collapsed matrices are then recomputed with `gl.collapse(tloc = )`,
default 0. `gl.fixed.diff()` does not record its `tloc`, so the two can
differ without anyone noticing.
Failure scenario: the roxygen example builds `fd` with `tloc = 0.05` and
calls `gl.collapse(fd, tpop = 1)`. The returned `fd` equals
`gl.fixed.diff(collapsed, tloc = 0)`, not `tloc = 0.05`. A second
`gl.collapse()` on that result groups on the stricter definition.
Proposed change, one of:
- (a) check: recompute `fd$fd` from `fd$gl` at the given `tloc` and stop
  with a message naming the mismatch when it differs. This costs one extra
  `gl.fixed.diff()` run.
- (b) record `tloc` as an attribute in `gl.fixed.diff()` (dartR.base) and
  make `gl.collapse(tloc = NULL)` use it. This is a change in two packages;
  objects made before it carry no attribute and would still need (a) or the
  argument.
Recommended: (a) now, (b) as a dartR.base follow-up (change 2).
**Consequence: calls whose `tloc` does not match the one used to build
`fd` (including the roxygen example) stop with an error instead of
returning mixed-definition matrices; runtime roughly doubles.**

**F3 [MEDIUM, confidence: high] — when every population falls into one
group the function errors instead of returning it (DOC5)**
`R/gl.collapse.r:200-203` — `gl.fixed.diff()` needs at least two
populations: "Distance calculation requires at least two populations, one
or none present".
Failure scenario: 37 of 400 random matrices; for real data, any set of
populations with no fixed differences at the chosen `tpop`. That result
("these are all one diagnosable unit") is lost.
Proposed change: return the `fd` object with the single-population
genlight and 1 x 1 zero matrices, and report it at `verbose >= 1`
(change 3).
**Consequence: this case returns a result instead of an error.**

**F4 [MEDIUM, confidence: high] — argument and matrix checks missing
(FS5)**
`R/gl.collapse.r:70-81, 141` — `tpop` and `tloc` are compared with
numbers without checking type or length (`tpop = "1"` compares as text).
An NA in `fd$fd` (population pairs with no shared loci) stops the loop
with "missing value where TRUE/FALSE needed".
Proposed change: `tpop` one non-negative number, `tloc` one number in
[0, 0.5]; an NA distance does not join two populations, with a warning at
`verbose >= 1` (change 4).

**F5 [LOW, confidence: high] — when nothing merges, part of the returned
object comes from the input and part from the recomputation (DOC5)**
`R/gl.collapse.r:205-225` — `sdfpos` is taken from the input `fd`, the
other matrices from a recomputation without simulation. If the input
was made with `test = TRUE`, the returned `expfpos` and `pval` are NA
while `sdfpos` keeps the simulated values.
Proposed change: when no populations merge, return the input `fd`
unchanged (change 5).

**F6 [LOW, confidence: high] — verbosity and messages (FS3, VRB1, VRB2)**
`R/gl.collapse.r:60-62, 171-190, 228-236, 256-274`.
- `utils.flag.start()` is called with `build = "Jody"`.
- Level 4 is used for the fd matrix; it is reserved.
- The group listing at level 3 uses bare `cat()` and `print()`.
- Two messages contain hard-wrapped line breaks inside strings
  ("with \n tolerance").
- The closing summary calls `$gl` the "input genlight object" (it is the
  collapsed one) and lists `$prob` (the element is `pval`).
Proposed change: drop `build`; move level 4 to 5; use `report()`; fix
the strings (change 6).

**F7 [LOW, confidence: high] — roxygen gaps (DOC1, DOC5, DOC7 (proposed
rule))**
`R/gl.collapse.r:1-48`.
- The title has a typo ("less that") and says "less than a threshold",
  but the code joins populations at `<=`.
- `@return` lists `$expfpos` twice, omits `$sdfpos` and names `$prob`
  instead of `$pval`.
- There is no `@family`, and `@author` has no Author(s) part.
- `@details` does not say that groups are formed transitively: A and C
  are merged if A–B and B–C are each within `tpop`, even when A–C is
  not. Users need to know this before reading a merged group as one
  unit.
Proposed change: rewrite the header (change 7).

## Proposed changes

1. Groups as connected components (F1).
   **Consequence: cases that stopped with "not present in the dataset"
   now return; results that worked are unchanged.**
2. `tloc` consistency, option (a): check `fd$fd` against a recomputation
   at the given `tloc` and stop on mismatch; dartR.base follow-up (b)
   noted (F2).
   **Consequence: mismatched `tloc` (including the current example) now
   errors; runtime about doubles.**
3. Return the single-group result (F3).
   **Consequence: returns instead of erroring.**
4. Argument checks; NA distances do not join (F4).
5. No merge: return the input `fd` unchanged (F5).
6. Verbosity and message fixes (F6).
7. Roxygen rewrite, including the example's `tloc` (F7). Docs only.

Callers: dartr2shiny generator copy (`input_generator/dartR.popgen/
gl.collapse.r`); the unmaintained dartR.popgenomics repo holds an old
copy; no `dartR.*` sibling calls it.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. DAT: the genlight
  is changed only by `gl.merge.pop()` (dartR.base, reviewed separately).
- Spec: testset.gl example; 400 random matrices vs independent connected
  components — run.
- `tloc` mismatch: returned matrix vs recomputation at 0 and 0.05 — run.
- NA distances: from code reading; not reproduced with real data (needs
  population pairs without shared loci).
- dartR Google Group / GitHub issues: not searched.
- FBM path (DAT6): not applicable in this function (delegated to
  `gl.fixed.diff()` and `gl.merge.pop()`).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | option (a) chosen; consequence approved; (b) dartR.base follow-up |
| 3 | approved | Luis | consequence approved |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

- Changes 1-7 applied (branch `review-collapse`): `R/gl.collapse.r`;
  grouping moved to a new internal helper `R/utils.collapse.groups.r` (FS1)
  so it can be tested on hand-made matrices; `man/gl.collapse.Rd`
  regenerated; NEWS entry added. `@family fixed difference analysis` (the
  family name `gl.fixed.diff` uses in dartR.base).
- Snapshot diffs against the pre-review baseline: 6, all mapped. The
  testset example without `tloc` now stops on the mismatch (change 2, two
  tests); the synthetic `fd` objects used for the grouping tests no longer
  pass the `tloc` check (change 2, three tests; grouping is now tested
  through the helper); the single-group case returns (change 3).
- Tests rewritten for the approved behaviour: 10 tests, 122 expectations,
  all pass, including 100 random matrices against independent connected
  components, the chain case that used to fail, NA distances, the
  single-group and no-merge cases, and the `tloc` of the returned matrices.
- testset.gl at `verbose = 3` with `tloc = 0.05`, `tpop = 1`: three groups
  (240, 5, 5 individuals), as before.
- `devtools::check()`: 0 errors, 1 warning and 2 notes, all present before
  this change; example runs in 0.47 s.
- Follow-up (not in this PR): record `tloc` as an attribute in
  `gl.fixed.diff()` (dartR.base) so `gl.collapse` can default to it.
- PR: dartR.popgen#114 (commit aad1c71, branch `review-collapse`).

## Machine block

```json
{
  "function": "gl.collapse",
  "package": "dartR.popgen",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "9ea1286",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "API1", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["NA distances not reproduced on real data", "forum and GitHub issues not searched"],
  "status": "done",
  "pr": 114
}
```
