# Review: gl.report.contamination + gl.filter.contamination (dartR.base)

Matched pair reviewed together: the filter runs the report's screen and
removes what it flags, so each finding names the file it applies to.

## Provenance

- Family mode: report (`gl.report.contamination`), modify
  (`gl.filter.contamination`)
- Date: 2026-09-24
- Reviewer: Claude (Opus 5.5, claude-opus-5-5, Claude Code),
  dartr-function-review v2.0.0
- Package commit: 1dd6a3a (origin/dev), worktree branch
  `review-contamination`
- Datasets: testset.gl, testset.gs, platypus.gl, bandicoot.gl (dartR.data);
  a plain adegenet genlight built from testset.gl; an FBM conversion of
  testset.gl; a simulated 1500-individual, 2000-locus set on 16 plates
- Baseline: characterization blocks appended to
  `tests/testthat/test-gl.report.contamination.R` and
  `tests/testthat/test-gl.filter.contamination.R`; all tests pass at the
  reviewed state (report 64 expectations, filter 27)
- Release state: both functions are on `dev` only (not on `main`/CRAN), so
  no released caller depends on current behaviour. No callers in
  dartR.popgen, dartR.sim, dartR.captive, dartR.spatial, dartR.sexlinked
  or dartr2shiny.

## Verdicts

**Standards: Needs work** — structure, history, flag resets and plotting
follow the house pattern; parameter validation is inconsistent (two
parameters are silently replaced, the rest stop, and NA or vector input
crashes with base-R messages), and the filter hides every warning the
screen raises.
**Spec: Needs work** — behaviour matches the documentation on packaged and
simulated data, including FBM-backed and plain-genlight input; plate
positions from lower-case wells or plate names containing "-" silently
yield no adjacency, and the pair search slows quadratically with sample
count.

## Findings

**F1 [MEDIUM, confidence: high] — screen caveats hidden by the filter (VRB3, VRB4 proposed rule)**
`R/gl.filter.contamination.r:131-143` — the filter calls the report at
`verbose = 0`, so none of the report's warnings reach the user: a
replaced `rare.freq`/`min.n`, populations with one individual (never
screened, so never removed), and missing plate positions.
Failure scenario: `gl.filter.contamination(testset.gl, flag = "adjacent",
verbose = 3)` prints "Number of individuals removed: 0" with no hint that
`testset.gl` has no plate positions, so "adjacent" can never be assigned.
`gl.filter.contamination(testset.gl, rare.freq = 0.9, verbose = 3)` runs
with 0.02 and says nothing.
Proposed change: in the filter, warn at `verbose >= 2` when populations
with one individual exist, and at `verbose >= 1` when `flag` includes
"adjacent" but no plate positions were found (the filter's output is then
partial).

**F2 [MEDIUM, confidence: high] — invalid `rare.freq` and `min.n` are replaced, not rejected (FS5, API1 proposed rule)**
`R/gl.report.contamination.r:237-248` — out-of-range `rare.freq` becomes
0.02 and `min.n < 2` becomes 2, with the warning printed only at
`verbose >= 2`. Every other invalid parameter (`z.flag`, `min.excess`,
`min.share`, `share.tol`, `depth.*`, `plate`) stops.
Failure scenario: a user who reads `rare.freq` as "common above" passes
0.9; at `verbose = 1`, or through the filter at any verbosity, the screen
runs with 0.02 and the user never learns their value was discarded.
Proposed change: `stop(error(...))` for both, matching the other
parameters.

**F3 [LOW, confidence: high] — NA or vector parameters crash with base-R errors (FS5)**
`R/gl.report.contamination.r:237-259` — the checks use `if (a || b)` on
the raw arguments.
Failure scenario: `rare.freq = NA` gives "missing value where TRUE/FALSE
needed"; `z.flag = c(3, 4)` gives "'length = 2' in coercion to
'logical(1)'". Neither names the parameter at fault.
Proposed change: check each numeric parameter is a single finite number
before the range checks, and stop naming the parameter.

**F4 [LOW, confidence: high] — some plate positions silently yield no adjacency (spec: DOC5 proposed rule)**
`R/gl.report.contamination.r:382-397` — wells are matched against
upper-case `LETTERS`, and `plate_location` is split at every "-".
Failure scenario: a `plate` table with wells "c4" gives 0 adjacent pairs
where "C4" gives 1206 (testset.gl, 96-well layout). A `plate_location` of
"PL-1-C4" is read as plate "PL", well "1", so no well parses and adjacency
is empty, with no warning because the well column is not NA.
`gl.read.dart()` writes "<plate>-<row><col>", so any plate name containing
"-" triggers this.
Proposed change: upper-case wells before parsing; split `plate_location`
at the last "-" only; warn when positions exist but no well parses as
row letter + column number.

**F5 [LOW, confidence: medium] — adjacency pair search is quadratic in individuals (STY2)**
`R/gl.report.contamination.r:405-406` — every one of the n(n-1)/2 pairs
calls the R closure `is.adj()` through `apply()`.
Failure scenario: 1500 individuals on 16 plates, 2000 loci: 2.9 s without
plate positions, 6.3 s with them, and 62 % of the run is this `apply()`
(Rprof). The cost grows with n², to about 25 s extra at 4000 individuals,
plus a 16-million-row pair index.
Proposed change: build a plate+well key per individual and look up the
four neighbour wells directly (at most 4n lookups). Output is identical.

**F6 [LOW, confidence: high] — `plot.theme` documentation points at non-existent options (DOC5 proposed rule)**
`R/gl.report.contamination.r:138` — "See Details for options", but
Details lists none.
Failure scenario: a user looking for theme options finds nothing.
Proposed change: say "a ggplot2 theme, for example `theme_dartR()`".
Docs only.

Notes (not findings):
- Printing the suspect table at `verbose >= 1` departs from VRB1, but most
  `gl.report.*` functions print results at `verbose >= 1`, so this report
  flags the rule, not the code.
- The verbose `@param` text uses the newer "adopting the global verbosity"
  wording (69 files) instead of the DOC2 text (30 files); DOC2 is out of
  date, and this report flags the rule, not the code.
- An individual with every genotype missing gets `het = NaN`, is never
  flagged, but still receives a `top.partner` and `kin.z`. This is harmless
  noise.
- The filter matches removals by name, so duplicate `indNames` would drop
  every copy of a flagged name. Duplicate names are already invalid for
  most dartR functions.

## Proposed changes

1. `rare.freq` outside (0, 0.5) and `min.n < 2` stop with an error naming
   the parameter instead of being replaced (F2).
   **Consequence: calls that currently run with a silently replaced value
   now error.**
2. Validate that every numeric parameter is a single finite number, with
   an error naming it (F3).
3. The filter warns when populations with one individual were not screened
   and when `flag` includes "adjacent" but no plate positions exist (F1).
4. Parse plate positions robustly: upper-case wells, split
   `plate_location` at the last "-", warn when no well parses (F4).
   **Consequence: adjacency results change for data with lower-case wells
   or dashed plate names, from "no adjacency" to real adjacency, so some
   "suspect" individuals become "adjacent" and the filter's
   `flag = "adjacent"` output changes.**
5. Replace the all-pairs adjacency search with a neighbour-well lookup
   (F5). Output identical; faster.
6. Fix the `plot.theme` description (F6). Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both files
- Spec: documented tiers 1-3, flags and plate handling against behaviour on
  packaged data and the existing simulated tier-3 fixtures — run
- Report leaves input untouched (digest before and after) and appends no
  history (FS8) — run, passes
- Results independent of plotting (PLT3) — run, passes; `plot.display =
  FALSE` and `verbose = 0` return the full list
- FBM path (DAT6): run — FBM-backed testset.gl gives an identical `ind`
  table; the filter keeps `ind.metrics` in sync (249/249). The report
  densifies with `as.matrix()` as documented; peak memory was not measured
  at scale
- Plain adegenet genlight input — run, both functions work
- Independent numerical check: heterozygosity against `rowMeans(m == 1)`
  (existing test); kinship formula checked by reading (VanRaden GRM / 2) —
  not recomputed independently
- Downstream callers (API3): grep of sibling packages and dartr2shiny —
  none
- Google Group / GitHub issues: SKIPPED — functions are unreleased (dev
  only, added September 2026), so no user reports can exist

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (invalid values now error) approved explicitly |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | consequence (adjacency results change for such data) approved explicitly |
| 5 | approved | Luis | |
| 6 | approved | Luis | |

## Outcome

- Change 1 (F2): `rare.freq`/`min.n` stop. Test "bad parameters stop"
  replaces "warn and coerce"; the only snapshot expectation flipped, and it
  maps to this approved change.
- Change 2 (F3): scalar checks. `NA`, a vector and a character value each
  stop naming the parameter (same test).
- Change 3 (F1): filter warnings. New test "the filter repeats the screen's
  caveats"; `verbose = 1` warns only when "adjacent" is the sole class.
- Change 4 (F4): plate parsing. New tests: lower-case wells and
  "PL-1-<well>" give the same pairs as upper-case plain wells (1206 pairs,
  previously 0); unparsable wells warn.
- Change 5 (F5): neighbour lookup. New test compares pairs with a
  brute-force all-pairs search (two plates, a shared well). On 1500
  individuals, `ind` and `pairs` are identical to the reviewed code;
  8.5 s to 2.4 s.
- Change 6 (F6): `plot.theme` text; `man/` regenerated.
- Baseline: packaged-data snapshots unchanged (testset.gl, platypus.gl,
  bandicoot.gl outputs identical to the reviewed code). Report 74
  expectations, filter 32, all pass. Both functions run end to end at
  `verbose = 3`.
- `devtools::document()` also rewrote 14 unrelated Rd files (cross-link
  drift covered by PR #421); reverted, not part of this change.
- PR: pending.

```json
{
  "function": "gl.report.contamination",
  "companion": "gl.filter.contamination",
  "package": "dartR.base",
  "family": "report+modify",
  "skill_version": "2.0.0",
  "commit": "1dd6a3a",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 3},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "STY2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["Google Group/issues: functions unreleased", "kinship not recomputed independently", "peak memory at scale not measured"],
  "status": "awaiting-approval",
  "pr": null
}
```
