# Review: gl.report.hamming, gl.filter.hamming, utils.hamming.engine (dartR.base)

Reviewed as one matched set: the report simulates the filter by running
the same compiled engine, so a defect in either can make the other wrong.

- Family mode: report (`gl.report.hamming`), modify (`gl.filter.hamming`),
  analysis (`utils.hamming.engine`)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 984ab05 (`dev_luis`, synced with `origin/dev` at def8e82)
- Datasets: platypus.gl (1000 loci, 81 individuals), testset.gl,
  testset.gs
- Baseline: tests/testthat/test-gl.filter.hamming.R (17 expectations),
  tests/testthat/test-gl.report.hamming.R (9 expectations), both captured
  before review; all pass

## Verdicts

**Standards: Needs work** — structure, guards and metadata handling
conform; argument validation is incomplete and differs between the two
functions.
**Spec: Needs work** — the core algorithm is correct (report counts equal
the filter's at every threshold tested; no two kept loci lie within the
threshold), but invalid arguments give silently wrong results, and the
report's main result is visible only as console text at `verbose >= 2`.

What works well: the Rcpp engine (already reviewed PASS/PASS, PR #318)
is exact against brute force, and the filter keeps `loc.metrics` in sync
with the loci and preserves ploidy for SNP and SilicoDArT data.

## Findings

**F1 [MEDIUM, confidence: high] — silent semantic trap (API1, proposed rule)**
`R/gl.filter.hamming.r:107-117` — no upper bound on `threshold`. When
`threshold >= min.length`, the engine treats every comparable pair as a
duplicate (`utils.hamming.blocks.r:87-91`) and keeps one comparable locus.
Failure scenario: `gl.filter.hamming(platypus.gl, threshold = 50)` removes
948 of 1000 loci with no warning. The report caps its own table at
`min.length - 1` (`gl.report.hamming.r:178`), so it never shows this.
Proposed change: both functions stop with an error when
`threshold >= min.length`.

**F2 [MEDIUM, confidence: high] — silent truncation (API1, proposed rule)**
`R/gl.filter.hamming.r:107-117`, `R/gl.report.hamming.r:166-175` — the
guards reject only proportions (0 < t < 1). Any other non-integer is
truncated by the Rcpp `int` conversion: `threshold = 2.9` behaves as 2
(10 loci removed on platypus.gl, not the 11 removed at 3). The report
highlights no bar because no simulated threshold equals 2.9. This is the
same defect class as the golden fixture fixed in a67aa07. `threshold = NA`
and `threshold = c(1, 3)` fail with base-R messages ("missing value where
TRUE/FALSE needed", "the condition has length > 1").
Failure scenario: a user who types `threshold = 2.5` believing "up to 2.5
bases" gets threshold 2 silently.
Proposed change: require `threshold` to be a single non-negative whole
number in both functions, with an informative error otherwise.

**F3 [MEDIUM, confidence: high] — precondition gap (FS5; VRB4, proposed rule)**
`R/gl.filter.hamming.r:95-117` — `rs` and `min.length` are not validated.
`rs = -3`, `min.length = 0`, or a `min.length` longer than every tag
leaves no locus comparable; the filter returns the object unchanged and
reports that only at `verbose >= 2`. `gl.report.hamming.r:152-160` checks
`rs`, but compares it with `min.length` (unrelated quantities), and its
message says "greater than zero" while accepting 0.
Failure scenario: `gl.filter.hamming(x, rs = -3, verbose = 1)` removes 0
loci and prints nothing about why.
Proposed change: in both functions, require `rs` to be a whole number
>= 0 and `min.length` a whole number >= 1, with one shared message; in the
filter, warn at `verbose >= 1` when fewer than two loci are comparable.

**F4 [MEDIUM, confidence: high] — result not returned (DOC5, proposed rule; report family)**
`R/gl.report.hamming.r:335-355` — the loci-removed table, the report's
main result, is printed only at `verbose >= 2` and never returned; the
function returns `x` unchanged. At `verbose = 0` or 1 the computation runs
and its result is discarded.
Failure scenario: `tab <- gl.report.hamming(x, verbose = 0)` gives the
genlight back; the table cannot be used in a script without parsing
console output. Reviewed siblings differ: `gl.report.hwe` and
`gl.report.heterozygosity` return their table; `gl.report.callrate`
returns `x`.
Proposed change: return the loci-removed table invisibly.
Callers: no sibling dartR package calls either function. dartr2shiny
(`shiny_fun/Fun_gl.report.hamming.R:162,175,191`) stores the result as a
report object (`Myreport.hamming`), not as the working dataset, so the
change does not break it.

**F5 [LOW, confidence: high] — duplicated logic (principle: one source of truth; STY1)**
`R/gl.filter.hamming.r:128-146` and `R/gl.report.hamming.r:188-242` hold
two copies of the sequence trimming, comparability test and worst-to-best
ordering. The report's claim of an "exact simulation" holds only while the
copies stay identical; validation already differs between them (F3).
Failure scenario: a future edit to one copy (for example, a change to how
ties in call rate are ordered) makes the report's counts disagree with the
filter without any test failing.
Proposed change: move the shared preparation into one internal helper in
`R/utils.hamming.blocks.r`, called by both functions, and keep the test
that compares report counts with filter results.

**F6 [LOW, confidence: high] — documentation gaps (DOC5, proposed rule; DOC2)**
- `gl.filter.hamming` does not say that several SNPs on the same tag
  (secondaries, same `CloneID`) are near-duplicates by construction and
  are removed: 9 of the 11 loci removed from platypus.gl by default are
  secondaries, which overlaps with `gl.filter.secondaries`.
- `gl.report.hamming` `@param x` says "SNP data"; SilicoDArT works.
- `gl.report.hamming` `@param verbose` wording differs from DOC2.
Failure scenario: a user runs `gl.filter.secondaries` and then
`gl.filter.hamming` and cannot account for the counts.
Proposed change: document the secondaries overlap and fix both
parameter texts.

**F7 [LOW, confidence: high] — structure and messaging (FS2/FS3, VRB2)**
`R/gl.report.hamming.r:116-135` sets the working directory and colours
before FLAG SCRIPT START and passes the outdated `build =` argument to
`utils.flag.start`; `gl.filter.hamming.r:89` keeps it as a comment. The
`verbose >= 3` summaries in both files use uncoloured `cat()`. The filter
body is indented inconsistently (two and four spaces).
Failure scenario: none for results; the start message prints after work
has begun, and summaries do not follow the house message style.
Proposed change: reorder to the FS sequence, drop `build =`, wrap summary
lines in `report()`, re-indent the filter body.

**F8 [INFO, confidence: high] — manifest duplicate**
`utils.hamming.engine` is defined in `R/utils.hamming.blocks.r`, which
was reviewed as `utils.hamming.blocks` (PASS/PASS, PR #318, 5 brute-force
tests). The manifest lists the function twice under both names.
Proposed change: mark the `utils.hamming.engine` row done with PR #318;
no code change.

## Proposed changes

1. Stop with an error when `threshold >= min.length`, in both functions
   (F1). **Consequence: calls that currently run and remove nearly all
   comparable loci now stop with an error.**
2. Require `threshold` to be a single non-negative whole number, in both
   functions (F2). **Consequence: `threshold = 2.9`, which currently runs
   as 2, now stops with an error.**
3. Validate `rs` and `min.length` in both functions with one shared
   message; the filter warns at `verbose >= 1` when fewer than two loci
   are comparable (F3). **Consequence: invalid `rs`/`min.length` values
   that currently return the object unchanged now stop with an error.**
4. Return the loci-removed table invisibly from `gl.report.hamming` (F4).
   **Consequence: the return value changes from the genlight object to a
   data frame; `x <- gl.report.hamming(x)` would overwrite the data.**
5. Move the shared preparation into one internal helper (F5).
6. Documentation: secondaries overlap, `@param x`, `@param verbose` (F6).
7. Structure and messaging housekeeping (F7).
8. Manifest: mark `utils.hamming.engine` done under PR #318 (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both
  functions.
- Spec: behaviour against roxygen on platypus.gl, testset.gl, testset.gs
  — run.
- Report–filter agreement: removal counts at thresholds 0, 3 and 10
  compared with actual filter runs — equal.
- Property check: minimum Hamming distance among kept comparable loci at
  threshold 3 = 4 — passes.
- DAT2 metadata sync and history: `loc.metrics` rows equal the original
  rows for the kept loci; one history entry added — passes.
- DAT1 ploidy: SNP 2 and SilicoDArT 1 preserved — passes.
- NA `TrimmedSequence`: loci are skipped, not an error — passes.
- Engine correctness: not re-run; covered by PR #318 tests.
- FBM path (DAT6): SKIPPED — not exercised; both functions use `glNA()`
  rather than densifying, per a67aa07.
- Candidate cap (5000): SKIPPED — no reference data reaches it; the
  warning path was read, not run.
- Google Group search: SKIPPED — not available: no browser session.
  GitHub issues: none open for either function.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (now errors) approved explicitly |
| 2 | approved | Luis | consequence (now errors) approved explicitly |
| 3 | approved | Luis | consequence (now errors) approved explicitly |
| 4 | approved | Luis | return-value change approved explicitly |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

| Change | Applied | Evidence |
|---|---|---|
| 1 | yes | `threshold = 50` now errors in both functions (tests: "invalid arguments stop with an error", "argument checks") |
| 2 | yes | `threshold` = 2.9, NA, c(1, 3), -1 now error; 0.2 keeps the migration message |
| 3 | yes | `rs = -3`, `min.length` = 0 or 2.5 now error; `min.length = 100` warns at `verbose = 1` and returns `x` |
| 4 | yes | returns a 11-row data frame on platypus.gl; `verbose = 1` prints only start and end |
| 5 | yes | new `R/utils.hamming.prepare.r` used by both functions. Deviation from the proposal: placed in its own file rather than `utils.hamming.blocks.r`, to keep one function per file (FS1) |
| 6 | yes | secondaries paragraph in `gl.filter.hamming`; `@param x`, `@param verbose`, `@return` in `gl.report.hamming` |
| 7 | yes | FS order restored, `build =` removed, summaries wrapped in `report()`, filter body re-indented |
| 8 | yes | manifest row `utils.hamming.engine` marked done, PR #318 |

Snapshot result: every baseline diff maps to an approved change. Unchanged
on platypus.gl: default filter removes 11 of 1000 loci; loci-removed table
0:10 = 8, 9, 10, 11, 12, 12, 15, 15, 17, 19, 20; report counts equal filter
results at thresholds 0, 3, 10; minimum Hamming distance among kept
comparable loci = 4 at threshold 3; SilicoDArT (testset.gs) 237 loci,
ploidy 1.

Tests: test-gl.filter.hamming.R 22/22, test-gl.report.hamming.R 16/16,
test-utils.hamming.blocks.R pass. `devtools::document()` run; both
functions run end to end at `verbose = 3`; examples run with
`run_donttest = TRUE`. Full R CMD check: not run locally (left to CI).

PR: #417

```json
{
  "function": ["gl.report.hamming", "gl.filter.hamming", "utils.hamming.engine"],
  "package": "dartR.base",
  "family": ["report", "modify", "analysis"],
  "skill_version": "2.0.0",
  "commit": "984ab05",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "API1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "API1", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "none", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["DAT6: FBM path not exercised", "candidate cap not reached by reference data", "Google Group: no browser session"],
  "status": "pr-open",
  "pr": 417
}
```
