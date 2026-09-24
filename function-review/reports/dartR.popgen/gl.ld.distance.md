# Review: gl.ld.distance (dartR.popgen)

- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 383279d (origin/dev, reviewed state)
- Datasets: platypus.gl (dartR.data) with the roxygen example set-up
  (complete loci, NCBIv1 chromosome and position), passed through
  `gl.report.ld.map(ld.max.pairwise = 1e7)`: 496 pairs, 3 populations;
  a two-population split of the same data; a synthetic 12-pair ld.report
  in the test
- Baseline: tests/testthat/test-gl.ld.distance.R (5 tests, snapshot
  captured pre-review; defects marked `BASELINE (F<n>)`)

**Standards: Needs work** — the structure is close to the house order, but
the result table prints at every verbosity, `fields` (Suggests) is used
without a guard and the start flag uses the outdated `build` argument.
**Spec: Needs work** — the binned means are correct, but some valid inputs
fail with unrelated errors, the documented palette option fails, and the
help page misdescribes the default and the returned columns.

What works well: the per-bin mean LD matches an independent `cut()` +
`mean()` computation on the platypus data to 5e-16 for every population
and bin.

## Findings

**F1 [MEDIUM, confidence: high] — some resolutions make the function fail
or add an empty bin (FS5, DOC5)**
`R/gl.ld.distance.r:78-96` — the breaks are
`c(seq(1, max(distance), ld.resolution), max(distance))`.
- When `ld.resolution` is at least the largest distance there is one bin.
  `stats.bin()$stats` is then a one-column matrix and
  `$stats[2, ]` reads past its end: "subscript out of bounds".
- When the sequence already ends on the maximum, that break appears twice
  and an empty, zero-width bin is added; its NA row shows as a gap and a
  warning in the plot.
- `ld.report` and `ld.resolution` are not checked.
Failure scenario: `gl.ld.distance(ld_res, ld.resolution = 1e7)` on the
roxygen example (largest distance 9,992,140) errors with "subscript out of
bounds". With `ld.resolution = max - 1` the output gains an NA row per
population.
Proposed change: check that `ld.report` has `pop`, `distance` and
`ld.stat`, and that `ld.resolution` is one positive number; build unique
breaks; read the bin statistics with `drop = FALSE` (change 1).

**F2 [MEDIUM, confidence: high] — verbosity not honoured, dependency not
guarded (VRB1, FS3, DEP1)**
`R/gl.ld.distance.r:70-74, 87, 146`.
- `print(bins_ld)` runs at every verbosity: 34 lines at `verbose = 0` for
  the roxygen example.
- `utils.flag.start()` is called with the outdated `build = "Jody"`.
- `fields::stats.bin()` is called without the `requireNamespace()` guard,
  although `fields` is in Suggests.
Failure scenario: a user without `fields` gets "there is no package called
'fields'" from inside the function; a quiet run still prints the table.
Proposed change: DEP1 guard for `fields`; print the table at
`verbose >= 3`; drop `build` (change 2).

**F3 [LOW, confidence: high] — `pop.colors` rejects the palette functions
its documentation offers (DOC5)**
`R/gl.ld.distance.r:14-16, 112-118, 136` — the value goes straight to
`scale_color_manual(values = )`. A palette function such as `rainbow`
errors; a vector shorter than the number of populations errors with
"Insufficient values in manual scale. 3 needed but only 1 provided."
Failure scenario: `pop.colors = rainbow`, as the help page suggests, fails.
Proposed change: call a function with the number of populations; stop
with a clear `error()` when a vector is shorter than that (change 3).

**F4 [LOW, confidence: high] — the red threshold line is unlabelled and
assumes r-squared (PLT1, DOC5)**
`R/gl.ld.distance.r:125-132` — `colour = "LD threshold for unlinked loci"`
is set inside `aes()` and then overridden by `color = "red"`, so the line
never reaches the legend. The line is drawn at 0.2, the r-squared
threshold cited in `@description`. `gl.report.ld.map()` can also return
D', LLR, OR, Q, covariance or R, and the report does not record which, so
for those the line has no meaning.
Failure scenario: a report made with `ld.stat = "D.prime"` shows a red
line at 0.2 labelled nowhere, which reads as a threshold for D'.
Proposed change: show the line in the legend as "r2 = 0.2 (unlinked
threshold)" and state in `@details` that it applies to R.squared
(change 4).

**F5 [LOW, confidence: medium] — bin means carry no pair counts (DOC5;
principle: results should allow the user to judge them)**
`R/gl.ld.distance.r:89-96` — the returned table has the mean LD per bin
but not how many pairs it rests on. In the platypus example bins hold 8 to
29 pairs per population, and every point is drawn the same way.
Failure scenario: a bin mean from 2 pairs and one from 300 look equally
reliable in the table and the plot.
Proposed change: add an `n.pairs` column (from `stats.bin()$stats["N", ]`)
(change 5).
**Consequence: the returned table gains a column; existing columns are
unchanged.**

**F6 [LOW, confidence: high] — roxygen gaps (DOC1, DOC5, DOC6 (proposed
rule), DOC7 (proposed rule))**
`R/gl.ld.distance.r:1-51`.
- The title is garbled ("…by population disequilibrium patterns").
- `ld.resolution` is documented as `[default NULL]` but defaults to 100,000.
- `plot.dir` says "working directory", but `gl.check.wd()` defaults to
  `tempdir()`.
- `@return` says "dataframe" (it is a data.table) and does not say that
  `distance` is the upper edge of each bin, with bins
  (1, res], (res + 1, 2 res + 1], ….
- `@author` has no Author(s) part.
- A reference contains a non-ASCII character ("André").
Failure scenario: a reader assumes each point sits at the start or
midpoint of its bin, or that the plot is saved in the working directory.
Proposed change: rewrite the header (change 6).

**A1 [addendum, Phase C] — `fields::stats.bin()` fails on a single bin**
Applying change 1 showed that the "subscript out of bounds" error is
raised inside `fields::stats.bin()` whenever `breaks` has two values, so
the single-bin case cannot be fixed while calling it. Proposed and
approved: bin with base R (`cut()` on the same bins, first bin closed;
mean and count per bin). The means are identical, and the function no
longer needs `fields`, which replaces the `fields` guard in change 2.

## Proposed changes

1. Input checks (columns of `ld.report`, one positive `ld.resolution`);
   unique breaks; single-bin case handled (F1). Results for inputs that
   work today are unchanged.
2. DEP1 guard for `fields`; table printed at `verbose >= 3`; `build`
   argument dropped (F2).
3. `pop.colors` accepts a palette function or a long-enough vector, with
   a clear error otherwise (F3).
4. Threshold line shown in the legend; `@details` says it applies to
   R.squared (F4). Plot only.
5. Add `n.pairs` to the returned table (F5).
   **Consequence: the returned table gains a column.**
6. Roxygen rewrite (F6). Docs only.

No change alters the signature. Callers: the dartr2shiny generator copy
(`input_generator/dartR.popgen/gl.ld.distance.r`, takes the roxygen
header); no `dartR.*` sibling calls the function.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run. DAT and FS8 not
  applicable (input is a data frame, output a new table).
- Spec: behaviour vs roxygen on the platypus example — run.
- Numerical check: binned means vs independent `cut()` + `mean()` — run
  (max difference 4.7e-16, all 30 bins).
- Plot: built and inspected via `ggplot_build()` and the saved RDS; not
  viewed — the visual check was not done.
- dartR Google Group: dartrverse Gmail searched for "gl.ld.distance", no
  messages. GitHub issues not searched.
- Other `ld.stat` values from `gl.report.ld.map()` (D.prime etc.): not
  run; the function treats `ld.stat` as a number whatever it measures.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis | consequence approved: returned table gains n.pairs |
| 6 | approved | Luis |  |
| A1 | approved | Luis | base R binning replaces fields::stats.bin; fields guard no longer needed |

## Outcome

- Changes 1-6 and addendum A1 applied in `R/gl.ld.distance.r` (branch
  `review-ld-distance`); `man/gl.ld.distance.Rd` regenerated; NEWS entry
  added. `fields` stays in DESCRIPTION Suggests although no dartR.popgen
  function now calls it (removal not in scope).
- Snapshot diffs against the pre-review baseline: 4, all mapped. Table no
  longer printed at `verbose = 0` (change 2); resolution above the largest
  distance returns one bin instead of failing (change 1 + A1); palette
  function accepted and short vector gives the new message (change 3, two
  expectations). The binned-means and plot tests passed unchanged.
- Tests rewritten for the approved behaviour: 10 tests, 27 expectations,
  all pass, including equality with `fields::stats.bin()` means and counts
  where that function works, empty bins, the single-bin and repeated-break
  cases, and the legend entry.
- platypus example at `verbose = 3`: the 30 bin means are identical to the
  pre-fix output; `n.pairs` 8-35 per bin; `ld.resolution = 1e7` returns one
  bin per population.
- `devtools::check()`: 0 errors, 1 warning and 2 notes, all present before
  this change.
- PR: dartR.popgen#112 (commit 09f0c7b, branch `review-ld-distance`).

## Machine block

```json
{
  "function": "gl.ld.distance",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "383279d",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["plot visual check", "GitHub issues not searched", "non-R.squared ld.stat not run"],
  "status": "pr-open",
  "pr": 112
}
```
