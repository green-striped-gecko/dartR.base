# Review: gl.run.popcluster + gl.plot.popcluster + gl.map.popcluster (dartR.popgen)
- Family mode: analysis (PopCluster wrapper, bar plot, map)
- Scope: one review, one report, one PR for the three functions (agreed with
  Luis, 2026-09-24); three manifest rows
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 7ddbca9 (origin/dev)
- Datasets: testset.gl, three populations (31 individuals, 18 loci);
  PopCluster 1.5.0.0 (build 20240925, `~/programs/PopClusterMac`); a
  hand-built run object for the plot and map
- Baseline: tests/testthat/test-gl.popcluster.R (run tests need
  `POPCLUSTER_DIR`; plot and map tests run anywhere)

## Verdict

**Standards: Needs work** — `gl.run.popcluster` writes its input files to
the working directory by default, never cleans up (the code is after
`return()`), ignores `plot_theme`, prints PopCluster's output at every
verbosity and checks arguments with `cat()` + empty `stop()`;
`gl.plot.popcluster` defaults `verbose` to 2 and causes the package's
standing R CMD check NOTE (bare column names in `aes()`).
**Spec: Needs work** — the likelihood table is kept as text, so the
LogL(K) plot is drawn upside down (the best K lowest) and the ΔLK plots
are scrambled; `gl.map.popcluster` is the same code as `gl.map.snmf`
before #102, with bars drawn at the wrong population's centre.

What works well: q-matrix rows carry the right individual and population
(checked against `indNames(x)`), the best run per K is PopCluster's own
choice, and in `gl.plot.popcluster` every bar sits above its own name.

## Findings

**F1 [HIGH, confidence: high] — likelihood table kept as text; LogL plot upside down (DOC5)**
`R/gl.run.popcluster.r:423-474` — the `.K` file is split into strings and
never converted, so `best_run` columns are character (`"-3.37E+02"`,
`"-"` for missing) and each plot has a discrete y axis sorted
alphabetically.
Failure scenario: K = 1:3 on testset.gl: LogL_Mean is −337.3, −283.5,
−259.4; the plot draws K = 1 at the top and K = 3 at the bottom, the
opposite of the values, which points a user to the wrong K. DLK1 shows
"-" as a category and 0.173 below 0.089.
Proposed change: convert the numeric columns (`"-"` → NA) and plot on
continuous axes. **Consequence: `best_run` columns K, LogL_*, DLK*,
FST.FIS become numeric; the four plots change (now correct).**

**F2 [MEDIUM, confidence: high] — files written to the working directory; cleanup never runs (FS7, DOC5)**
`R/gl.run.popcluster.r:119, 272-346, 566-578` — `output.path = getwd()`
receives the `.dat` and `.PcPjt` input files; the `cleanup` code comes
after `return()`.
Failure scenario: a default call leaves `output.popcluster.dat` and
`output.popcluster.PcPjt` in the working directory, and every call leaves
a temporary folder with a copy of the binary and all run files.
Proposed change: `output.path` defaults to `tempdir()`; remove the
temporary folder when `cleanup = TRUE`. **Consequence: nothing is
written to the working directory by default.**

**F3 [LOW, confidence: high] — plot arguments ignored (PLT1, PLT2, DEP1)**
`R/gl.run.popcluster.r:140, 445-490` — `plot_theme` is never used (each
panel adds `theme_dartR()`); `plot.dir` is not resolved with
`gl.check.wd()`; the panels are combined with `gridExtra` (Suggests,
unguarded).
Failure scenario: `plot_theme = theme_bw()` has no effect; without
`gridExtra` the call stops after all runs.
Proposed change: apply `plot_theme` to the four plots; combine them with
patchwork (Imports); resolve `plot.dir`. The returned `plots` list keeps
its four names.

**F4 [LOW, confidence: high] — PopCluster output ignores `verbose` (VRB1, VRB3)**
`R/gl.run.popcluster.r:393-417` — `system()` sends PopCluster's log and
the `chmod` calls straight to the terminal.
Failure scenario: `verbose = 0` prints about 150 lines for K = 2.
Proposed change: show program output only at `verbose >= 3`.

**F5 [MEDIUM, confidence: high] — argument and binary checks (FS5, VRB2)**
`R/gl.run.popcluster.r:150-235, 348-388` — `stringr` and `pillar` (both
Imports) are guarded with `cat()` + `return(-1)`; the two migration-model
checks overlap and the first uses `&&` where either condition should
fail; a missing binary prints with `cat()` and stops with an empty error;
`minK > maxK` is not checked; an unknown OS leaves the binary name
undefined; `build = "Jody"`.
Failure scenario: `minK = 3, maxK = 2` runs PopCluster, then stops with
"cannot open the connection"; a wrong `popcluster.path` gives an error
with an empty message.
Proposed change: one set of `stop(error(...))` checks before any file is
written: binary present (with the download source), OS supported,
`minK <= maxK` whole numbers, `rep >= 1`, migration-model conditions;
drop the Imports guards and `build =`.

**F6 [LOW, confidence: medium] — population assigned by row position (DAT2)**
`R/gl.run.popcluster.r:538-540` — `Label` is taken from the individual
index in the file but `Pop` from `x$pop` by row position.
Failure scenario: none with PopCluster 1.5.0 (it lists individuals in
index order); a version that lists them in another order would pair
individuals with the wrong population.
Proposed change: take `Pop` by the same index as `Label`.

**F7 [LOW, confidence: high] — gl.plot.popcluster inputs, verbosity, CRAN NOTE (FS2, FS5, VRB1)**
`R/gl.plot.popcluster.r:65-135` — `verbose = 2` default ignores
`gl.set.verbosity()`; `plot.K` is not checked (a K not run gives
"subscript out of bounds", several values a recycling warning and the same
error); a palette function in `color_clusters` fails with "Insufficient
values"; `aes(x = factor(Order), y = values, fill = K)` uses bare column
names (the "no visible binding for global variable" NOTE in every
`R CMD check`); `build = "Jody"`.
Failure scenario: `gl.set.verbosity(0)` still prints start and end lines.
Proposed change: `verbose = NULL`; `plot.K` must be one K of the run;
palette function and colour count; `.data` in `aes()`.

**F8 [HIGH, confidence: high] — gl.map.popcluster is the pre-#102 gl.map.snmf (DOC5)**
`R/gl.map.popcluster.r` — the body is identical to `gl.map.snmf` before
#102 (bars at another population's centre for non-alphabetical levels,
`movepops` by column position, extra populations in `x` fail).
Failure scenario (baseline): `pop(x)` levels `EmmacBurnBara,
EmmacBrisWive, EmmacBurdMist`: EmmacBrisWive's bars at EmmacBurdMist's
centre; `movepops` lon shift moves latitude on testset.gl.
Proposed change: call `gl.map.snmf` (fixed in #102; the PopCluster
q-matrix has the same `Pop_1..K` and `Label` columns), keeping the
signature and `list(Q_name, map)`; adds `plot.out` and `verbose`. This
PR is then stacked on #102.
**Consequence: as in #102 — maps change for non-alphabetical levels and
`lat`-first `movepops`; default colours change; extra populations in `x`
work.**

**F9 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7)**
All three files — `output.path` and `plot.dir` described as the working
directory; `@return` does not describe `best_run` columns or `matrix`
tables; the best run per K is PopCluster's choice, replicates are not
averaged; `gl.map.popcluster` `@return` duplicated; no `@family`; authors
not in Author(s)/Custodian form (DOC7, proposed rule). The PopCluster
binary prints an expiry date (2027-09-25 for the build tested), which the
help does not mention.
Proposed change: rewrite the three headers; regenerate Rd.

## Proposed changes

1. Numeric likelihood table, continuous plots (F1). **Consequence:
   `best_run` columns become numeric; the plots change (now correct).**
2. `output.path` default `tempdir()`; `cleanup` honoured (F2).
   **Consequence: no files in the working directory by default.**
3. `plot_theme` applied; patchwork instead of gridExtra; `plot.dir` via
   `gl.check.wd()` (F3).
4. PopCluster output only at `verbose >= 3` (F4).
5. Argument and binary checks with `stop(error())`; drop Imports guards
   and `build =` (F5).
6. `Pop` by individual index (F6).
7. `gl.plot.popcluster`: `verbose = NULL`, `plot.K` and colour checks,
   palette function, `.data` in `aes()` (F7).
8. `gl.map.popcluster` calls `gl.map.snmf`, stacked on #102 (F8).
   **Consequence: as in #102.**
9. Documentation for the three, Rd regenerated, NEWS entry (F9).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on all three.
- Spec: real PopCluster runs (K = 1:3 x 2, single K, wrong path,
  minK > maxK, default output path); likelihood plots read from the ggplot
  build; labels and populations checked against `x`; plot and map on a
  hand-built run object.
- Linux, Windows and `parallel = TRUE` (MPI): SKIPPED — macOS only.
- Models 1, 3, 4, `PopFlag = 1`, relatedness/kinship options: SKIPPED —
  default model 2 only.
- dartR Google Group / GitHub issues search: not run (no search access).
- Downstream callers: dartr2shiny runs all three (`input_Tabs.csv`,
  `variables_matrix.csv`: `Mypopcluster` passed to the plot, `Myqmat` to
  the map, a `checkKpopcluster` step on the run result). Every proposed
  change keeps the returned list names; change 1 makes `best_run$K`
  numeric, which compares equal to the same K given as number or text in
  R. No callers in other `dartR.*` packages.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |
| 9 | approved | Luis | |

## Outcome

- Changes 1–9 applied in commit e487727 on `review-popcluster`, PR #103 to
  `dev`; the branch also contains #102 (merged in, commit 140896e), so #102
  must merge first.
- Characterization test: 6 tests, 29 expectations pass with PopCluster
  1.5.0; every diff from baseline is tagged `[approved n]`.
- Old vs new where no fix applies: Q matrices, best runs, likelihood values,
  bar plot data and map rectangles identical.
- `verbose = 0`: no terminal output; default call writes nothing to the
  working directory.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; the global
  variables NOTE from `gl.plot.popcluster` is gone; `NAMESPACE` unchanged.

## Machine block

```json
{
  "function": "gl.run.popcluster",
  "also_covers": ["gl.plot.popcluster", "gl.map.popcluster"],
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "7ddbca9",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "DAT2", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 9}
  ],
  "coverage_skipped": [
    "Linux, Windows, MPI: macOS only",
    "models other than 2, PopFlag = 1: not exercised",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 103
}
```
