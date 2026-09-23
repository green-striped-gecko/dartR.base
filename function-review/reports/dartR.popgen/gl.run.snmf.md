# Review: gl.run.snmf + gl.plot.snmf + gl.map.snmf (dartR.popgen)
- Family mode: analysis (sNMF through LEA, bar plot, map)
- Scope: one review, one report, one PR for the three functions (agreed with
  Luis, 2026-09-23); three manifest rows
- Custodian of all three: Ching Ching Lau (STY5: this report is the
  discussion record; changes approved by Luis)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: a16da26 (origin/dev)
- Datasets: testset.gl, three populations (31 individuals, 18-20 loci); LEA
  sNMF runs (no external binary)
- Baseline: tests/testthat/test-gl.snmf.R (runs anywhere with LEA)

## Verdict

**Standards: Needs work** — `gl.run.snmf` prints LEA's output at every
verbosity, ignores `plot.out`, saves into the working directory and never
cleans up; `gl.plot.snmf` defaults `verbose` to 2 and uses `aes_()`;
`gl.map.snmf` has no `verbose`, checks or plot switch.
**Spec: Needs work** — `gl.map.snmf` is a copy of `gl.map.structure` before
#95 and #99: with non-alphabetical population levels every population's
bars are drawn at another population's centre, `movepops` shifts by column
position, and an `x` with extra populations fails with an opaque error.
`gl.run.snmf`'s `cleanup` does nothing.

What works well: `gl2geno` writes individuals in `indNames(x)` order, so
every q-matrix row carries the right label (checked against the genotype
matrix); the best run per K by cross-entropy is chosen correctly; in
`gl.plot.snmf` every bar sits above its own individual's name.

## Findings

**F1 [HIGH, confidence: high] — gl.map.snmf repeats the gl.map.structure bugs fixed in #95 and #99 (DOC5, FS5)**
`R/gl.map.snmf.r:76-198` — population centres are reordered with
`centers[match(rownames(centers), levels(...)), ]` (the inverse
permutation); `movepops` is added by column position; a population
mismatch is a `message()` followed by an opaque failure; colours are
`rainbow()`; no `verbose`, no checks, the map is always printed.
Failure scenario (baseline): `pop(x)` levels `EmmacBurnBara,
EmmacBrisWive, EmmacBurdMist`: EmmacBrisWive's bars are drawn at
EmmacBurdMist's centre under the right label. `movepops = data.frame(lon
= c(1, 0, 0), lat = 0)` on testset.gl (latlon stored `lat`, `lon`) moves
the first population north instead of east. An `x` with a fourth
population stops with "row names contain missing values" after many
leaflet warnings.
Proposed change: convert the sNMF q-matrix to the `gl.plot.structure`
layout (Label, cluster1..K, K, orig.pop, ord) and draw with
`gl.map.structure`, keeping the signature (`color_clusters` passed as
`plot.colors`) and the return `list(Q_name, map)` (per-population tables of
the original q-matrix rows, in map order); add `plot.out` and `verbose`.
**Consequence: maps change for `x` with non-alphabetical population
levels and for `movepops` on genlights storing `lat` first; default
colours change from `rainbow()` to the dartR palette; `x` with extra
populations now works.**

**F2 [MEDIUM, confidence: high] — gl.run.snmf: `cleanup` never runs (DOC5)**
`R/gl.run.snmf.r:173-183` — the cleanup code comes after `return()`.
Failure scenario: `cleanup = TRUE` (default) leaves the LEA project
folder with every run of every K in `tempdir()`; each call adds another.
Proposed change: remove the folder before returning when `cleanup = TRUE`;
return `best_run` as the run names (e.g. `"K2/run1"`) when the files are
removed, and as full paths when `cleanup = FALSE`.
**Consequence: with the default, `best_run` holds run names instead of
paths to folders that are removed.**

**F3 [MEDIUM, confidence: high] — gl.run.snmf ignores `plot.out` and saves into the working directory (PLT2, FS7)**
`R/gl.run.snmf.r:150-164` — the cross-entropy plot is always drawn;
`plot.dir` is not resolved with `gl.check.wd()`, so `plot.file` writes to
the working directory.
Failure scenario: `plot.out = FALSE, plot.file = "ce"` still draws the
plot and writes `ce.RDS` into the working directory (baseline).
Proposed change: resolve `plot.dir` with `gl.check.wd()`; with
`plot.out = FALSE` draw off-screen so `cross_entropy` is still a
`recordedplot` (dartr2shiny uses the returned object).

**F4 [LOW, confidence: high] — gl.run.snmf prints LEA's output at every verbosity (VRB1, VRB3)**
`R/gl.run.snmf.r:104-121` — `LEA::snmf` and `LEA::cross.entropy` print
their logs.
Failure scenario: `verbose = 0` prints 544 lines for K = 1:3 with three
replicates.
Proposed change: capture LEA's output below `verbose = 3`.

**F5 [LOW, confidence: high] — gl.run.snmf house structure (FS3, FS5, DEP1)**
`R/gl.run.snmf.r:76-97` — `build = "Jody"`; the LEA guard returns -1
after a `cat()` although LEA is in Imports; `minK`, `maxK`, `rep` are not
checked.
Failure scenario: `minK = 3, maxK = 2` reaches LEA with `K = 3:2`.
Proposed change: drop `build =` and the guard; check `minK <= maxK`, both
whole numbers >= 1, `rep >= 1`.

**F6 [LOW, confidence: high] — gl.plot.snmf input handling and verbosity (FS2, FS5, VRB1)**
`R/gl.plot.snmf.r:73-123` — `verbose = 2` default ignores
`gl.set.verbosity()`; `gl.colors()` is called without `verbose`, so it
prints at `verbose = 0`; `plot.K` with several values silently uses the
first; a K not in the run stops with "no applicable method for
'pivot_longer'"; a palette function in `color.clusters` fails with
"Insufficient values".
Failure scenario: `gl.set.verbosity(0); gl.plot.snmf(r, plot.K = 3)`
prints "Starting gl.plot.snmf", the gl.colors lines and "Completed".
Proposed change: `verbose = NULL`; quiet `gl.colors()`; `plot.K` must be a
single K present in the run (clear `error()`); accept a palette function;
check colour count.

**F7 [LOW, confidence: high] — gl.plot.snmf dendrogram and deprecated aes_() (DOC5, PLT1)**
`R/gl.plot.snmf.r:137, 306` — `hclust(dist(res))` clusters Euclidean
distances between rows of the distance matrix (the defect fixed in
`gl.plot.structure` in #94); `aes_()` is deprecated.
Failure scenario: `den = TRUE` orders individuals by a tree that is not
the clustering of their distances; the first call prints a ggplot2
deprecation warning asking users to report to the dartR group.
Proposed change: `hclust(res)` directly; `aes()` with `.data`.
**Consequence: individual order in `den = TRUE` plots changes.**

**F8 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7)**
All three files — `gl.run.snmf`: `@name`/`@title` split over lines,
`plot.dir` "working directory", `@return` does not describe `matrix`
columns or that only the best run per K is kept (no averaging over
replicates), `@export` twice; `gl.plot.snmf`: `plot.K` must be one value,
`@return` "Q-matrix"; `gl.map.snmf`: `@return` mentions PopCluster, two
`@return` blocks, `movepops` order not stated. No `@family`; authors not in
Author(s)/Custodian form (DOC7, proposed rule).
Failure scenario: a user expects replicate averaging (as in
`gl.plot.structure`) and gets the single best run.
Proposed change: rewrite the three headers; regenerate Rd.

## Proposed changes

1. `gl.map.snmf` draws with `gl.map.structure` (F1). **Consequence: maps
   change for non-alphabetical population levels and `lat`-first
   `movepops`; default colours change; extra populations in `x` work.**
2. `cleanup` honoured; `best_run` = run names when cleaned (F2).
   **Consequence: default `best_run` holds names, not paths.**
3. `plot.out` honoured (off-screen recording), `plot.dir` via
   `gl.check.wd()` (F3).
4. LEA output only at `verbose >= 3` (F4).
5. House structure and argument checks in `gl.run.snmf` (F5).
6. `gl.plot.snmf`: `verbose = NULL`, quiet `gl.colors()`, `plot.K` and
   colour checks, palette function (F6).
7. `gl.plot.snmf`: `hclust` on the distances, `aes()` (F7).
   **Consequence: `den = TRUE` order changes.**
8. Documentation for the three, Rd regenerated, NEWS entry (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on all three.
- Spec: real LEA runs (K = 1:3, replicates, seed); row-to-label mapping
  checked against the `.geno` file; bar/label alignment read from the
  ggplot build; map positions read from the leaflet widget.
- Replicate averaging: not a defect — the functions keep the best run per
  K by cross-entropy by design; documented under F8.
- `ploidy_lv` other than 2: SKIPPED — no polyploid fixture.
- dartR Google Group / GitHub issues search: not run (no search access).
- Downstream callers: dartr2shiny runs all three (`input_Tabs.csv`,
  `variables_matrix.csv`, `template_report.csv`: `Mysnmf`, `Myqmat`,
  `Mymap.snmf`); every proposed change keeps the returned list names and
  the `recordedplot` class of `cross_entropy`. No callers in other
  `dartR.*` packages.

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

## Outcome

- Changes 1–8 applied in commit 1654a62 on `review-snmf`, PR #102 to `dev`.
- Characterization test: 10 tests, 36 expectations pass (real LEA runs);
  every diff from baseline is tagged `[approved n]`.
- Old vs new where no fix applies: seeded run matrices and best runs,
  `gl.plot.snmf` returned Q and plot data, and maps (alphabetical levels,
  with and without `movepops`) are identical.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; `NAMESPACE`
  unchanged; no new NOTE.

## Machine block

```json
{
  "function": "gl.run.snmf",
  "also_covers": ["gl.plot.snmf", "gl.map.snmf"],
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "a16da26",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "PLT2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8}
  ],
  "coverage_skipped": [
    "ploidy other than 2: no fixture",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 102
}
```
