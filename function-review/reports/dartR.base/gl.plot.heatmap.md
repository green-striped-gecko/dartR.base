# Review: gl.plot.heatmap (dartR.base)

## Provenance

- Model: Claude Fable 5.1 (claude-fable-5-1)
- Skill: dartr-function-review v2.0.0 (Phase A, read-only review)
- Package commit: 250b840 (dev_luis after merging origin/dev f9ee087;
  `R/gl.plot.heatmap.r` identical on both)
- Date: 2026-09-15
- Family mode: graphics (pure plotting wrapper; report-mode checks
  applied: input untouched, no history append, nothing computed that
  plotting could hide)
- Datasets: testset.gl[1:12, ] (7 populations) for individual distances,
  testset.gl for `gl.dist.pop()` and `gl.fixed.diff()`; hand-built
  matrices (non-zero diagonal, asymmetric cell, non-square, renamed
  individuals) for the matrix path
- Baseline: tests/testthat/test-gl.plot.heatmap.R (new file; 13 tests,
  34 expectations, all passing at 250b840; defects pinned as-is and tagged
  with finding IDs)
- Related review: utils.heatmap (Arthur, PR #327, merged 2026-09-14) is
  the vendored `gplots::heatmap.2` fork this function calls. Its report
  flagged this function's verbose-0 leak for this review.
- Checks skipped: dartR Google Group not searched (no browser session);
  GitHub issue search across the org found no report naming this
  function.

## Verdicts

**Standards: Needs work** -- the datatype gate, flag start and flag end
conform, but the default palette prints three lines at `verbose = 0`, the
legend leaves `par(mar)` at 1 line on every side, the dendextend guard
returns `-1` instead of stopping, the population colour vector is not
validated, and the documentation describes a different function (a
`gplots::heatmap.2` wrapper returning `NULL`).

**Spec: Needs work** -- a `dist` object is drawn correctly, but a
`matrix` is silently coerced through `as.dist()`, which discards the
diagonal and the upper triangle; two populations sharing a colour get a
wrong legend swatch; an `fd` object with `x` crashes; and `x` whose
individuals do not match `D` is dropped without a message or fails with
an error from the helper.

What works well: the `dist` path reproduces `hclust` ordering exactly,
the input object is never modified, `plot.out = FALSE` draws nothing,
and `x` is used only through `indNames()`/`pop()`, so an FBM-backed
object is never densified.

## Findings

**F1 [HIGH, confidence: high] -- matrix input loses its diagonal and
upper triangle (FS5; DOC5 proposed rule)**
`R/gl.plot.heatmap.r:190` -- `D <- as.dist(D)` runs for every `matrix`
input. `as.dist()` keeps the lower triangle only and sets the diagonal to
zero; a non-square matrix only warns ("non-square matrix") and the plot
is drawn anyway.
Failure scenario: `dartR.captive::gl.run.EMIBD9()` and
`gl.relatedness()` pass full relatedness matrices whose diagonal holds
self-relatedness (or `NA`); the heatmap shows a zero diagonal in the
colour of "unrelated", and `diag.na = TRUE` has nothing left to hide. A
matrix with a cell of 99 in the upper triangle plots without it
(baseline test "matrix input"). `gl.report.fstat()` matrices are full
and symmetric, so they are unaffected.
Proposed change: treat a matrix as given. Stop with a clear message
when it is not square or its row and column names differ; when one
triangle is entirely `NA`, mirror the other so half-filled matrices
still plot; otherwise draw the values as supplied. `dist` input is
unchanged.
**Consequence: the plot changes for matrix inputs with a non-zero
diagonal or asymmetric entries.**

**F2 [MEDIUM, confidence: high] -- `verbose = 0` prints three lines
(VRB1, VRB3)**
`R/gl.plot.heatmap.r:82` -- the default `palette.divergent =
gl.colors("div")` is evaluated inside the function at the session
verbosity, so "Starting gl.colors", "Selected color type div" and
"Completed: gl.colors" print at every verbosity, including 0. Callers
that forward `verbose = 0` (`gl.report.fstat()`, both dartR.captive
functions) inherit the leak. A palette passed explicitly is silent
(baseline test "verbose = 0 leaks").
Proposed change: evaluate the default silently
(`gl.colors("div", verbose = 0)`); the palette is the same.

**F3 [MEDIUM, confidence: high] -- legend swatches are wrong when two
populations share a colour (correctness; no catalogue rule)**
`R/gl.plot.heatmap.r:184-185` -- `legend_text <- unique(pop)` and
`legend_color <- unique(color)` are built independently. With colours
`red, red, blue, ...` for seven populations the legend has 7 labels and
6 swatches, and `legend()` recycles: the seventh population is labelled
orange while its side bar is red (baseline test "legend swatches").
Proposed change: build the legend from the population-to-colour table,
one row per population, so labels and swatches always pair.

**F4 [MEDIUM, confidence: high] -- `fd` input with `x` crashes (FS5,
FS6)**
`R/gl.plot.heatmap.r:177` -- the population-colour block calls
`as.matrix(D)` before the datatype dispatch. For an `fd` object (a
list) this gives `data.frame(): arguments imply differing number of
rows: 0, 1`.
Failure scenario: `gl.plot.heatmap(gl.fixed.diff(gl), x = gl)`.
Proposed change: extract the plotted matrix first (`D$fd` for `fd`),
then run the colour block only for individual-level input; for `fd`
ignore `x` with a note at `verbose >= 2` (fixed differences are
population-level, so population side bars do not apply).

**F5 [MEDIUM, confidence: high] -- `x` that does not match `D` is
dropped silently or fails in the helper (FS5; VRB4 proposed rule)**
`R/gl.plot.heatmap.r:178-196` -- individuals are matched by name through
`merge()`, then the only check is `ncol(m) != nInd(x)`. Fewer
individuals in `D` than in `x`: colours and legend vanish with no
message at any verbosity. Same count but different names: `merge()`
returns no rows and `utils.heatmap()` stops with "ColSideColors must be
a character vector of length ncol(x)" (baseline test "x whose
individuals do not match D").
Proposed change: match `D`'s dimnames against `indNames(x)` explicitly;
when every column has a match, colour by population (regardless of
`nInd(x)`); otherwise drop the colours with a warning at `verbose >= 1`
naming how many columns were not found.

**F6 [LOW, confidence: high] -- `par(mar)` is changed and never restored
(STY1; graphics hygiene)**
`R/gl.plot.heatmap.r:238` -- `par(mar = c(1, 1, 1, 1))` before
`legend()` persists after the call, so the user's next base-graphics
plot has one-line margins (baseline test "legend leaves par(mar)").
Proposed change: `op <- par(mar = c(1, 1, 1, 1)); on.exit(par(op), add
= TRUE)`.

**F7 [LOW, confidence: high] -- `palette_discrete` is not validated
(FS5)**
`R/gl.plot.heatmap.r:164-171` -- a colour vector of the wrong length
fails on `names(colors_pops) <-` with "'names' attribute [7] must be
the same length as the vector [2]" (baseline test "palette_discrete of
the wrong length").
Proposed change: stop with `error()` when a vector is not `nPop(x)`
long; a function is called with `nPop(x)` as now.

**F8 [LOW, confidence: high] -- dependency guard returns `-1` instead of
stopping (DEP1)**
`R/gl.plot.heatmap.r:137-144` -- `cat(error(...)); return(-1)`. A
caller that assigns the result (`gl.report.fstat()`,
`gl.relatedness()`) receives `-1` as its plot and carries on.
Proposed change: `stop(error(...))`, the house idiom.

**F9 [LOW, confidence: high] -- documentation describes a different
function (DOC1, DOC2, DOC5 proposed rule, DOC7 proposed rule)**
`R/gl.plot.heatmap.r:5-78`:
- `@description` and `@param ...` say the function wraps
  `gplots::heatmap.2`; it calls `utils.heatmap()`, the vendored fork,
  and gplots is not used anywhere in the package.
- `@return` says `NULL`; the function returns the `utils.heatmap()` list
  invisibly (`rowInd`, `colInd`, `carpet`, dendrograms, `breaks`, ...)
  when `plot.out = TRUE`.
- `legendy` documented as `[default 1]`, code default `0.5`;
  `legendx`/`legendy` are 0-1 coordinates inside the colour-key panel,
  which is not stated.
- `margins` uses `[Default = ...]`; `revC`, `cexRow`, `cexCol`,
  `srtRow`, `key.title`, `key.xlab`, `key.ylab`, `main`, `xlab`, `ylab`
  have no `[default]`; `cexRow`/`cexCol`/`srtRow` are "Integer value"
  but take numerics.
- `verbose` text is not the DOC2 standard; no `@details`; `@author`
  has no `Author(s):` part (DOC7).
- `@examples`: the second block is guarded on gplots (not the
  dependency; dendextend is) and repeats the `possums.gl` example from
  the `\donttest` block.
- `@importFrom gtools invalid` imports a function this file never
  calls; it is the only gtools reference in the package, so gtools sits
  in Imports for nothing (dropping it from DESCRIPTION is the
  custodian's call).
Proposed change: docs-only rewrite of the header (description, return,
defaults, coordinate system, DOC2 verbose text, `@details`, author
block, examples) and removal of the gtools import line; run
`devtools::document()`.

**F11 [LOW, confidence: high] -- addendum: gtools declared in Imports
but unused after F9 (DEP2; R CMD check NOTE)**
`DESCRIPTION:41` -- once `@importFrom gtools invalid` is removed, no
file in the package references gtools, and check reports "Namespace in
Imports field not imported from: 'gtools'".
Proposed change: drop gtools from Imports.

**F10 [INFO, confidence: high] -- argument naming (PLT1)**
`palette_discrete` uses an underscore beside `palette.divergent`; the
house names are `plot.colors`-style dotted. dartr2shiny passes
`palette_discrete` by name, so a rename is an API2 change; no change
proposed, recorded for the API pass.

## Proposed changes

1. Handle matrix input as a matrix: no `as.dist()`; stop on non-square
   or mismatched dimnames; mirror when one triangle is all `NA` (F1).
   **Consequence: the plot changes for matrix inputs with a non-zero
   diagonal or asymmetric entries.**
2. Evaluate the default palette silently so `verbose = 0` prints nothing
   (F2).
3. Build the legend one row per population so labels and swatches pair
   (F3).
4. Extract the plotted matrix before the colour block; ignore `x` for
   `fd` input with a note at `verbose >= 2` (F4).
5. Match `D` columns to `indNames(x)` by name; warn at `verbose >= 1`
   and drop the colours when any column is unmatched (F5).
6. Restore `par(mar)` after the legend (F6).
7. Validate the length of a `palette_discrete` vector (F7).
8. Make the dendextend guard `stop()` (F8). **Consequence: a session
   without dendextend now errors instead of returning `-1`.**
9. Documentation rewrite and removal of the unused gtools import; docs
   only (F9).
10. Addendum: drop gtools from DESCRIPTION Imports (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: behaviour vs roxygen on `dist`, `matrix`, `fd` input, with and
  without `x` -- run
- Family checks (graphics/report mode): input untouched (run, identical);
  no history append (run, none); `plot.out = FALSE` (run, returns NULL,
  draws nothing)
- FBM path (DAT6): `x` is read only through `indNames()` and `pop()`;
  no densification of `x` -- checked by reading, not run
- SilicoDArT `x`: not run; the function reads only individual and
  population names from `x`, which do not depend on datatype
- Sibling callers (API3): `gl.report.fstat()` (dartR.base),
  `gl.relatedness()` and `gl.run.EMIBD9()` (dartR.captive), dartr2shiny
  `slot_exceptions.csv` -- read; none passes `x` with `fd`; both
  dartR.captive callers pass full matrices (affected by change 1)
- Google Group: SKIPPED -- no browser session

## Report notes (other functions, not fixed here)

- `utils.heatmap()` emits 13 warnings per call ("'x' is NULL so the
  result will be NULL") whenever `colRow`/`colCol` are `NULL`, its own
  default: `dendextend::set("labels_col", NULL)` at lines 743 and 766.
  Every `gl.plot.heatmap()` call without `x` inherits them. Belongs to
  the utils.heatmap custodian (PR #327 is merged; the manifest row still
  says `pr-open`).
- `dartR.captive::gl.relatedness()` passes this function's return value
  (a list) to `utils.plot.save()`, which expects a ggplot; recorded for
  that function's review.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | behaviour change approved as stated |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |
| 8 | approved | Luis | behaviour change approved as stated |
| 9 | approved | Luis |  |
| 10 | approved | Luis | addendum, accepted with the pre-push OK |

## Outcome

Branch `review-gl.plot.heatmap` from origin/dev f9ee087 (an open PR, #402,
already sits on dev_luis).

- Change 1 (F1): `as.dist()` coercion removed; a matrix is drawn as
  supplied, mirrored when one triangle is all NA, and stops when
  non-square or when row and column names differ. Evidence: test "matrix
  input is plotted as supplied" (diagonal of 5 and the 99 cell present,
  upper-NA matrix mirrored equals the full distance matrix); test
  "non-square or name-mismatched matrix stops". NEWS entry added.
- Change 2 (F2): default `palette.divergent = gl.colors("div", verbose =
  0)`. Evidence: `capture.output()` at verbose 0 is empty with and without
  `x`; verbose 2 prints exactly the start, datatype and end lines;
  `gl.report.fstat(platypus.gl, verbose = 0, plot.display = TRUE)` prints
  0 lines.
- Change 3 (F3): legend built from the population-to-colour table.
  Evidence: mocked `legend()` receives 7 labels and 7 fills equal to the
  supplied palette for `red, red, blue, ...`; side colours (mocked
  `utils.heatmap()`) equal the palette indexed by each individual's
  population.
- Change 4 (F4): plotted matrix extracted before the colour block; `x`
  ignored for `fd` with a note at verbose 2. Evidence: `gl.plot.heatmap(fd,
  x = gl)` returns the list; verbose 2 prints "D is a population-level
  matrix; x is ignored".
- Change 5 (F5): columns matched to `indNames(x)` by name. Evidence: an
  8-individual subset of a 12-individual object is coloured (8 side
  colours, legend limited to the populations present); renamed
  individuals give the verbose-1 warning "12 of 12 columns of D are not
  individuals of x" and a plot with no colours; nothing at verbose 0.
- Change 6 (F6): `par(mar)` restored via `on.exit()`. Evidence:
  `par("mar")` identical before and after a call with legend.
- Change 7 (F7): `palette_discrete` vector length checked. Evidence: a
  2-colour vector for 7 populations stops with "...vector of 7 colours,
  one per population of x; got 2".
- Change 8 (F8): guard is `stop(error(...))`. Evidence: by reading; a
  missing dendextend cannot be simulated in the test session.
- Change 9 (F9): roxygen header rewritten, `@importFrom gtools invalid`
  removed (NAMESPACE loses `importFrom(gtools,invalid)`; gtools stays in
  DESCRIPTION Imports, custodian's call), `devtools::document()` run;
  the example runs clean via `tools::Rd2ex()` including the `\donttest`
  block.

Characterization test: 14 tests, 53 expectations, 0 failures. Every
flipped pin carries an `[approved diff, change n]` tag: changes 1-7. No
unexplained diff. Caller tests unchanged: test-gl.report.fstat.R 27/27,
test-utils.heatmap.R 4/4. Function run end to end at verbose 3 on `dist`
+ `x`, matrix + `x`, `fd` + `x` and `gl.dist.pop()` + `x`.

PR #406 (review-gl.plot.heatmap -> dev).

Package check (`R CMD check --no-tests`): no finding from the package
code except one new NOTE, "Namespace in Imports field not imported from:
'gtools'", caused by change 9 removing the only gtools import. Addendum
F11 (approved by Luis with the pre-push OK, recommended option): gtools
dropped from DESCRIPTION Imports; no file in R/, tests/ or inst/ uses it.
The remaining check output (non-portable file names, top-level files,
"built under R 4.4.3" install warnings) comes from local scratch files
and the R installation, not from the package.

```json
{
  "function": "gl.plot.heatmap",
  "package": "dartR.base",
  "family": "graphics",
  "skill_version": "2.0.0",
  "commit": "250b840",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1",
     "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "correctness",
     "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "STY1",
     "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS5",
     "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DEP1",
     "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1",
     "status": "approved", "change": 9},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "PLT1",
     "status": "noted", "change": null},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DEP2",
     "status": "approved", "change": 10}
  ],
  "coverage_skipped": ["Google Group: no browser session",
                       "DAT6: checked by reading only",
                       "SilicoDArT x: not run"],
  "status": "pr-open",
  "pr": 406
}
```
