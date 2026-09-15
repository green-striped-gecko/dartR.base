# Review: gl.map.interactive (dartR.base)

## Provenance

- Model: Claude Fable 5.1 (claude-fable-5-1)
- Skill: dartr-function-review v2.0.0 (Phase A, read-only review)
- Package commit: be7cb2e (origin/dev; `git diff origin/dev --
  R/gl.map.interactive.r` empty on dev_luis at 39c7ed2, so `load_all()`
  exercised the reviewed code)
- Date: 2026-09-15
- Family mode: report (plotting; the input object must come back
  untouched, no history append)
- Datasets: platypus.gl (3 pops, lat,lon column order), bandicoot.gl
  (5 pops, lon,lat column order), testset.gs (SilicoDArT),
  platypus.gl[1:6,] with hand-built individual matrices,
  gl.dist.pop(platypus.gl), gl.dist.ind(platypus.gl[1:6,]),
  gl.gen2fbm(platypus.gl), a synthetic 20x20 EPSG:4326 GeoTIFF for the
  raster branch
- Baseline: tests/testthat/test-gl.map.interactive.R (new file; 46
  assertions, all passing at be7cb2e; defects pinned as-is and tagged
  with finding IDs). Introspection is on the leaflet widget's call list
  (`m$x$calls`), not on rendered tiles.

## Verdicts

**Standards: Needs work** -- the entry scaffold (verbosity, flag start,
datatype gate, flag end) conforms and the function never touches the
input object, but the dependency guards return `-1` instead of stopping,
`terra` and `scales` run unguarded on the raster path, and `gl.colors`
prints three lines at `verbose = 0`.

**Spec: Rework** -- the map of individuals and population labels is
correct on multi-population data, but every path through the `matrix`
argument is broken in some way: `dist` objects (what `gl.dist.pop` and
`gl.dist.ind` return) crash, individual-level links attach to the wrong
individuals, zero-valued pairs and self-links are drawn because the
`> 0` filter is dead, and a single-population object places its label
with latitude and longitude swapped.

What works well: circle coordinates, popups, colours and population
centres track `@other$latlon`, `indNames`, `pop` and `popNames` exactly,
on SNP, SilicoDArT and FBM-backed objects alike.

## Findings

**F1 [HIGH, confidence: high] -- single-population label placed with
lat/lon swapped (DOC5; positional indexing of named columns)**
`R/gl.map.interactive.r:165-168` -- with one population `apply` returns a
named vector and the code rebuilds `centers` as
`data.frame(lon = centers[1], lat = centers[2])`, assuming the first
column of `latlon` is `lon`. dartR's own datasets disagree on column
order: `platypus.gl` and `testset.gl` store `lat,lon`, `bandicoot.gl` and
`possums.gl` store `lon,lat`.
Failure scenario: `gl.map.interactive(gl.keep.pop(platypus.gl,
pop.list = "SEVERN_ABOVE"))` emits the label marker at
`lat = 151.52, lng = -29.48` (reproduced; baseline test "F1"). Latitude
151 is outside the valid range, so the label is never seen.
Proposed change: index by name (`centers["lon"]`, `centers["lat"]`), or
build `centers` with `drop = FALSE` semantics for every `nPop`.

**F2 [HIGH, confidence: high] -- `dist` objects crash with an opaque
error (FS5, DOC5)**
`R/gl.map.interactive.r:139` -- `nrow()` of a `dist` object is `NULL`, so
the dimension check fails with "argument is of length zero". `@param
matrix` says "A distance matrix", and the natural producers inside dartR
(`gl.dist.pop`, `gl.dist.ind`) return class `dist`.
Failure scenario: `gl.map.interactive(platypus.gl,
matrix = gl.dist.pop(platypus.gl))` errors before drawing anything
(reproduced; baseline test "F2"). Nothing tells the user to call
`as.matrix()`.
Proposed change: coerce `dist` input with `as.matrix()` at the top of the
matrix branch; stop with a clear message for anything that is neither
`dist` nor a square matrix.

**F3 [HIGH, confidence: high] -- individual-level matrix rows are
reordered without the columns, so links join the wrong individuals
(DAT2)**
`R/gl.map.interactive.r:219` -- `matrix <- matrix[order(indNames(x)), ]`
permutes rows into alphabetical individual order while columns stay in
object order, and the coordinates (`df`) also stay in object order.
Failure scenario: on `platypus.gl[1:6,]` (individual order T27, T35,
SDS4, SDS12, SUS20, SUS28) a matrix with a single non-zero cell between
individuals 1 and 2 draws two lines, one from individual 1 to 6 and one
from individual 2 to 5 (reproduced; baseline test "F3"). Any dataset
whose individuals are not already sorted alphabetically shows this;
population-level matrices skip the reorder and are drawn positionally,
so a matrix whose rows are in a different order than `popNames(x)` also
misattaches silently.
Proposed change: remove the row reorder; when the matrix carries row and
column names that match `indNames(x)` (or `popNames(x)`), align it with
`matrix[names, names]`; otherwise assume object order and say so in
`@param matrix`.

**F4 [MEDIUM, confidence: high] -- the `> 0` filter never fires, so every
pair including self-links is drawn (STY3; operator precedence)**
`R/gl.map.interactive.r:248-249` -- the condition is
`!is.null(v) | !is.na(v) & v > 0`. `&` binds tighter than `|`, and
`is.null()` of a matrix element is always `FALSE`, so the whole test is
always `TRUE`.
Failure scenario: three populations produce six polylines, three of them
zero-length self-links; six individuals with one non-zero cell produce
21 polylines instead of one (reproduced; baseline tests "population
matrix" and "F3"). With `standard = TRUE` a zero distance is rescaled to
1 and drawn in the palette's lowest colour, so the dead filter is masked
on the default path but not with `standard = FALSE`.
Proposed change: `if (i != ii && !is.na(v) && v > 0)`.

**F5 [LOW, confidence: high] -- `gl.colors` prints at `verbose = 0`
(VRB3)**
`R/gl.map.interactive.r:238` -- `gl.colors("div")` is called without
`verbose`, so it prints "Starting gl.colors", "Selected color type div",
"Completed: gl.colors" whenever a matrix is given.
Failure scenario: `gl.map.interactive(x, matrix = m, verbose = 0)`
prints three lines (reproduced; baseline test "F5").
Proposed change: `gl.colors("div", verbose = 0)`.

**F6 [MEDIUM, confidence: high] -- asymmetric branch crashes on `NA`
(FS5)**
`R/gl.map.interactive.r:277-278` -- the guard is `!is.null(matrix[i, ii])`
where `!is.na()` was clearly intended (the `else` branch assigns a grey
"missing" colour that is unreachable). An `NA` cell then reaches
`if (matrix[i, ii] > matrix[ii, i])`.
Failure scenario: any asymmetric matrix with a missing cell, `symmetric
= FALSE`, errors with "missing value where TRUE/FALSE needed"
(reproduced; baseline test "F6").
Proposed change: test `is.na()` on both cells; skip the pair (or draw it
in the existing grey) when either is missing.

**F7 [LOW, confidence: high] -- colour vector shorter than `nPop`
silently yields `NA` colours (DOC5)**
`R/gl.map.interactive.r:154-157` -- `cols[as.numeric(pop(x))]` indexes
past the end of a short vector; `@param ind.circle.cols` promises "as
many colors as there are populations".
Failure scenario: `ind.circle.cols = "red"` on `platypus.gl` colours one
population red and passes `NA` for the other 76 individuals (reproduced;
baseline test "F7"); leaflet falls back to its default blue without a
message.
Proposed change: stop with an informative error when
`length(cols) < nPop(x)` (recycling with a `verbose >= 2` warning is the
alternative).

**F8 [LOW, confidence: medium] -- `latlon` stored as a matrix crashes
(DAT5)**
`R/gl.map.interactive.r:159,176` -- `df$lon` fails with "$ operator is
invalid for atomic vectors" when `@other$latlon` is a matrix rather than
a data frame. `gl.compliance.check` renames `latlong`/`long` but does not
coerce the class.
Failure scenario: an object assembled outside `gl.read.dart` with
`latlon <- cbind(lat, lon)` errors at the first `$` (reproduced; baseline
test "F8").
Proposed change: `df <- as.data.frame(x@other$latlon)`.

**F9 [MEDIUM, confidence: high] -- dependency guards return `-1`; `terra`
and `scales` are unguarded (DEP1)**
`R/gl.map.interactive.r:106-124, 320, 88` -- the `leaflet` and
`leaflet.minicharts` guards `cat()` the error and `return(-1)`, so a
caller that assigns the result gets a number instead of an error;
`terra` (Suggests) is called on the raster path with no guard, and the
`raster.colors` default evaluates `scales::viridis_pal` (Suggests) on
that path.
Failure scenario: without `terra`, `raster.image = "x.tif"` fails with a
namespace-load error after the map has been built.
Proposed change: use the DEP1 `stop(error(...))` idiom for `leaflet` and
`leaflet.minicharts`; guard `terra` and `scales` the same way when
`raster.image` is not `NULL`.

**F10 [MEDIUM, confidence: high] -- `standard` promises line width, the
symmetric branch maps the value to colour only (DOC5)**
`R/gl.map.interactive.r:10-11, 250-256` -- `@param standard` says "line
width will be standardised to be between 1 to 10", but `addPolylines` is
called without `weight`, so every symmetric link has leaflet's default
width 5; the standardised value drives `qpal()` (colour) and the legend.
Only the asymmetric branch maps the value to thickness (via `addFlows`).
Failure scenario: a user sets `standard = FALSE` expecting raw distances
as widths and sees identical widths (reproduced: all polylines carry
`weight = 5`).
Proposed change: document what actually happens (colour scale 1-10 in
the symmetric branch, arrow thickness in the asymmetric branch); a
separate optional change maps the value to `weight` as well.

**F11 [LOW, confidence: high] -- roxygen gaps and errors (DOC1, DOC2,
DOC5)**
`R/gl.map.interactive.r:1-70` -- no `@description` tag; `@return` says
"plots a map" while the function returns a leaflet htmlwidget visibly;
the `verbose` text deviates from the DOC2 standard ("progress log");
`ind.circle.transparency` ends "Defaults to 0.8" instead of `[default
0.8]`; typos "vectot", "Size or circles", "Should individuals plotted";
`palette.links` and `legend.title` are documented "in case a matrix is
provided" but the asymmetric branch ignores both (colours are
hard-coded, no legend); `@param matrix` does not state that a `dist`
object is (after F2) accepted or that rows must follow object order
(after F3).
Proposed change: docs only, in one pass, followed by
`devtools::document()` (DOC4).

**F12 [LOW, confidence: high] -- `@author` has no Custodian label (DOC7,
proposed rule)**
`R/gl.map.interactive.r:57` -- "Bernd Gruber -- Post to ..." names an
author but no custodian.
Proposed change: `Author(s): Bernd Gruber. Custodian: Bernd Gruber --
Post to \url{...}` (custodian to be confirmed by the team).

**F13 [INFO, confidence: high] -- housekeeping (FS3, FS9, PLT1, STY3)**
- `:97` `build = "v.2023.2"` is the outdated FS3 idiom.
- `:170-171, 318` `addTiles()` (OpenStreetMap) and `addProviderTiles()`
  both load, so two base tile sets download for every map.
- `:312-314` "Completed:" prints before the provider tiles and raster
  are added; a `terra::rast()` failure appears after "Completed".
- `:73` the argument named `matrix` shadows `base::matrix` inside the
  body (no call to `matrix()` exists today).
- `:79` default `rainbow` bypasses the house palettes in `R/zzz.r`
  (PLT1); a change of default is an API1 matter and is not proposed.
- when every population has one individual, `nInd(x) == nPop(x)` and a
  matrix is treated as population-level; harmless because both index the
  same points, noted for completeness.
Proposed change: drop `build`, drop `addTiles()`, move the flag-end
block after the raster step.

## Proposed changes

1. Place the single-population label from named columns (`centers["lon"]`,
   `centers["lat"]`) (F1).
2. Accept `dist` objects via `as.matrix()`; stop with a clear message for
   input that is neither `dist` nor a square matrix (F2).
3. Drop the row-only reorder; align the matrix to `indNames(x)` /
   `popNames(x)` by name when dimnames are present, otherwise use object
   order (F3). **Consequence: links drawn from individual-level matrices
   move to the individuals the matrix actually names; population-level
   matrices with named rows are reordered to `popNames(x)`.**
4. Draw a symmetric link only when `i != ii`, the cell is not `NA`, and
   the value is above 0 (F4). **Consequence: self-links and zero-valued
   pairs are no longer drawn; with `standard = TRUE` (default) zero
   distances still draw because they rescale to 1.**
5. Call `gl.colors("div", verbose = 0)` (F5).
6. Test `is.na()` on both cells in the asymmetric branch; skip or grey
   missing pairs instead of crashing (F6).
7. Input and dependency guards: error when `ind.circle.cols` has fewer
   colours than populations; coerce `latlon` with `as.data.frame()`;
   `stop(error(...))` for `leaflet`/`leaflet.minicharts`; guard `terra`
   and `scales` when `raster.image` is given (F7, F8, F9). Callers that
   tested for a `-1` return now receive an error.
8. Documentation pass: `@description`, `@return`, DOC2 verbose text,
   `standard`/`palette.links`/`legend.title`/`matrix` param texts, typos,
   Custodian label; then `devtools::document()` (F10, F11, F12).
9. Optional: also map the (standardised) value to `weight` in the
   symmetric branch so line width follows the value as documented (F10).
   **Consequence: symmetric links change width on screen for every
   existing call with a matrix.**
10. Housekeeping: drop `build`, drop the redundant `addTiles()`, print
    "Completed:" after the raster step (F13).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: behaviour vs roxygen on platypus.gl, bandicoot.gl, testset.gs --
  run (default map, switches, pop-level matrix, ind-level matrix,
  symmetric and asymmetric, `standard` on/off, NA cells, short colour
  vector, matrix `latlon`, missing/misnamed `latlon`)
- Single-population path: run on platypus.gl (lat,lon) and bandicoot.gl
  (lon,lat)
- FBM path (DAT6): run on `gl.gen2fbm(platypus.gl)`; the function reads
  no genotypes, so FBM is safe by construction
- Raster path: run with a synthetic EPSG:4326 GeoTIFF; works when `terra`
  is installed
- Read-only check: input identical after the call -- run
- Rendered output (tile downloads, browser rendering): SKIPPED -- not
  reproducible in a test; introspection stops at the widget call list
- Google Group / GitHub issues: searched `gl.map.interactive` in
  dartR.base and dartR issue trackers -- no open or closed issue found
- dartr2shiny: `input_generator/dartR.base/gl.map.interactive.r` is a
  copy of this file's header; `config/functions.csv` lists the function
  under Data_Overview with a leaflet output slot. No signature change is
  proposed, so the generator is unaffected.

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
| 10 | approved | Luis | |

## Outcome

Applied on branch `review-gl.map.interactive` (from origin/dev be7cb2e);
PR #403.

- Change 1: `centers` built from `tapply` per column with `popNames`
  row names; single-population label now at the centre for both column
  orders (test "F1", both platypus.gl and bandicoot.gl).
- Change 2: `dist` coerced with `as.matrix`; non-square or non-matrix
  input stops with "must be a square matrix or a dist object" (test
  "F2": gl.dist.pop and gl.dist.ind inputs draw 3 and 15 links).
- Change 3: row-only reorder removed; alignment by dimnames with a
  `verbose >= 2` warning when names are present but do not match (test
  "F3": named, permuted, unnamed and mismatched matrices all draw the
  1-2 link; a reversed population matrix draws the same links).
- Change 4: symmetric loop draws only `i != ii && !is.na(v) && v > 0`
  (3 pops: 3 links instead of 6; 6 individuals with one cell: 1 link
  instead of 21).
- Change 5: `gl.colors("div", verbose = 0)` (0 output lines at
  `verbose = 0` with a matrix).
- Change 6: asymmetric loop tests `is.na` on both cells; one NA cell
  gives 29 flows and one grey arrow instead of a crash.
- Change 7: short colour vector stops with "colours but the dataset
  has"; `latlon` coerced with `as.data.frame`; DEP1 `stop(error())` for
  leaflet/leaflet.minicharts, plus terra/scales when `raster.image` is
  given.
- Change 8: roxygen updated (`@description`, `@return`, DOC2 verbose
  text, `standard`/`palette.links`/`legend.title`/`matrix` texts, typos,
  `Author(s)`/`Custodian`); `devtools::document()` run, only
  `man/gl.map.interactive.Rd` changed; `devtools::check_man()` clean.
- Change 9: `weight = v` in `addPolylines`; widths equal the
  standardised values (default) or the raw values (`standard = FALSE`).
- Change 10: `build` argument dropped, `addTiles()` dropped (call list
  starts with `addProviderTiles`), "Completed:" printed after the raster
  step.
- Characterization test: 46 baseline assertions at be7cb2e; after the
  changes 57 assertions pass, every flipped assertion tagged
  `[approved diff, change n]` in the file.
- End-to-end at `verbose = 3`: platypus.gl with `gl.dist.pop` (symmetric)
  and a 6-individual asymmetric `gl.dist.ind`; both examples; raster
  branch with a synthetic GeoTIFF; widget saved with `htmlwidgets::saveWidget`.
- NEWS.md entry added under 1.2.3 (development).
- Callers: no sibling `dartR.*` package calls the function; dartr2shiny
  keeps a header copy under `input_generator/` and a leaflet output slot
  in `config/functions.csv`; the signature is unchanged.

```json
{
  "function": "gl.map.interactive",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "be7cb2e",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "STY3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "DAT5", "status": "approved", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 7},
    {"id": "F10", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "approved", "change": 8},
    {"id": "F13", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "approved", "change": 10}
  ],
  "coverage_skipped": ["rendered output: not reproducible in a test"],
  "status": "pr-open",
  "pr": 403
}
```
