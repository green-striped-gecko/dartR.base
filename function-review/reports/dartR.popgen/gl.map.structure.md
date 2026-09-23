# Review: gl.map.structure (dartR.popgen)
- Family mode: analysis (map front end over `gl.plot.structure` output)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 95fde36 (origin/dev)
- Datasets: testset.gl (three populations, 31 individuals, `latlon` stored
  as columns `lat`, `lon`; all 250 individuals for timing); hand-built
  q-matrix tables in the layout `gl.plot.structure` returns
- Baseline: tests/testthat/test-gl.map.structure.R (snapshot captured
  pre-review; needs `leaflet`, not STRUCTURE)

## Verdict

**Standards: Needs work** — the function has none of the house anatomy (no
`verbose`, no start/end flags, no argument checks, no `leaflet` guard), and
its example calls `gl.plot.structure` with an argument that does not exist.
**Spec: Needs work** — when the population levels of `x` are not in
alphabetical order, each population's bars are drawn at another
population's centre under the wrong label; `movepops` shifts by column
position, so it moves latitude on genlights that store `lat` first.

What works well: with alphabetical population levels and a `latlon` in
`lon`, `lat` order, bars, labels and the returned per-population tables are
correct.

## Findings

**F1 [HIGH, confidence: high] — population centres matched in the wrong direction (DOC5)**
`R/gl.map.structure.r:123-132` — `sc <- match(rownames(centers),
levels(factor(qmat$orig.pop)))` gives, for each genlight population, its
position in the q-matrix levels; `centers[sc, ]` then applies that as a
reordering, which is the inverse permutation. Bars are drawn per q-matrix
level (`levels(bb$orig.pop)[p]`) at `cx[p]`, while labels use
`rownames(centers)`.
Failure scenario: `pop(x)` levels `EmmacBurnBara, EmmacBrisWive,
EmmacBurdMist` (not alphabetical): EmmacBrisWive's bars are drawn at
EmmacBurdMist's centre, EmmacBurdMist's at EmmacBurnBara's, and each label
sits at its own coordinates, so the map shows the wrong ancestry for every
population without any message. With alphabetical levels the permutation is
the identity and the map is right, which hides the defect.
Proposed change: order centres by the q-matrix populations,
`centers[match(levels(bb$orig.pop), rownames(centers)), ]`, and use the
same order for bars and labels.

**F2 [MEDIUM, confidence: high] — population mismatch is reported with `message()` and then fails (FS5, VRB2)**
`R/gl.map.structure.r:125-131` — a mismatch prints a red message and
continues.
Failure scenario: STRUCTURE was run on three populations and `x` has four
(a common case: `x` is the full dataset). The call prints the message, then
120+ leaflet warnings ("missing or invalid lat/lon values"), then stops with
"row names contain missing values". The same happens when `x` lacks a
population that is in `qmat`.
Proposed change: compute centres only for the populations in `qmat`
(populations of `x` that are not in `qmat` are ignored, with a note at
`verbose >= 2`); stop with `error()` naming any `qmat` population that has
no coordinates in `x`.

**F3 [MEDIUM, confidence: high] — `movepops` is added by column position (DOC5)**
`R/gl.map.structure.r:118-119` — `centers[, 1] + movepops[, 1]` assumes the
first column of `x@other$latlon` is longitude. testset.gl (and any genlight
built with `lat` first) stores `lat`, `lon`.
Failure scenario: `movepops = data.frame(lon = c(1, 0, 0), lat = 0)` on
testset.gl moves the first population 1 degree north instead of east.
Proposed change: add `movepops$lon` to the `lon` column and `movepops$lat`
to the `lat` column by name (unnamed columns read as lon, lat, as
documented); match rows to populations by row name when they are
population names, otherwise in the order of `levels(pop(x))` as today.
Document the row order.

**F4 [MEDIUM, confidence: high] — inputs are not checked (FS5)**
`R/gl.map.structure.r:89-107` — nothing checks `x` (class, `latlon`
present, `lon`/`lat` columns, coordinates present) or `qmat` (a list of
tables with `cluster*` and `orig.pop` columns).
Failure scenario: a genlight whose `latlon` is all NA returns a map with no
drawable bars and 60+ leaflet warnings, and no error.
Proposed change: check up front with `stop(error(...))`: `x` is a genlight
(`utils.check.datatype`), `x@other$latlon` has `lon` and `lat` columns, no
mapped population has only missing coordinates, `qmat` has the columns
`gl.plot.structure` returns, `K` is a single value.

**F5 [LOW, confidence: high] — one population or K = 1 fails (principle: output integrity)**
`R/gl.map.structure.r:105-108, 146, 135` — `apply()` over `tapply()` returns
a vector, not a matrix, when there is one population, and one cluster
column drops to a vector; with one population the bar width
`range(lon) / 100` is 0.
Failure scenario: `gl.map.structure(q, x, K = 1)` and a single-population
`x` both stop with "incorrect number of dimensions".
Proposed change: keep matrix shape (`drop = FALSE`); when all centres share
one longitude, use a bar width of 0.01 degrees times `scalex`.

**F6 [LOW, confidence: high] — map colours differ from the bar plot (PLT1)**
`R/gl.map.structure.r:172` — clusters are coloured with `rainbow(K)`, while
`gl.plot.structure` uses `gl.select.colors()`, so cluster 1 is red on the
map and a different colour in the bar plot.
Failure scenario: a user reads the map legend from the bar plot and assigns
the wrong cluster to each colour.
Proposed change: add an optional `plot.colors` argument (vector or palette
function) defaulting to `gl.select.colors(ncolors = K)`, the
`gl.plot.structure` default. The two match exactly when `gl.plot.structure`
was called for that K alone (it sizes its palette to the largest K).

**F7 [LOW, confidence: high] — only the first mode of a K can be mapped (DOC5)**
`R/gl.map.structure.r:89-99` — `qmat[eq.k][[1]]` silently takes the first
table with K clusters.
Failure scenario: `gl.plot.structure` returned modes `2.1` and `2.2`;
`K = 2` maps `2.1`, and there is no way to map `2.2`.
Proposed change: accept `K` as a panel label (`"2.2"`) as well as a number;
when a number matches several modes, map the first and say so at
`verbose >= 2`.

**F8 [LOW, confidence: high] — no verbose argument, flags, or plot switch (FS2, FS3, FS9, PLT1)**
`R/gl.map.structure.r:80-88, 194` — no `verbose`, no
`utils.flag.start`/end, and the map is printed on every call.
Failure scenario: `gl.set.verbosity(0)` does not silence the function's
messages; a script that only wants the returned tables cannot skip drawing.
Proposed change: add `plot.out = TRUE` and `verbose = NULL` at the end of
the signature, with the standard start/end flags.

**F9 [LOW, confidence: high] — `leaflet` from Suggests used without a guard (DEP1)**
`R/gl.map.structure.r:149` — no `requireNamespace("leaflet")`.
Failure scenario: without `leaflet` installed, the call stops with "there
is no package called 'leaflet'".
Proposed change: add the DEP1 guard.

**F10 [LOW, confidence: high] — documentation errors (DOC1, DOC5, DOC7)**
`R/gl.map.structure.r:1-78` —
- The example calls `gl.plot.structure(sr, k = 2:4)`: there is no `k`
  argument, so the example stops with "unused argument"; the line ends in a
  stray `#' #head(qmat)`.
- Two `@return` blocks (the second is correct: a list with `qmats` and
  `map`).
- `movepops` does not say which row belongs to which population, or that
  columns are read by position.
- `K` does not say that the table is chosen by its number of cluster
  columns.
- Typos ("plotstructure", "parplots", "corret"); `@family` missing;
  `@author` lacks the Author(s)/Custodian structure (DOC7, proposed rule).
Failure scenario: a user who copies the example gets an error at the
`gl.plot.structure` line.
Proposed change: fix each item; regenerate the Rd.

**F11 [INFO, confidence: high] — one leaflet layer per bar segment (STY2)**
`R/gl.map.structure.r:161-176` — `addRectangles` is called once per
individual and cluster.
Failure scenario: 250 individuals at K = 3 produce 752 widget calls, a
3.5 MB widget, 0.8 s; this grows linearly and is acceptable at typical
sizes.
Proposed change: none now; vectorising per population would give the same
rectangles if map size becomes a complaint.

## Proposed changes

1. Order population centres by the q-matrix populations so bars, labels and
   coordinates agree (F1). **Consequence: the map changes whenever the
   population levels of `x` are not alphabetical — bars move to their own
   population's centre.**
2. Centres only for populations in `qmat`, extra populations of `x`
   ignored with a note, missing coordinates a clear error (F2).
   **Consequence: calls with extra populations in `x`, which errored, now
   return a map.**
3. `movepops` added by column name, rows matched by population name when
   named (F3). **Consequence: `movepops` shifts change for genlights whose
   `latlon` stores `lat` before `lon`.**
4. Up-front checks of `x`, `latlon`, `qmat` and `K` with `stop(error())`
   (F4).
5. One population and K = 1 work; fallback bar width (F5).
6. New `plot.colors` argument with the `gl.plot.structure` default palette
   (F6). **Consequence: default map colours change from `rainbow()` to
   the dartR palette.**
7. `K` accepts a mode label such as `"2.2"`; note when several modes match
   (F7).
8. `plot.out` and `verbose` arguments, start/end flags (F8).
9. DEP1 guard for `leaflet` (F9).
10. Documentation corrections in F10, regenerate Rd, NEWS entry.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run.
- Spec: behaviour vs roxygen on testset.gl with hand-built q-matrix tables —
  run; bar positions and labels read from the leaflet widget's call list.
- Rendered map appearance: SKIPPED as a snapshot (htmlwidget; no
  webshot/vdiffr in Suggests).
- Real STRUCTURE output: not used; `gl.plot.structure` output was
  reproduced by hand in its documented layout (Label, cluster1..K, K,
  orig.pop, ord).
- DAT1–DAT6: not applicable — `x` is read, not modified.
- dartR Google Group / GitHub issues search: not run (no search access in
  this session).
- Downstream callers: dartr2shiny (`input_Tabs.csv`, `variables_matrix.csv`,
  `template_report.csv`) — grepped; new arguments are added at the end of
  the signature, so positional calls are unaffected. No callers in other
  `dartR.*` packages.
- Interaction with PR #94 (`gl.plot.structure`): with `den = TRUE` the
  current `gl.plot.structure` blanks `orig.pop`, which makes this function
  fail at F2; PR #94 fixes that on the producer side.

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

- Changes 1–10 applied in commit 8653e41 on `review-gl.map.structure`, PR #95 to `dev`.
- Characterization test: 41 expectations pass; every diff from baseline is
  tagged `[approved n]` (1–8). Old vs new where no fix applies (alphabetical
  levels; lon/lat order with `movepops`; all 30 testset.gl populations):
  identical `qmats`, rectangles and labels.
- `verbose = 3` end to end on four populations of testset.gl with a
  three-population `qmat`: extra population reported, three maps drawn.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; no new NOTE.
- Not checked: rendered map in a browser.

## Machine block

## Addendum (2026-09-23)

Found during the gl.read.structure review: a q-matrix with `orig.pop` all NA
stopped with "attempt to select less than one element in integerOneIndex",
and partly NA added empty rows and leaflet warnings. Now: all NA is a clear
error, partial NA individuals are dropped with a warning. Commit 0215099,
PR #99. Requested by Luis.

```json
{
  "function": "gl.map.structure",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "95fde36",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "principle: output integrity", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 10},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "STY2", "status": "no-change", "change": null}
  ],
  "coverage_skipped": [
    "rendered map snapshot: htmlwidget, no webshot/vdiffr",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 95
}
```
