# Review: gl.plot.structure (dartR.popgen)
- Family mode: analysis (plotting front end over CLUMPP/Clumpak averaging)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 95fde36 (origin/dev)
- Datasets: testset.gl (three populations, 31 individuals); a real
  `structure.result` from STRUCTURE 2.3.4 (`~/programs/structure`, K = 1:4,
  five replicates each) for exploration; hand-built `structure.result`
  objects (one mode per K, and two distinct modes at K = 2) for the baseline
- Baseline: tests/testthat/test-gl.plot.structure.R (snapshot captured
  pre-review; runs without the STRUCTURE binary)

## Verdict

**Standards: Needs work** — the structure follows the house anatomy, but
arguments are not checked before work starts, two Suggests packages are
used unguarded, and the documentation disagrees with the code in several
places.
**Spec: Needs work** — when Clumpak finds more than one mode, each mode is
represented by a single replicate instead of the documented mode average,
and three smaller output defects (duplicated K = 1 results, population
labels lost under `den = TRUE`, a dendrogram built on distances of
distances) change what the user sees.

What works well: with a single mode per K (the common case), the CLUMPP
alignment, averaging, per-population sorting and the shared individual
order across K panels are correct.

## Findings

**F1 [HIGH, confidence: high] — mode averaging uses only the first replicate (DOC5)**
`R/gl.plot.structure.r:212-220` — inside the multi-mode branch, `x[1]` is a
one-element sub-list, so `length(x[1])` is always 1 and every mode returns
`x[[1]]`, the first replicate of that mode. The `@details` promise that the
function "averages the replicates within each mode" holds only when there
is exactly one mode.
Failure scenario: six K = 2 replicates falling into two modes of three
(hand-built `structure.result`): both returned q-matrices are identical to
a single replicate (max difference 0), not to the mean of their mode.
Proposed change: average all replicates of each mode
(`Reduce("+", x) / length(x)`).

**F2 [MEDIUM, confidence: high] — K = 1 duplicates earlier results (principle: output integrity)**
`R/gl.plot.structure.r:169` — the K = 1 branch appends
`c(res, as.matrix(Q_list_tmp[1]))`, i.e. the whole result list built so far
plus the K = 1 matrix. When K = 1 is not the first K processed, the earlier
Ks appear twice.
Failure scenario: `K = c(2, 1)` returns three q-matrices labelled
`2.1`, `2.2`, `1.1` — the K = 2 panel is plotted twice and relabelled as if
it were two modes. `K = NULL` and `K = 1:n` are unaffected because K = 1
comes first.
Proposed change: append only the K = 1 matrix.

**F3 [MEDIUM, confidence: high] — `den = TRUE` blanks `orig.pop` in the returned q-matrices (PLT3)**
`R/gl.plot.structure.r:263-265` — the population column is overwritten
with `" "` in the tables that are returned, not only in the plotting copy.
Failure scenario: `q <- gl.plot.structure(sr, K = 2, den = TRUE, x = x)`
returns `orig.pop == " "` for every individual; passing `q` to
`gl.map.structure()` then fails its population-name match against the
genlight.
Proposed change: keep `orig.pop` in the returned tables; blank it only in
the data used to draw the plot.

**F4 [MEDIUM, confidence: medium] — dendrogram is built on distances between rows of the distance matrix (DOC5)**
`R/gl.plot.structure.r:319` — `dist(res)` is applied to an object that is
already a `dist` (from `dis.mat` or `gl.dist.ind`), so `hclust` clusters
Euclidean distances between the rows of the distance matrix rather than the
distances themselves. The documentation says `dis.mat` is "used to order"
the plot.
Failure scenario: on testset.gl (three populations), the leaf order from
`hclust(dist(d))` differs from `hclust(d)` at 28 of 31 positions; the
cophenetic correlation with the supplied distances drops from 0.65 to 0.59.
A user who supplies their own `dis.mat` gets a tree that is not the
clustering of that matrix.
Proposed change: call `hclust(as.dist(res))` directly.

**F5 [MEDIUM, confidence: high] — arguments are not checked before work starts (FS5)**
`R/gl.plot.structure.r:133-160, 310-326` — only the class of `sr` is
checked. `den = TRUE` needs `x` even when `dis.mat` is supplied (labels are
taken from `indNames(x)`), `met_clumpp` is validated deep inside `clumpp()`
with an uncoloured `stop()`, and a short `color_clusters` fails in ggplot2.
Failure scenario: `den = TRUE` without `x` stops with "inappropriate object
passed to function, found NULL expecting SNP or SilicoDArT"; with `dis.mat`
but no `x` it stops with "unable to find an inherited method for function
'indNames' for signature x = NULL". Neither names the missing argument.
Proposed change: check up front — `den = TRUE` needs `dis.mat` or `x`
(take labels from `labels(dis.mat)` when supplied, and check they match the
individuals in `sr`); `met_clumpp` in the three allowed values; every
requested K present in `sr`; enough colours for the largest K. Each check
fails with `stop(error(...))` naming the argument.

**F6 [LOW, confidence: high] — a palette function is documented but rejected (DOC5)**
`R/gl.plot.structure.r:294-296, 356` — `@param color_clusters` says "A color
palette for clusters (K) or a list with as many colors", but the value goes
straight to `scale_fill_manual(values = )`.
Failure scenario: `color_clusters = rainbow` stops with "Insufficient values
in manual scale. 2 needed but only 1 provided."
Proposed change: if `color_clusters` is a function, call it with the
largest K.

**F7 [LOW, confidence: high] — Suggests packages used without a guard (DEP1)**
`R/gl.plot.structure.r:191, 304` — `proxy::simil` and `reshape2::melt` come
from Suggests with no `requireNamespace()` check.
Failure scenario: on an install without `proxy`, any call with more than
one replicate per K stops with "there is no package called 'proxy'"
instead of the dartR install message.
Proposed change: add the DEP1 guard for `proxy` and `reshape2` at the top.

**F8 [LOW, confidence: high] — deprecated `aes_()` (PLT1)**
`R/gl.plot.structure.r:347` — ggplot2 has deprecated `aes_()` since 3.0.0.
Failure scenario: the first call in a session prints a lifecycle warning
asking the user to report the issue to the dartR group; a future ggplot2
release that removes `aes_()` breaks the function.
Proposed change: use `aes(x = factor(.data$ord), y = .data$value,
fill = .data$Cluster)`.

**F9 [LOW, confidence: high] — documentation disagrees with the code (DOC1, DOC2, DOC5, DOC7)**
`R/gl.plot.structure.r:1-98` —
- `border_ind` documented default 0.25; code default 0.15.
- `plot.dir` documented "[default = working directory]"; `gl.check.wd()`
  resolves NULL to `tempdir()`.
- `@return` says "List of Q-matrices"; the function returns a list of
  `data.table`s with columns `Label`, `cluster1..K`, `K`, `orig.pop`, `ord`,
  rows sorted by `Label`, names `"1"`, `"2"`, ...
- `k_name` does not say its values are the panel labels (`"3"`, or `"2.1"`
  when a K has several modes).
- `x` text is garbled ("used in gl.run.structure description").
- `verbose` text is not the standard DOC2 wording; `@family` is missing;
  `@author` lacks the Author(s)/Custodian structure (DOC7, proposed rule);
  `@seealso` lists the function itself; the CLUMPP permutation search is
  random and this is not mentioned (colour order can differ between calls
  without `set.seed()`).
Failure scenario: a user who reads the help page expects 0.25-width borders
and files in the working directory, and cannot find `k_name` values.
Proposed change: correct each item above; regenerate the Rd.

**F10 [LOW, confidence: high] — outdated start flag and a misleading error message (FS3, VRB2)**
`R/gl.plot.structure.r:125-129, 156-158` — `utils.flag.start()` is called
with the obsolete `build = "Jody"`; the missing-K message pastes the whole
`K` vector.
Failure scenario: `K = c(2, 7)` reports "No entries for K = 2 found in
'sr'. No entries for K = 7 found in 'sr'." although K = 2 exists.
Proposed change: drop `build =`; name only the missing K values (folded
into the up-front check of change 5).

**F11 [INFO, confidence: high] — argument names follow the old underscore style (PLT1)**
`R/gl.plot.structure.r:100-116` — `plot_theme`, `color_clusters`,
`met_clumpp`, `iter_clumpp`, `ind_name`, `k_name`, `border_ind` versus the
current `plot.theme`/`plot.colors` idiom. The same names are shared by
`gl.plot.faststructure`, `gl.plot.snmf` and `gl.plot.popcluster`, and
dartr2shiny calls this function.
Failure scenario: none at run time; inconsistent API across the suite.
Proposed change: none here — a rename is an API change (API2, API3) best
made once for the whole structure/admixture plotting family.

## Proposed changes

1. Average every replicate within each Clumpak mode instead of returning the
   first one (F1). **Consequence: numerical output (returned q-matrices and
   bar heights) changes whenever Clumpak finds more than one mode at a K.**
2. Append only the K = 1 q-matrix in the K = 1 branch (F2).
   **Consequence: when K = 1 is requested after another K (e.g.
   `K = c(2, 1)`), the returned list and the plot lose the duplicated
   panels, and panel labels change from `2.1`/`2.2` to `2`.**
3. Keep `orig.pop` in the returned tables under `den = TRUE`; blank it only
   for plotting (F3). **Consequence: the returned object under
   `den = TRUE` now carries population names; the plot is unchanged.**
4. Build the dendrogram with `hclust(as.dist(res))` instead of
   `hclust(dist(res))` (F4). **Consequence: individual order in
   `den = TRUE` plots changes.**
5. Up-front argument checks with `stop(error())`: `den = TRUE` needs
   `dis.mat` or `x`, labels taken from `dis.mat` when supplied and checked
   against `sr`; `met_clumpp` valid; requested K present (message names only
   the missing K); enough colours (F5, F10 message part).
6. Accept a palette function in `color_clusters` (F6).
7. DEP1 guards for `proxy` and `reshape2` (F7).
8. Replace `aes_()` with `aes()` and `.data` (F8).
9. Documentation corrections listed in F9, drop `build =` from
   `utils.flag.start()` (F10), regenerate Rd, NEWS entry.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run.
- Spec: behaviour vs roxygen on a real STRUCTURE run (K = 1:4, 5 replicates)
  and hand-built `structure.result` objects — run.
- Multi-mode averaging: tested on a hand-built two-mode fixture; the real
  STRUCTURE run on testset.gl converged to one mode at every K, so it could
  not exercise F1.
- Plot appearance: SKIPPED as a snapshot (pure ggplot/patchwork output, no
  vdiffr in Suggests); checked by rendering to a null device without error.
- `aes_()` deprecation warning: asserted from the source, because ggplot2
  emits it once per session and a test cannot trigger it reliably.
- DAT1–DAT6 (genlight integrity): not applicable — the function does not
  modify a genlight.
- dartR Google Group / GitHub issues search for known complaints: not run
  (no search access in this session).
- Downstream callers: `gl.map.structure` (consumes the returned list, selects
  by `ncol - 4 == K` and matches `orig.pop`), dartr2shiny (`input_Tabs.csv`,
  `template_report.csv`) — grepped; no caller passes `K` with 1 out of first
  position.

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

- Changes 1–9 applied in commit a749667 on `review-gl.plot.structure`, PR #94 to `dev`.
- Characterization test: 35 expectations pass; every diff from baseline is
  tagged `[approved n]` (1, 2, 3, 4, 5, 6, 8). Old vs new on the real
  STRUCTURE run: identical returned objects for all unaffected calls;
  identical bars (only ggplot's internal group numbering changed with the
  aesthetic order).
- Change 3 side effect: under `den = TRUE` the returned tables are now the
  same as with `den = FALSE` (per-population sort in `ord`); the plot order
  still comes from the dendrogram.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; install
  warning is local packages built under R 4.4.3.

## Machine block

```json
{
  "function": "gl.plot.structure",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "95fde36",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "principle: output integrity", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 9},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "PLT1", "status": "no-change", "change": null}
  ],
  "coverage_skipped": [
    "plot snapshot: no vdiffr",
    "real multi-mode STRUCTURE run: fixture converged to one mode",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 94
}
```
