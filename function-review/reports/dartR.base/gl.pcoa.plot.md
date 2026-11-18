# Review: gl.pcoa.plot (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5), via dartr-dev agent
- Skill: dartr-function-review v2.0.0 (Phase A only — read-only review)
- Package commit: ddaed27 (upstream/dev; `git diff upstream/dev --
  R/gl.pcoa.plot.r` empty on the working tree, so `load_all()` exercised
  the reviewed code)
- Branch reviewed from: integration-local (ed99203; file identical to
  upstream/dev)
- Date: 2026-09-07
- Family mode: plotting/report — the product is a plot of an ordination;
  the input objects must come back untouched
- Datasets: testset.gl, testset.gs, testset.gl[1:30,] +
  gl.dist.ind (individual distance, 9 negative eigenvalues),
  gl.dist.pop(testset.gl) (pop-level distance), gl.pcoa fixtures at
  nfactors 1/2/3/5, correction none/cailliez, shuffled-pop fixture,
  as.pop = "sex"
- Baseline: tests/testthat/test-gl.pcoa.plot.R (new file; 67 assertions,
  all passing at ddaed27 — defects pinned as-is and tagged with finding
  IDs). Introspection is on the ggplot/plotly object (labels, layer data,
  mappings), not rendered images.

## Verdicts

**Standards: Needs work** — the entry scaffold (verbosity, flag start,
datatype gates, flag end) conforms, but a dozen messages print at
`verbose = 0`, the dependency guards return `-1` instead of stopping, and
the save path skips `gl.check.wd` so `plot.file` writes into the working
directory.

**Spec: Needs work** — axis truth and data truth both hold (the plot's
strongest properties), but `hadjust`/`vadjust` are silently ignored,
`pop.labels = "ind"` crashes, the axis-bounds "clamp" writes out-of-bounds
values that crash on small ordinations, and `@return` documents NULL while
a visible ggplot is returned.

What works well: every plotted coordinate is exactly the score column it
claims to be, on SNP, SilicoDArT, dist and 3D branches alike, and every
axis label's percentage equals `100*eig/sum(eig >= 0)` for the object's own
eigenvalues.

## Axis truth and data truth (spec axis, empirical)

- **Axis truth: PASS.** For all tested `xaxis`/`yaxis`/`zaxis` choices
  (1/2, 2/5, z=3) the label percentage equals
  `round(100 * eig[k] / sum(eig[eig >= 0]), 1)` computed from the input
  object's `$eig`. The denominator is the positive eigenvalues only, which
  matches `gl.pcoa`'s own dist-branch convention (`eig.raw.pos.pc`) and is
  the honest choice when negative eigenvalues exist; on a 30-entity
  distance ordination with 9 negative eigenvalues the labels remained
  consistent with that convention. The SNP branch of `gl.pcoa` reports
  percentages over `sum(pca$eig)` (all eigenvalues), but glPca eigenvalues
  are non-negative there, so the two conventions coincide.
- **Data truth: PASS.** `p$data$PCoAx/PCoAy` are identical to
  `glPca$scores[, xaxis/yaxis]` (2D, all label modes), the 3D plotly
  traces carry the score columns exactly, `ind` equals `indNames(x)`, and
  `pop` tracks `pop(x)` 1:1 — verified with a shuffled-pop fixture and
  with `as.pop = "sex"`.

## Findings

**F1 [HIGH, confidence: high] — `pop.labels = "ind"` accepted, then
crashes (DOC5, FS5)**
`R/gl.pcoa.plot.r:247-256` (gate) vs `:496-804` (branches) — the
validation accepts `none|ind|pop|legend`, but no `"ind"` plotting branch
exists; `plott` is never assigned and `show(plott)` fails.
Failure scenario: `gl.pcoa.plot(pca, gl, pop.labels = "ind")` errors with
"object 'plott' not found" (reproduced). The `@details` text promises
specimens "shown optionally with adjacent labels"; the `@param` list omits
`"ind"`.
Proposed change: either implement an individual-labels branch (the
interactive branch already maps `label = ind`) or remove `"ind"` from the
gate and route it to the fallback warning; align `@param`/`@details`.

**F2 [HIGH, confidence: high] — `hadjust` and `vadjust` are silently
ignored (DOC5)**
`R/gl.pcoa.plot.r:138-139, 265-277` — both parameters are documented
("Horizontal/Vertical adjustment of label position"), validated, then
never used; no later code references them (grep: only signature,
validation, examples).
Failure scenario: a user tuning label positions gets byte-identical plots
for `hadjust = 0` vs `3` (reproduced: identical `ggplot_build()` data).
Proposed change: wire them into the label placement or remove the
parameters and their docs (removal is API2 — see change 2).

**F3 [HIGH, confidence: high] — axis-bounds "clamp" writes out-of-bounds
values and crashes (FS5)**
`R/gl.pcoa.plot.r:290-312` — an out-of-range `yaxis` is reset to the
constant 2 and `zaxis` to 3, without checking those against
`ncol(glPca$scores)`.
Failure scenario: any 1-axis ordination — `gl.pcoa(gl, nfactors = 1)`, or
a 2-individual genlight — dies with "subscript out of bounds" on the
default call; `zaxis = 5` on a 2-factor ordination likewise (all three
reproduced). PR #369 makes small ordinations more common by clamping
`nfactors` on the FBM and dist paths instead of crashing in `gl.pcoa`, so
the crash now surfaces here.
Proposed change: clamp to valid columns (`min(., ncol(scores))`), require
distinct axes, and error informatively when fewer than 2 (or 3, for
`zaxis`) axes exist.

**F4 [MEDIUM, confidence: high] — `plot.file` with default `plot.dir`
writes to the working directory, not the documented tempdir (FS7, PLT2)**
`R/gl.pcoa.plot.r:882-887` — `plot.dir` is passed straight to
`utils.plot.save`, whose `dir = NULL` default is `getwd()`; the function
never calls `gl.check.wd`. The roxygen (`:40-41`) promises "otherwise to
the tempdir()".
Failure scenario: `gl.pcoa.plot(pca, gl, plot.file = "x")` drops `x.RDS`
into the user's current directory (reproduced).
Proposed change: add `plot.dir <- gl.check.wd(plot.dir, verbose = 0)` to
the preamble, per the house idiom.

**F5 [MEDIUM, confidence: high] — corrected PCoA is labelled "PCA Axis"
(DOC5; axis-label class)**
`R/gl.pcoa.plot.r:170-176` — the PCA/PCoA classification keys on
`glPca$loadings` being NULL. `gl.pcoa` on a distance matrix with
`correction = "cailliez"`/`"lingoes"` stores `vectors.cor` in `$loadings`,
so the object classifies as "PCA".
Failure scenario: `gl.pcoa.plot(gl.pcoa(D, correction = "cailliez"), gl)`
titles the axes "PCA Axis 1 (...)" for a PCoA (reproduced). With
`correction = "none"` the labels are correctly "PCoA Axis".
Proposed change: classify from a positive signal (e.g. the `$call`, or a
flag set by `gl.pcoa`) rather than the absence of loadings; minimally,
also treat objects whose call used a dist as PCoA.

**F6 [MEDIUM, confidence: high] — `verbose = 0` is not silent (VRB3,
VRB5)**
Multiple sites — all six parameter-validation warnings
(`:250-311`) are ungated; the interactive branch prints three ungated
lines (`:585-590, 648-650`); the "none" branch's message is gated at
`verbose >= 0`, which is always true (`:741`).
Failure scenario: `verbose = 0` with `pop.labels = "none"`, any invalid
parameter, or `interactive = TRUE` produces console output (all
reproduced; the default call is clean — pinned at 0 lines).
Proposed change: gate the warnings at `verbose >= 2` (VRB3) and fix the
`>= 0` typo to `>= 2`.

**F7 [MEDIUM, confidence: high] — the plot always renders; no
`plot.display` gate (VRB5, PLT idiom)**
`R/gl.pcoa.plot.r:809-814, 854` — `show(plott)` runs unconditionally on
every path; the function has no `plot.display` parameter, so `verbose = 0`
still pops a graphics window/viewer.
Failure scenario: batch scripts calling the function for the returned
object (or the RDS side effect) cannot suppress the display.
Proposed change: add the standard `plot.display` gate (`if (verbose == 0)
plot.display <- FALSE`), keeping the returned object independent of
display (PLT3 already holds).

**F8 [MEDIUM, confidence: high] — `@return` documents NULL; a visible
ggplot/plotly is returned (DOC5, FS10)**
`R/gl.pcoa.plot.r:84` vs `:898` — "@return returns no value (i.e. NULL)"
but the function ends `return(plott)`; the object prints on call
(reproduced: `withVisible()$visible` is TRUE), so an unassigned call
renders the plot twice.
Proposed change: document the returned plot object; returning
`invisible(plott)` would match the family convention but is a visible
behaviour change (see change 8's consequence line).

**F9 [MEDIUM, confidence: high] — dependency guards return -1 instead of
stopping (DEP1)**
`R/gl.pcoa.plot.r:196-240` — the directlabels/plotly/gganimate/tibble
guards use `cat(error(...)); return(-1)` rather than `stop(error(...))`.
Failure scenario: on a machine without directlabels,
`p <- gl.pcoa.plot(...)` yields `p == -1` and downstream code fails
obscurely instead of the call halting.
Proposed change: use the DEP1 idiom `stop(error(...))`.

**F10 [LOW, confidence: high] — `vadjust` range check tests `hadjust`
(FS5)**
`R/gl.pcoa.plot.r:272` — `if (vadjust < 0 | hadjust > 3)`; a too-large
`vadjust` (e.g. 5) passes silently (reproduced), and when the block does
fire it resets `vadjust` to 1.5, not the documented default 1.
Proposed change: test `vadjust > 3` and reset to 1 — moot if F2 removes
the parameters.

**F11 [MEDIUM, confidence: high] — input pairing unguarded: opaque crashes
for pop-level dist ordinations and data.frame `glPca` (DAT5, FS5)**
`R/gl.pcoa.plot.r:164-240, 442-458` — nothing checks
`nrow(glPca$scores) == nInd(x)`, and the `accept = c("glPca","list")` gate
admits any list-classified object (utils.check.datatype classifies
data.frame as "list" — its review, F8) and does not require `glPca` and
`x` to be lists together.
Failure scenario: `gl.pcoa.plot(gl.pcoa(gl.dist.pop(gl)), gl)` dies with
"arguments imply differing number of rows: 30, 250"; a data.frame passed
as `glPca` routes into the animation branch and dies with "this S4 class
is not subsettable" (both reproduced). Genuinely wrong classes (matrix,
dist) are rejected informatively by the gate.
Proposed change: fail fast with clear messages — require
`nrow(glPca$scores) == nInd(x)` (pointing pop-level ordinations at the
commented-out entity-labelling idea or at documentation), and require both
inputs to be lists for the animation branch.

**F12 [MEDIUM, confidence: high] — `pt.colors`/`pt.shapes` silently
ignored in the interactive branch; `pt.shapes` in 3D (DOC5)**
`R/gl.pcoa.plot.r:592-651` (no `scale_color_manual`/`scale_shape_manual`),
`:826-853` (`pt.shapes` unused; `pt.colors` is passed to `plot_ly`) — the
`@details` text promises colour/shape control without qualification.
Failure scenario: `interactive = TRUE` with `pt.colors` set returns the
plotly default palette (reproduced: 30 default colours, none from the
supplied vector).
Proposed change: add the manual scales to the interactive branch (built
pre-`ggplotly`, they carry over); document that `pt.shapes` does not apply
in 3D.

**F13 [LOW, confidence: high] — false NOTE in the interactive branch
(DOC5, VRB2)**
`R/gl.pcoa.plot.r:586-590` — "NOTE: Returning the ordination scores, not a
ggplot2 compatable object" — the branch returns a plotly htmlwidget, not
scores (reproduced); "compatable" is also a typo.
Proposed change: correct the message (or drop it once F8 settles the
return contract).

**F14 [LOW, confidence: high] — `as.pop` error message names the wrong
slot (VRB2)**
`R/gl.pcoa.plot.r:335-339` — an unknown `as.pop` metric errors with
"Check names(gl@other$loc.metrics)"; the lookup is in `ind.metrics`
(reproduced).
Proposed change: say `ind.metrics`.

**F15 [LOW, confidence: high] — `scale = TRUE` fixes the ratio at 1; the
computed ratio is dead code (DOC5, STY1)**
`R/gl.pcoa.plot.r:566-575` (and twin blocks `:629-639, :721-731,
:789-799`) — `s1`, `s2`, `r` are computed then discarded;
`coord_fixed(ratio = 1)` is applied regardless of the axes chosen
(reproduced: ratio 1 for axes 1/2 and 1/5 alike). Equal data-unit scaling
is a defensible reading of "scaled to represent the proportion of
variation explained", but the commented-out `coord_fixed(ratio = r)`
lines show an abandoned different intent.
Proposed change: delete the dead computation and commented alternatives
(custodian to confirm ratio-1 is the intended semantics); tighten the
`scale` doc text.

**F16 [LOW, confidence: high] — legend titled "pop", defeating the
`Population` mapping (STY1)**
`R/gl.pcoa.plot.r:660-684` — the legend branch maps
`color = Population` for a capitalised legend title, but
`geom_point(aes(color = pop))` overrides it; the legend renders titled
"pop" (reproduced).
Proposed change: drop the redundant layer-level `aes(color = pop)` in
that branch.

**F17 [LOW, confidence: medium] — house plot bundle absent (PLT1)**
Whole function — no `plot.theme` parameter (theme elements are hardcoded
bold-italic/black at four sites), and colour control is `pt.colors`
rather than the `plot.colors` bundle. The function predates the bundle;
aligning it is an API addition, not a bug fix.
Proposed change: adopt `plot.theme` (default `theme_dartR()`) when the
function is next reworked; keep `pt.colors`/`pt.shapes` as the
documented per-population interface.

**F18 [INFO, confidence: high] — documentation and hygiene bundle (DOC1,
DOC2, DOC7 (proposed rule), STY1)**
Roxygen: tag order deviates from the house order (`@family`,
`@seealso`, `@export` trail the examples); the `verbose` text is the old
wording (DOC2); `@author` has a Custodian but no `Author(s):` line
(DOC7, proposed rule). Body: ~27-line commented-out dist block
(`:467-493`), commented-out `save2tmp` block (`:860-879`), commented-out
as.pop restore (`:891-892`), an unreachable "Plotting entities from the
Distance Matrix" message (`:512-514` — `datatype2` cannot be a dist), and
`hold_x`/`hold_glPca` copies made for all inputs but used only by the
animation branch.
Proposed change: docs-only tidy plus dead-code removal.

### Notes on other functions (scope rule: one line each, no action)

- `utils.plot.save` — `dir = NULL` defaults to `getwd()`, not `tempdir()`,
  the root of F4's cwd write and at odds with sibling functions' "saved to
  tempdir" doc language.
- `utils.check.datatype` — data.frame classifies as "list" (already F8 of
  its own review); it is what routes a data.frame `glPca` into the
  animation branch here.
- `gl.pcoa` — on corrected dist ordinations, stores uncorrected vectors as
  `$scores` but corrected `vectors.cor` as `$loadings`; that asymmetry is
  what F5's classification heuristic trips over (gl.pcoa reviewed
  separately; PR #369 does not touch this).
- `gl.dist.pop` — produces pop-level distance matrices whose ordinations
  no current dartR.base function can plot with pop labels (the code that
  once did is the commented-out block in this function).

## Proposed changes

1. Restore or remove the `pop.labels = "ind"` mode; align gate, branches
   and docs (F1).
2. Wire `hadjust`/`vadjust` into label placement, or remove both
   parameters (F2, F10). **Consequence if removed: signature change —
   callers passing them by name error (API2).**
3. Clamp axis choices to `ncol(glPca$scores)` with informative errors for
   ordinations holding too few axes (F3).
4. Add `gl.check.wd(plot.dir, verbose = 0)` to the preamble so `plot.file`
   saves to the documented location (F4).
5. Classify PCA vs PCoA from a positive signal, not `is.null($loadings)`
   (F5).
6. Gate all warnings/messages per VRB3 and fix the `verbose >= 0` typo
   (F6); add a `plot.display` gate so `verbose = 0` shows no plot (F7).
   **Consequence: a new parameter is added; `verbose = 0` calls stop
   displaying plots (API1).**
7. Replace `return(-1)` dependency guards with `stop(error(...))` (F9).
   **Consequence: code that tested for `-1` now sees an error.**
8. Correct `@return` to describe the returned plot object (F8); optionally
   return `invisible(plott)`. **Consequence if invisible: unassigned calls
   no longer print/render the returned object a second time (API1).**
9. Guard input pairing: `nrow(scores) == nInd(x)`, and both-or-neither
   lists for the animation branch (F11).
10. Honour `pt.colors`/`pt.shapes` in the interactive branch; document the
    3D limitation (F12).
11. Message fixes: interactive NOTE (F13), `as.pop` slot name (F14),
    legend title override (F16).
12. Dead-code and docs tidy: dead `scale` ratio computation and commented
    blocks (F15, F18); roxygen order/DOC2/DOC7 (F18). Docs-only plus
    deletions.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Axis truth: SNP (axes 1/2, 2/5, z=3), SilicoDArT, dist with 9 negative
  eigenvalues, corrected dist — run (PASS; F5 label wording aside)
- Data truth: 2D/3D coordinates vs scores, ind/pop mapping, shuffled-pop
  fixture, `as.pop` — run (PASS)
- Parameter sweep: `xaxis/yaxis/zaxis`, `pop.labels` all four values plus
  invalid, `ellipse`/`plevel`, `scale`, `pt.size/pt.colors/pt.shapes`,
  `label.size` (code-read; see below), `hadjust/vadjust`, `interactive`,
  `as.pop`, `plot.file/plot.dir`, `verbose` 0/2 — run
- Input variants: glPca from SNP, from SilicoDArT, from ind-dist
  (corrected and not), from pop-dist, fd unwrap path (code-read),
  data.frame/matrix/dist wrong-class probes, nfactors 1/2, 2-individual
  object, single-pop object — run
- Contract: return class/visibility, input untouched, cwd/sink balance,
  text silence at `verbose = 0` — run
- `label.size` end value inside `geom_dl`: code-read only — the method
  list is not cleanly introspectable from the built object under
  ggplot2 4.x; the parameter is passed as `cex` in the `smart.grid`
  method list (`:540-541`)
- Animation ("list") branch with a genuine simulation fixture: SKIPPED —
  requires dartR.sim output carrying `$other$sim.vars$generation`; only
  the gate behaviour (F11) was exercised
- FBM-backed genlight as `x`: SKIPPED — metadata-only use of `x` here;
  the FBM concerns live in `gl.pcoa` (PR #369)
- Rendered-image comparison: not attempted by design — object
  introspection per the family mode

## Approval

Approved by Arthur Georges on 2026-09-07 via the formal approval boxes,
acknowledging the consequences recorded against changes 2, 6, 7 and 8.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1) | approved | Arthur Georges, 2026-09-07 | Implement the missing `pop.labels = "ind"` branch; the documented mode starts working |
| 2 (F2, F10) | approved | Arthur Georges, 2026-09-07 | Wire `hadjust`/`vadjust` into label placement rather than remove them; plots with non-zero values change, and the defaults (1.5, 1) are non-zero. Fix the `vadjust` range check that tests `hadjust` |
| 3 (F3) | approved | Arthur Georges, 2026-09-07 | Replace the constant-writing clamp with bounds checking; small ordinations plot or stop informatively instead of "subscript out of bounds" |
| 4 (F4) | approved | Arthur Georges, 2026-09-07 | Route saves through `gl.check.wd`; files that landed in `getwd()` now land in `tempdir()` unless `plot.dir` is given |
| 5 (F5) | approved | Arthur Georges, 2026-09-07 | Corrected-distance ordinations labelled "PCoA Axis" |
| 6 (F6, F7) | approved | Arthur Georges, 2026-09-07 | `verbose = 0` fully silent; new `plot.display` argument, default TRUE, so the plot is not forced |
| 7 (F9) | approved | Arthur Georges, 2026-09-07 | DEP guards become `stop(error(...))` |
| 8 (F8) | approved | Arthur Georges, 2026-09-07 | `@return` documents the visible ggplot the function actually returns; the return stays visible |
| 9 (F11) | approved | Arthur Georges, 2026-09-07 | Informative guard where `nrow(scores) != nInd(x)`, and both-or-neither lists for the animation branch |
| 10 (F12) | approved | Arthur Georges, 2026-09-07 | `pt.colors`/`pt.shapes` honoured in the interactive branch; 3D shape limitation documented |
| 11 (F13, F14, F16) | approved | Arthur Georges, 2026-09-07 | Message fixes |
| 12 (F15, F17) | approved | Arthur Georges, 2026-09-07 | Dead `scale` ratio computation removed; `plot.theme` adopted |
| 12 (F18) | no action | Arthur Georges, 2026-09-07 | INFO: roxygen tag order, DOC2 wording, DOC7 author line and the commented-out blocks left as they are |

## Outcome

Phase C applied on branch `review-gl.pcoa.plot` off `upstream/dev` (ddaed27);
PR #381 into `dev`.
Seventeen findings applied (F1-F17); F18 left as approved.

**Verification** (ggplot/plotly object introspection under `pdf(NULL)`; no
rendered-image comparison):

- **Axis truth and data truth preserved.** Every Phase A equality still holds
  exactly: axis label percentages equal
  `round(100 * eig[k] / sum(eig[eig >= 0]), 1)` on the 2D (axes 1/2 and 2/5),
  SilicoDArT, distance (9 negative eigenvalues) and 3D branches; plotted
  coordinates are identical to the score columns on all four 2D label modes
  and in the 3D traces; `ind` equals `indNames(x)` and `pop` tracks `pop(x)`
  1:1 on the shuffled-pop fixture and under `as.pop = "sex"`.
- **F1**: `pop.labels = "ind"` builds a 4-layer ggplot whose `GeomDl` layer
  data carries `indNames(x)` as labels.
- **F2**: `hadjust`/`vadjust` of 0 versus 3 give different `hjust`/`vjust`
  in the built label-layer data, and those values reach `grid.text()` at
  draw time (traced). The point layer is unchanged. `vadjust = 5` now warns
  and resets to the documented 1; `hadjust` is no longer what the `vadjust`
  check tests.
- **F3**: `gl.pcoa(gl, nfactors = 1)`, a two-individual object and
  `zaxis = 5` on a two-factor ordination all stop with "a 2D plot requires
  at least 2 axes" / "a 3D plot requires at least 3 axes"; `xaxis = yaxis`
  stops with "must differ". No "subscript out of bounds" on any of them.
- **F4**: `plot.file` with the default `plot.dir`, called from a working
  directory that is not `tempdir()`, writes into `tempdir()` and not into
  the working directory. An explicit `plot.dir` is still honoured.
- **F6, F7**: `verbose = 0` produces zero lines on the default call, the
  `"none"` branch, the invalid-parameter paths and the interactive branch,
  and draws nothing (device display list empty). At `verbose = 2` the plot
  renders, and `plot.display = FALSE` suppresses it while the returned
  object is unaffected.
- **F5**: a Cailliez-corrected distance ordination is titled "PCoA Axis 1".
- **F11**: an ordination of `gl.dist.pop()` output stops with a message
  naming the entity mismatch; a data.frame passed as `glPca` stops with the
  animation-pairing message instead of "this S4 class is not subsettable".
- **F12**: the interactive branch returns a plotly object whose marker
  colours are exactly the supplied `pt.colors` and whose symbols are the
  supplied `pt.shapes`.
- **Baseline**: 15 assertions flipped, each mapping to an approved finding
  and tagged `# [approved Fn]` in `tests/testthat/test-gl.pcoa.plot.R`
  (F1 x1, F2/F10 x1, F3 x2, F4 x1, F5 x1, F6 x4, F11 x2, F13 x1, F14 x1,
  F16 x1). No unexplained diff. The updated file passes in full.

**Not verified**: F9 (dependency guards) is a code change only -- removing
`directlabels`, `plotly`, `gganimate` or `tibble` from the library was not
attempted.

**Caller impact**: `gl.assign.pca()` in dartR.captive and dartR.popgen calls
`gl.pcoa.plot(pcoa, x, ellipse = TRUE, plevel = plevel, verbose = 0)` for
the display side effect and discards the result. Under change 6 that call no
longer displays a plot. Nothing errors and no returned value changes; those
two functions need their own follow-up if the plot is wanted. The
dartRstartup tutorial calls use the default verbosity and still display,
with label positions moved by change 2. No other caller in the eight clones.

```json
{
  "function": "gl.pcoa.plot",
  "package": "dartR.base",
  "family_mode": "plotting/report",
  "commit": "ddaed27",
  "skill_version": "2.0.0",
  "reviewer": "Claude Fable 5 via dartr-dev agent",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "axis_truth": "pass",
  "data_truth": "pass",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5/FS5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS7/PLT2", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3/VRB5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 6},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5/FS10", "status": "applied", "change": 8},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 7},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 2},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5/FS5", "status": "applied", "change": 9},
    {"id": "F12", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 10},
    {"id": "F13", "severity": "LOW", "confidence": "high", "rule": "DOC5/VRB2", "status": "applied", "change": 11},
    {"id": "F14", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 11},
    {"id": "F15", "severity": "LOW", "confidence": "high", "rule": "DOC5/STY1", "status": "applied", "change": 12},
    {"id": "F16", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "applied", "change": 11},
    {"id": "F17", "severity": "LOW", "confidence": "medium", "rule": "PLT1", "status": "applied", "change": 12},
    {"id": "F18", "severity": "INFO", "confidence": "high", "rule": "DOC1/DOC2/DOC7/STY1", "status": "no-action", "change": 12}
  ],
  "other_function_notes": [
    "utils.plot.save: dir=NULL defaults to getwd(), not tempdir()",
    "utils.check.datatype: data.frame classifies as 'list' (its own F8)",
    "gl.pcoa: corrected dist path stores uncorrected $scores with corrected $loadings",
    "gl.dist.pop: pop-level ordinations currently unplottable with labels"
  ],
  "coverage_skipped": [
    "animation/list branch with real simulation fixture (needs dartR.sim output)",
    "FBM-backed x (metadata-only use here; FBM lives in gl.pcoa PR #369)",
    "label.size end value in geom_dl (code-read only)"
  ],
  "datasets": ["testset.gl", "testset.gs", "gl.dist.ind(testset.gl[1:30,])", "gl.dist.pop(testset.gl)", "nfactors 1/2/3/5 fixtures", "shuffled-pop fixture"],
  "baseline_test": "tests/testthat/test-gl.pcoa.plot.R",
  "approved_by": "Arthur Georges",
  "approved_date": "2026-09-07",
  "applied_branch": "review-gl.pcoa.plot",
  "axis_truth_after": "pass",
  "data_truth_after": "pass",
  "baseline_flips": 15,
  "status": "pr-open",
  "pr": 381
}
```
