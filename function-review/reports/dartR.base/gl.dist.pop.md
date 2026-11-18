# Review: gl.dist.pop (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  (`git diff upstream/dev -- R/gl.dist.pop.r` empty — the loaded code is
  the reviewed code); working branch integration-local at ed99203.
- Reviewed in the population-distance chain wave (gl.allele.freq,
  gl.dist.pop, gl.dist.ind, utils.collapse.matrix v2). The SilicoDArT
  simple/jaccard/sorensen paths route through gl.dist.ind +
  utils.collapse.matrix, whose engines were loaded at their post-PR
  #315/#316/#325 state on this branch; the SNP methods do not touch
  those engines, so their verification is valid for ddaed27 as such.
- Datasets: testset.gl, testset.gs (4-pop subsets for anchors), a
  constructed clean 3-pop/40-locus fixture with zero missing data, a
  2-pop fully-fixed-opposite fixture; adegenet::dist.genpop as the
  independent reference implementation.
- Family mode: analysis.
- Baseline: tests/testthat/test-gl.dist.pop.R (12 tests, 40 assertions,
  all pass at the reviewed state).
- Checks skipped: Google Group not searched (not available: no browser
  session); FBM path (DAT6) not exercised (no FBM fixture); StAMPP
  comparison abandoned (stamppNeisD rejected the converted testset —
  dist.genpop used instead, which is exact rather than approximate).

## Verdicts

**Standards: Needs work** — the DEP1 guard returns -1 instead of
stopping, an unfatal "Fatal Error" prints ungated at `verbose = 0`, and
the plot block references an object that may not exist.
**Spec: Rework on one axis, verified on the other** — every distance
formula is exactly right (see the verification matrix), but the matrix
is silently mislabelled whenever the population factor's levels are not
in alphabetical order, which makes the numbers wrong under the labels
users read. Formula quality deserves stating: this is the rare review
where all five SNP methods reproduce independent references to machine
precision.

### Method verification matrix

| Method | Hand recomputation (code's own NA policy) | Independent reference | Result |
|---|---|---|---|
| euclidean | exact (max diff 0) | — | verified |
| euclidean, scale=TRUE | exact | see F7: max attainable is 0.5, not the documented 1 | verified, doc wrong |
| nei | exact | == adegenet dist.genpop(1) Nei 1972 to 2.3e-16 on a clean fixture | verified |
| reynolds | exact | == -log(1 - dist.genpop(3)) to 2.3e-15 — a linearised (divergence-time) transform of standard Reynolds, undocumented | verified, doc gap |
| chord | exact | == (2*sqrt(2)/pi) x dist.genpop(2) Edwards, ratio 0.9003 exactly | verified, doc gap |
| fixed-diff | == gl.fixed.diff $pcfd/100 exactly (a dist, label-safe) | gl.fixed.diff already reviewed | verified |
| SilicoDArT euclidean (+scale) | exact from presence frequencies | — | verified |
| SilicoDArT simple/jaccard | == collapse(gl.dist.ind) route exactly | engines verified in #315/#316 | verified |
| SilicoDArT sorensen | returns SIMPLE matching (see F2) | — | defect |

Empirical NA policy: frequencies are per-cell `na.rm = TRUE` (from
gl.allele.freq); each pair then drops loci with a NaN frequency in
either population (pairwise complete cells); unscaled euclidean does not
rescale for the dropped loci. No NA pairs arise on testset.gl (0 of 435).
Cross-note (gl.tree.nj, PR #370): the fix's `na.rm = TRUE` alignment
claim is correct at the frequency level, but `stats::dist` (which
gl.tree.nj uses) rescales for the 545 all-NA cells while gl.dist.pop
sums complete cells unrescaled — the two euclidean matrices still differ
by up to 14.7% (median 5.1%) on testset.gl, so the tree is close to,
not identical to, `gl.dist.pop(method = "euclidean")`.

## Findings

**F1 [BLOCKER, confidence: high] — distances land on the wrong population labels when pop levels are not alphabetical (spec; DAT2 lens)**
`R/gl.dist.pop.r:188-192,366` — `reshape2::dcast` orders rows by the
alphabetical re-factored `popn` from gl.allele.freq, but `dimnames(dd)`
is assigned positionally from `popNames(x)` (the object's own level
order). Affects euclidean, nei, reynolds, and chord on both datatypes;
fixed-diff and the SilicoDArT collapse methods are name-safe (verified).
Failure scenario: reverse the level order of `pop(testset.gl)` (contents
unchanged) — D[pop1, pop2] silently changes from 1.586 to 4.137; hand
check against the object's own labels shows errors up to 4.29. Any
object whose population factor was built with custom level ordering
returns wrong numbers under every label, with no warning. The default
datasets happen to have alphabetical levels, which is why this has not
surfaced.
Proposed change: index the frequency matrix rows by `popNames(x)` after
the dcast (or take dimnames from the dcast's own row names), so labels
and rows always travel together.
**Consequence: numerical output changes (becomes correct) for any object
with non-alphabetical population levels.**

**F2 [HIGH, confidence: high] — method 'sorensen' returns simple matching distances (spec; chain: root cause in gl.dist.ind)**
`R/gl.dist.pop.r:391-398` — the sorensen block calls
`gl.dist.ind(x, method = "sorensen", verbose = 0)`; gl.dist.ind's
validation list omits "sorensen" and coerces it to "simple" (its
finding F1). Failure scenario: `gl.dist.pop(testset.gs, method =
"sorensen")` returns values identical to `method = "simple"` (verified)
and leaks a one-line warning at `verbose = 0`. sorensen is in this
function's own `available_methods`.
Proposed change: none here — fixed by gl.dist.ind change 1; this report
records the user-visible consequence at the gl.dist.pop level.
**Consequence: numerical output changes for `method = "sorensen"` once
the chain fix lands (values become true Sorensen/Dice).**

**F3 [HIGH, confidence: high] — `plot.file` with plotting off crashes after the computation (PLT3)**
`R/gl.dist.pop.r:459-470` — `p3` is built only inside
`if (plot.display)`, but the save block references it whenever
`plot.file` is set. Failure scenario:
`gl.dist.pop(x, plot.file = "f", plot.display = FALSE)` — and equally
any `verbose = 0` call with `plot.file`, since `verbose == 0` forces
`plot.display <- FALSE` — dies with "object 'p3' not found" (verified)
after all distances are computed; the user loses the result.
Proposed change: build the plot objects whenever
`plot.display || !is.null(plot.file)` (the gl.dist.ind pattern, which
verifies clean).

**F4 [MEDIUM, confidence: high] — unknown method: ungated non-fatal "Fatal Error", silent euclidean (VRB5, VRB2, FS5)**
`R/gl.dist.pop.r:157-162` — `cat(error("Fatal Error: ... set to
Euclidean"))` prints at `verbose = 0` (16 lines leaked, verified) and
the function continues with euclidean (verified identical output).
Failure scenario: `method = "neii"` returns euclidean distances; in a
pipeline the leaked text is the only clue, and it claims to be fatal.
Same class as gl.tree.nj's "ugpma" precedent.
Proposed change: either stop with `stop(error(...))`, or warn via
`warn()` gated at `verbose >= 2` — pick one; the current hybrid is the
worst of both.

**F5 [MEDIUM, confidence: high] — `type = "matrix"` output shape depends on the method (spec; DOC5 (proposed rule))**
`R/gl.dist.pop.r:481-487` — for the SNP frequency methods only the lower
triangle is ever filled, so `type = "matrix"` returns a matrix whose
entire upper triangle and diagonal are NA (verified: 435/435 upper cells
NA); fixed-diff and the SilicoDArT collapse methods return full
symmetric matrices with zero diagonals (verified). The function's own
commented-out validation code manually symmetrises before use.
Failure scenario: any code consuming `type = "matrix"` (heatmaps,
`D[i, j]` with i < j) reads NA for half the pairs, method-dependently.
Proposed change: symmetrise and zero the diagonal before the matrix
return path.

**F6 [MEDIUM, confidence: high] — DEP1 guard returns -1 instead of stopping (DEP1)**
`R/gl.dist.pop.r:68-77` — the reshape2 guard does `cat(error(...));
return(-1)`. Failure scenario: without reshape2 installed, callers
receive `-1` where a dist is expected and fail later, far from the
cause.
Proposed change: the DEP1 idiom — `stop(error(...))`.

**F7 [MEDIUM, confidence: high] — `scale = TRUE` documented range [0,1] is unattainable for SNP data (DOC5 (proposed rule))**
`R/gl.dist.pop.r:16-17,203-208` — the SNP scaling is
`0.5*sqrt(mean(sq))`; two populations fixed for opposite alleles at
every locus score 0.5, the true maximum (verified on a constructed
fixture). SilicoDArT scaling (`sqrt(mean(sq))`) does reach 1.
Failure scenario: a user reading 0.31 against a documented [0,1] scale
underestimates divergence by half.
Proposed change: docs-only — state the SNP range as [0, 0.5] (or note
the 0.5 factor's intent, Rogers-type scaling), and while in the header
fix the garbled `@details` sentence (F8).

**F8 [LOW, confidence: high] — `@details` is garbled and incomplete (DOC1, DOC5 (proposed rule))**
`R/gl.dist.pop.r:30-34` — the SilicoDArT sentence breaks off mid-list
("can be one of 'Refer to the documentation..."), never naming
euclidean/simple/jaccard/sorensen; the reynolds -log transform and the
chord scaling constant (verification matrix above) are undocumented.
Proposed change: docs-only — complete the method lists and state the
two formula variants.

**F9 [LOW, confidence: high] — a single population fails with a bare subscript error (FS5)**
`R/gl.dist.pop.r:196` — `1:(nP-1)` with nP = 1 iterates and indexes out
of bounds ("subscript out of bounds", verified). Two populations work
(verified). Proposed change: fail fast with a clear "at least two
populations" message.

**F10 [LOW, confidence: high] — `as.pop` error message points at the wrong slot (VRB2)**
`R/gl.dist.pop.r:138` — the message says "Check
names(gl@other$loc.metrics)" for an individual metric; it should say
`ind.metrics` (verified message text). Proposed change: docs/message
fix.

**F11 [INFO, confidence: high] — cosmetic and efficiency notes (STY1, STY2)**
"Czfordi-Edwards" for Cavalli-Sforza-Edwards in the chord stop message;
the fixed-diff progress message prints after the computation it
announces; `gl.allele.freq` is computed for every method including
fixed-diff and the collapse paths, which never use it; frequency input
carries gl.allele.freq's 2 dp percentage rounding (max 1.1e-05 effect,
benign). The direct SilicoDArT euclidean (frequency-based) and a
collapsed individual-level euclidean are different statistics (cor 0.95,
not equal) — not a defect; worth one doc line if members expect them to
agree.

Propagation notes: #315 (utils.dist.ind.snp) does not reach gl.dist.pop
— its SNP methods are computed from frequencies, not from the individual
engine. #316 (utils.dist.binary) reaches the SilicoDArT collapse paths
only via gl.dist.ind with `scale = FALSE`, so the known upstream scale
warning leak is not triggered; jaccard/simple values were verified
correct through the chain. The one chain defect that does surface here
is gl.dist.ind's own sorensen coercion (F2).

## Proposed changes

1. Re-anchor frequency-matrix rows to `popNames(x)` (F1).
   **Consequence: numerical output changes (corrects silent
   mislabelling) for objects with non-alphabetical population levels.**
2. Build plot objects whenever a plot is displayed or saved (F3).
3. Replace the unknown-method hybrid with a proper stop or a gated warn
   (F4). **Consequence: `method` typos stop (or warn) instead of
   silently running euclidean.**
4. Symmetrise `type = "matrix"` output for the SNP frequency methods
   (F5). **Consequence: matrix output gains upper-triangle values and a
   zero diagonal where callers previously saw NA.**
5. Make the reshape2 guard stop per DEP1 (F6).
6. Docs-only: scaled range, `@details` method lists and formula
   variants, `as.pop` message, chord typo (F7, F8, F10, F11 text parts).
7. Guard nPop < 2 with a clear error (F9).

(Change for F2 lives in gl.dist.ind's report, change 1.)

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run.
- Spec: all 5 SNP methods + 4 SilicoDArT methods hand-verified with the
  code's own NA policy; nei/reynolds/chord/scaled-euclidean checked
  against adegenet::dist.genpop on a clean fixture; fixed-diff against
  gl.fixed.diff; label-safety probe with reversed factor levels;
  method/datatype refusal pairs; 1-pop and 2-pop edges; verbose-0
  silence; plot bundle decoupling; input-object immutability — run.
- FBM path (DAT6): SKIPPED — no FBM fixture in this wave.
- StAMPP cross-check: SKIPPED — stamppNeisD rejected the converted
  object; dist.genpop served as the (exact) independent reference.
- Google Group search: SKIPPED — not available, no browser session.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges (2026-09-06) | BLOCKER consequence acknowledged: output becomes correct for non-alphabetical pop levels; alphabetical datasets unchanged |
| 2 | approved | Arthur Georges (2026-09-06) | |
| 3 | approved | Arthur Georges (2026-09-06) | consequence acknowledged: method typos now stop (stop chosen over gated warn) |
| 4 | approved | Arthur Georges (2026-09-06) | consequence acknowledged: matrix return full-symmetric, zero diagonal; caller grep run first |
| 5 | approved | Arthur Georges (2026-09-06) | |
| 6 | approved | Arthur Georges (2026-09-06) | |
| 7 | approved | Arthur Georges (2026-09-06) | |

## Outcome

All 7 approved changes applied on branch `review-gl.dist.pop`
(base `ddaed27`, upstream/dev) and submitted as
[PR #376](https://github.com/green-striped-gecko/dartR.base/pull/376),
covering findings F1, F3-F10 and the
F11 text parts; F2 is fixed in gl.dist.ind (its change 1, PR #375) and
F11 is otherwise INFO/no action. For F4 the stop option was chosen (a
method typo is an error, and the fallback claimed to be fatal without
being so). Verification:

- F1 verified against hand computation (frequencies recomputed from
  `as.matrix` under the object's own labels, the code's own NA policy):
  max |code - hand| = 0 on the alphabetical 4-pop testset.gl subset, on
  the REVERSED-level fixture, and on a SHUFFLED-level fixture (seed 7:
  EmmacBurdMist, EmmacClarJack, EmmacBurnBara, EmmacBrisWive); nei
  checked shuffled-vs-alphabetical label-aligned equal to 1e-9.
- Alphabetical-level anchors for euclidean, nei, reynolds, chord,
  fixed-diff and scaled euclidean are unchanged to 6 dp (baseline).
- F3: `plot.file` with `plot.display = FALSE` at `verbose = 0` returns
  the dist and writes the RDS (file existence checked).
- Baseline characterization test: 41 assertions pass; flipped
  expectations tagged `[approved F1]`, `[approved F3]`, `[approved F4]`,
  `[approved F5]`, `[approved F9]`; the sorensen assertion is
  chain-state aware (see below); no unexplained diffs.
- `verbose = 0` empirically silent for valid methods; `verbose = 3`
  end-to-end run clean.
- Caller grep across the 8 dartRverse clones: no code caller of
  gl.dist.pop anywhere (documentation references only:
  gl.plot.heatmap examples, gl.tree.nj and utils.collapse.matrix doc
  text). The `type = "matrix"` symmetrisation and the label fix
  therefore break no sibling caller.
- Sorensen on this branch still routes through the unfixed gl.dist.ind
  (this branch bases on upstream/dev): it returns simple matching with
  a leaked warning until PR #375 merges, at which point real Sorensen
  arrives here with no code change on this branch. The baseline test
  pins both states.

Integration probe (branches review-gl.allele.freq + review-gl.dist.ind
+ review-gl.dist.pop loaded together in a scratch state):

- `gl.dist.pop(testset.gs[pops 1-4], method = "sorensen")` equals the
  hand-computed Sorensen collapse (mean over all cross-population
  individual pairs of 1 - 2a/(2a+b+c)) exactly -- real Sorensen now.
- euclidean/nei/reynolds/chord on alphabetical data identical to the
  single-branch (review-gl.dist.pop only) state, and labels correct on
  the reversed-level fixture (max |code - hand| = 0).
- `gl.allele.freq(by = 'popxloc')` output identical between the
  combined state and each single-branch state, feeding gl.dist.pop
  identically (dist matrices byte-identical).

```json
{
  "function": "gl.dist.pop",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "spec", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "fixed-in-chain", "change": null, "note": "fixed by gl.dist.ind change 1, PR #375"},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "applied", "change": 2},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 6},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 7},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 6},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "StAMPP: conversion rejected, dist.genpop used", "Google Group: no browser session"],
  "baseline_test": "tests/testthat/test-gl.dist.pop.R",
  "status": "pr-open",
  "pr": 376
}
```
