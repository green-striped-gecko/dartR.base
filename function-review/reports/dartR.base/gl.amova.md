# Review: gl.amova (dartR.base)

- Family mode: analysis
- Date: 2026-09-07
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203 (integration-local HEAD; `R/gl.amova.r` is identical to upstream/dev ddaed27 — `git diff upstream/dev -- R/gl.amova.r` empty, so `load_all` exercised the reviewed code)
- Datasets: testset.gl (250x255), testset.gs (218x255), bandicoot.gl (example), derived fixtures (two-pop and single-pop subsets via `gl.keep.pop`, relabelled/shuffled copies, all-NA individual, all-monomorphic subset, random-pop weak-structure subset)
- Baseline: tests/testthat/test-gl.amova.R (new; 12 tests / 30 assertions, all passing pre-review)

## What the function actually computes

One-level AMOVA. Unless `distance` is supplied, it computes an
individual-by-individual Nei's D matrix with `StAMPP::stamppNeisD(x, FALSE)`
(row order = `indNames(x)`, verified), then delegates to
`pegas::amova(distance ~ pop.names, nperm = permutations)` with `pop.names =
factor(pop(x))`, via a private environment to dodge the amova-name clash the
header mentions. There is no hierarchy argument — populations are the only
stratum; no region/pop nesting. Significance is pegas's permutation test:
one-tailed, `p = sum(randomised sigma2 >= observed)/(nperm + 1)` (no +1 in the
numerator, so p = 0 is attainable; confirmed against pegas 1.3 source).
Determinism comes only from the caller's `set.seed()` — there is no seed
parameter. Phi statistics are not returned in `$varcomp` but appear when the
object is printed (pegas `print.amova`).

## Independent verification (spec axis 1)

Four-way agreement on testset.gl (30 pops, 250 ind), seed 42, nperm 99:

| Quantity | gl.amova | direct pegas on stamppNeisD | hand-computed (Excoffier 1992) | StAMPP::stamppAmova |
|---|---|---|---|---|
| SSD among/within/total | 0.111257 / 0.005057 / 0.116315 | identical | identical (1e-9) | identical |
| sigma2 among/within | 4.595238e-04 / 2.298811e-05 | identical | identical | identical |
| phi_ST (derived) | 0.952357 | — | 0.952357 | — |
| P.value | 0 | 0 | — | 0 |

The delegation passes the right distance with the right labels: components are
invariant to individual order and to non-alphabetical pop labels (verified on a
shuffled, `zzz_`/`aaa_`-relabelled copy). The gl.dist.pop F1 label class does
not occur here. gl.amova does not consume `gl.dist.pop`, `utils.dist.ind.snp`
or the sorensen path by default, so none of the open distance-chain defects
(#376, #315, #375) propagate into the default path.

## Verdict

**Standards: Needs work** — the preamble conforms (FS2–FS4, verbose-0 fully
silent), but the pegas guard uses `cat` + `return(-1)` instead of `stop`, the
function does no input validation of its own, and two thirds of the body is
dead code that densifies the whole genotype matrix to read pop names.

**Spec: Needs work** — the default computation verifies exactly against three
independent implementations, but the `distance` argument is applied
positionally with labels ignored (silently wrong results for any re-ordered
matrix), and SilicoDArT input is admitted although the docs describe SNP
genotypes only.

## Findings

**F1 [HIGH, confidence: high] — supplied `distance` used positionally, labels ignored (FS5, DAT5)**
`R/gl.amova.r:69-72` — `dd <- distance` with no validation: dimension is not
checked against `nInd(x)`, and row/column names are never consulted.
Failure scenario: a distance computed elsewhere with individuals in a
different order (e.g. sorted output from another tool) but correct labels is
accepted silently; on testset.gl a row-permuted copy of the very same NeisD
matrix returns sigma2 among = 8.0e-06 instead of 4.6e-04 — no warning, no
error. A wrong-dimension matrix dies inside pegas with bare
"subscript out of bounds".
Proposed change: validate `distance` — require size `nInd(x)`; if dimnames
(or `Labels` for dist) are present, check them against `indNames(x)` and
reorder or stop with an informative message; document that unlabelled input
is taken to be in `indNames(x)` order.

**F2 [MEDIUM, confidence: high] — dependency guard does not stop (DEP1)**
`R/gl.amova.r:54-62` — the pegas guard prints via `cat(error(...))` and
`return(-1)` instead of `stop(error(...))`.
Failure scenario: on a machine without pegas, a scripted pipeline receives
the numeric `-1`, and downstream code (`res$varcomp`) fails later with an
unrelated message.
Proposed change: use the DEP1 idiom `stop(error(...))`. (Not executed
empirically — pegas is installed; static finding.)

**F3 [MEDIUM, confidence: high] — no guard for fewer than two populations (FS5)**
`R/gl.amova.r:69-104` — a single-population object runs to completion and
returns a table with df = 0 and MSD/sigma2 = NaN, silently.
Failure scenario: an object whose pop slot collapsed (e.g. after
`gl.keep.pop` with one name) yields an all-NaN AMOVA the user may pass on.
Proposed change: `stop(error(...))` when `nPop(x) < 2` (matching the checks
in other analysis functions).

**F4 [MEDIUM, confidence: high] — SilicoDArT admitted against the documented contract (DAT7)**
`R/gl.amova.r:51` — `utils.check.datatype(x)` uses the default `accept`, so
testset.gs runs and returns components (sigma2 = 0.01461 / 0.00195), while
`@description`/`@param x` describe SNP genotypes. Nei's D on ploidy-1
presence/absence scores is an undocumented, unvalidated basis.
Failure scenario: a SilicoDArT user gets numbers with no warning that the
distance basis was never designed for P/A data.
Proposed change: either `accept = "SNP"` or document the P/A behaviour as
intended — custodian's call.

**F5 [MEDIUM, confidence: high] — dead-code densification of the full genotype matrix (STY2; DAT6 proposed)**
`R/gl.amova.r:74-95` — the block builds a full StAMPP-format data frame
(`as.matrix(x)`, ploidy scaling, NaN substitution, per-individual pop-number
loop, cbind of the whole genotype matrix into a data.frame) but only column 2
— the pop names — is ever read (line 96). `pop(x)` already provides this.
Failure scenario: for a large object the function materialises and copies the
entire genotype matrix twice for no computational effect; on FBM-scale data
this is the difference between instant and out-of-memory.
Proposed change: replace lines 74-96 with
`pop.names <- factor(as.character(pop(x)))`. No numerical change (the
baseline test suite pins this).

**F6 [LOW, confidence: high] — NA in the distance propagates to an all-NaN result with only a generic warning (DAT5, VRB4 proposed)**
`R/gl.amova.r:69-104` — an all-NA individual makes every stamppNeisD entry
involving it NaN; the entire result table becomes NaN. The only signal is
pegas's "at least one missing value in the distance object" warning, which
does not name the individuals or suggest `gl.filter.allna`/ind-callrate
filtering.
Proposed change: after computing/receiving `dd`, check for non-finite
entries and stop (or warn at verbose >= 1) naming the offending individuals.

**F7 [LOW, confidence: high] — verbose levels 3/5 promise a results summary that never prints; `@return` inaccurate (DOC5 proposed)**
`R/gl.amova.r:19-21,32-34,107-111` — output is identical at verbose 2, 3 and
5 (begin/end plus the datatype line only); no results summary exists. The
`@return` text says "a vector of variance components" — with permutations the
`varcomp` element is a data.frame of sigma2 and P.value, and the object also
carries `varcoef` and `call`; phi statistics appear only via pegas's print
method.
Proposed change: print the component/phi summary at verbose >= 3, and correct
`@return` to describe the pegas `amova` object as returned.

**F8 [LOW, confidence: high] — roxygen order and verbose wording off-canon (DOC1, DOC2)**
`R/gl.amova.r:1-34` — `@return` sits after `@export` (house order puts it
before `@author`); no `@details` tag; the `verbose` text uses the outdated
"[default 2, unless specified using gl.set.verbosity]" instead of the DOC2
default clause.
Proposed change: reorder tags to the ratified order and adopt the DOC2
wording; run `devtools::document()` in the same change (DOC4).

**F9 [LOW, confidence: high] — author block lacks the Author(s)/Custodian structure (DOC7) (proposed rule)**
`R/gl.amova.r:23-24` — "Bernd Gruber (bugs? Post to ...)" has neither an
`Author(s):` nor a `Custodian:` label.
Proposed change: `Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
\url{https://groups.google.com/d/forum/dartr}`.

**F10 [INFO, confidence: high] — permutation characterisation**
p-values are one-tailed, `k/(nperm + 1)` with no numerator correction (p = 0
attainable); reproducibility requires the caller's `set.seed()` (verified:
identical seeded runs match, different seeds differ). `permutations` is
honoured (p granularity 1/(nperm+1) at nperm 9 and 999); `permutations = 0`
runs and returns components without P.values. No seed parameter exists —
worth an `@details` sentence, not a defect.

**F11 [INFO, confidence: high] — FBM path densifies via interim fix**
`R/gl.amova.r:64-65` — the flagged `#!# intermediate fbm fix` converts
FBM-backed objects with `gl.fbm2gen` (full densification, DAT6 proposed).
Acknowledged interim measure; superseded when F5's fix lands, since only
pop names and stamppNeisD need the data.

**F12 [INFO, confidence: medium] — `print.amova` S3 method is session-dependent**
The header's namespace-conflict note is real: with ade4 loaded after pegas,
ade4's `print.amova` overrides pegas's (observed in a bare Rscript session);
under `load_all(dartR.base)` pegas's method wins and prints components plus
phi correctly. Not a gl.amova defect; the returned data are identical either
way.

## Proposed changes

1. Add `distance` validation: size check against `nInd(x)`, label
   check/reorder against `indNames(x)` when dimnames are present, informative
   `stop(error(...))` otherwise; document the positional assumption for
   unlabelled input (F1).
2. Replace the pegas guard with the DEP1 `stop(error(...))` idiom (F2).
3. Add an `nPop(x) >= 2` guard with an informative error (F3).
4. Restrict to `accept = "SNP"` or explicitly document SilicoDArT support —
   custodian decision (F4). **Consequence: SilicoDArT input that currently
   returns a result would error.**
5. Replace the dead StAMPP-format block (lines 74-96) with
   `pop.names <- factor(as.character(pop(x)))` (F5, F11). No numerical
   change; baseline tests pin equality.
6. Check the computed/supplied distance for non-finite entries and report the
   offending individuals (F6).
7. Print a component/phi summary at verbose >= 3; correct `@return`; adopt
   DOC2 verbose wording; reorder roxygen tags; re-document (F7, F8).
8. Add the Author(s)/Custodian structure to `@author` (F9, proposed rule).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec 1 independent verification: direct pegas + hand-computed Excoffier
  one-level components + StAMPP::stamppAmova, all agreeing to 1e-9 — run
- Spec 2 distance basis: stamppNeisD row order vs indNames; gl.dist.ind
  Euclidean supply path (labels in matching order; accepted) — run;
  distance-chain propagation (#376/#315/#375): not consumed by default — noted
- Spec 3 hierarchy: single mode only (pop); no nesting to exercise — run
- Spec 4 permutations: seed determinism, count honoured, p definition from
  pegas 1.3 source — run
- Spec 5 NA handling: all-NA individual (NaN result), testset.gl's all-NA
  loci (datatype-check warning only) — run
- Spec 6 edge cases: 2 pops, single-pop, single-individual pops (shipped in
  testset.gl — runs), all-monomorphic (zero SSD, runs), SilicoDArT,
  verbose-0 silence (capture.output empty), `permutations = 0` — run
- Spec 7 label integrity: non-alphabetical labels + shuffled individuals — run
- Spec 8 parameters: `x`, `distance` (NULL/matrix/dist/permuted/wrong-dim),
  `permutations` (0/1/9/99/999), `verbose` (0/2/3) — run; packaged example
  (bandicoot.gl, permutations = 1) — run, passes
- pegas-absent guard path: SKIPPED — pegas installed; would require
  library removal (static reading only, F2)
- FBM fixture (DAT6): SKIPPED — no FBM fixture; `dartR_fbm` off in the test
  session; interim fix noted as F11
- poppr.amova cross-check: SKIPPED as redundant — three independent
  computations already agree exactly

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges, 2026-09-07 | `distance` dimension-checked and aligned to `indNames(x)` by labels when present; informative fatal when labels are absent or don't match; correct-order inputs unchanged |
| 2 | Approved | Arthur Georges, 2026-09-07 | House `stop(error())` per DEP1 |
| 3 | Approved | Arthur Georges, 2026-09-07 | Informative fatal instead of the silent all-NaN table |
| 4 | Approved | Arthur Georges, 2026-09-07 | SilicoDArT fatal gate per the documented SNP-only contract; silico runs that returned undefined numbers now stop |
| 5 | Approved | Arthur Georges, 2026-09-07 | Dead StAMPP-format densification block removed; no output change |
| 6 | Approved | Arthur Georges, 2026-09-07 | Non-finite distances fatal, naming the individuals involved |
| 7 | Approved | Arthur Georges, 2026-09-07 | Docs batch: one-level-only stated plainly, p-value definition k/(nperm+1), set.seed guidance; verbose >= 3 summary; `@return` corrected |
| 8 | Approved | Arthur Georges, 2026-09-07 | DOC7 Author(s)/Custodian structure |

INFO items (F10-F12): no action per approval; F10's p-value and seed
characterisation folded into the change-7 `@details` text.

## Outcome (Phase C)

Applied 2026-09-07 on branch `review-gl.amova` (base `upstream/dev`
ddaed27).

- Changes 1-8 applied to `R/gl.amova.r`; INFO items F10-F12 not applied
  (no action per approval).
- Change 1 disposition: a supplied `distance` (matrix or dist) must be
  square, sized `nInd(x)`, and labelled with the individual names
  (dimnames for a matrix, Labels for a dist); it is aligned to
  `indNames(x)` by those labels before use. Unlabelled input is a fatal
  error (the report's fallback of documenting a positional assumption was
  not adopted — the approval requires labels). Duplicate labels on either
  side are also fatal, since name alignment would be ambiguous.
- Change 4 disposition: explicit SilicoDArT gate with a bespoke message
  after `accept = c("SNP", "SilicoDArT")`, following the
  `gl.dist.phylo` precedent, rather than `accept = "SNP"` with the
  generic datatype message.
- Change 6 disposition: implemented as a fatal error (the report offered
  stop or a gated warning); it fires on both computed and supplied
  distances and names the individuals in any non-finite row.
- Change 5 note: the `#!# intermediate fbm fix` densification (F11) is
  retained — `StAMPP::stamppNeisD()` still needs the genotypes; only the
  dead StAMPP-format data-frame block was removed.

Verification (all rerun at the applied state):

- (a) Default-path components match the pinned four-way-verified values
  exactly (seed 42, nperm 99): SSD 0.111257/0.005057/0.116315, sigma2
  4.595238e-04/2.298811e-05, phi_ST 0.952357, P 0/NA; `$tab` and
  `$varcomp` identical (tolerance 1e-15) to the pre-change function run
  side by side in the same session.
- (b) A correctly-ordered labelled matrix gives identical results
  (`identical()` on sigma2); the row-permuted labelled matrix now gives
  the same correct results via name alignment (equal at 1e-12;
  previously sigma2 among = 8.0e-06 instead of 4.6e-04); an unlabelled
  matrix and a wrong-dimension matrix stop with informative messages
  ("carries no individual names...", "covers 10 individuals; the
  genlight object has 250").
- (c) Single-population input stops ("AMOVA requires at least two
  populations; this object has 1"); testset.gs stops ("works only with
  SNP data; Nei's genetic distance is not defined for SilicoDArT
  presence-absence data").
- (d) pegas-absent guard exercised empirically by masking
  `requireNamespace` in the function's environment: stops with the DEP1
  message ("Package pegas needed for this function to work...").
- (e) Dead-block removal changes no output (pinned components hold; old
  vs new `$tab`/`$varcomp` identical) and reduces runtime: 5 default-path
  runs (nperm 0) 3.94 s before, 3.32 s after on testset.gl.
- (f) Original 30-assertion baseline rerun: exactly 4 flips, all mapped
  to approved findings — permuted-matrix pin and "subscript out of
  bounds" pin (F1), SilicoDArT pin (F4), single-pop NaN pin (F3).
  Updated characterization file (32 assertions, flips tagged
  `[approved Fn]`) all green.
- (g) `verbose = 3` end-to-end clean (begin, datatype, AMOVA table,
  variance components, Phi_ST, P.value, end); `verbose = 0` fully silent
  (`capture.output()` empty).
- Caller grep across all 8 dartRverse clones (`R/` and `vignettes/`):
  no callers of `gl.amova` outside its own file — all clear.

```json
{
  "function": "gl.amova",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "upstream_dev_commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5,DAT5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "STY2,DAT6", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DAT5,VRB4", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1,DOC2", "status": "applied", "change": 7},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 8},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "note", "change": null},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "DAT6", "status": "note", "change": 5},
    {"id": "F12", "severity": "INFO", "confidence": "medium", "rule": "DEP2", "status": "note", "change": null}
  ],
  "coverage_skipped": [
    "DEP1 pegas-absent path: pegas installed, static only",
    "DAT6 FBM fixture: none available",
    "poppr.amova cross-check: redundant, three implementations already agree"
  ],
  "status": "pr-open",
  "pr": 378
}
```
