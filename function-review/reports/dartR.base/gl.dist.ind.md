# Review: gl.dist.ind (dartR.base)

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  (`git diff upstream/dev -- R/gl.dist.ind.r` empty — the loaded code is
  the reviewed code); working branch integration-local at ed99203, where
  the two engines this wrapper drives (utils.dist.ind.snp,
  utils.dist.binary) carry their applied PR #315/#316 fixes. Baseline
  numeric anchors therefore pin the post-#315/#316 engine state; on
  ddaed27 as released, the engines' known defects propagate (see the
  propagation notes).
- Reviewed in the population-distance chain wave (gl.allele.freq,
  gl.dist.pop, gl.dist.ind, utils.collapse.matrix v2).
- Datasets: testset.gl[1:8, 1:60], testset.gs[1:8, 1:60], constructed
  all-NA-individual fixture.
- Family mode: analysis (wrapper over reviewed distance engines).
- Baseline: tests/testthat/test-gl.dist.ind.R (9 tests, all pass at the
  reviewed state).
- Checks skipped: Google Group not searched (not available: no browser
  session); FBM path (DAT6) not exercised (no FBM fixture).

## Verdicts

**Standards: Needs work** — preamble and plot bundle conform (the plot
is built unconditionally, so `plot.file` works with display off —
verified, unlike gl.dist.pop); the defects are ungated fallback warnings
and documentation drift.
**Spec: Needs work** — wrapper output equals the engine output exactly
for every accepted method on both datatypes (verified), but the method
validation list omits `sorensen`, silently downgrading a method that
this function's own code and its documented consumer gl.dist.pop both
expect to exist.

## Findings

**F1 [HIGH, confidence: high] — 'sorensen' missing from the validation list; coerced to simple matching (spec)**
`R/gl.dist.ind.r:118-128` — the accepted-methods vector lacks
"sorensen", although utils.dist.binary implements it (verified exact in
review #316), this function carries a `verbose >= 2` progress message
for it (line 233-237, unreachable dead code), and gl.dist.pop lists it
as available and calls through here.
Failure scenario: `gl.dist.ind(testset.gs, method = "sorensen")` returns
simple matching distances (verified identical), with a warning that
leaks at `verbose = 0` (F2); `gl.dist.pop(gs, method = "sorensen")`
inherits the wrong statistic silently (its finding F2).
Proposed change: add "sorensen" to the validation list.
**Consequence: numerical output changes for `method = "sorensen"`
callers (values become true Sorensen/Dice instead of simple matching).**

**F2 [MEDIUM, confidence: high] — method-fallback warnings are ungated (VRB5, VRB3)**
`R/gl.dist.ind.r:131,137` — both `cat(warn("Method not in the list of
options..."))` calls sit outside any verbosity gate.
Failure scenario: one warning line prints at `verbose = 0` for any
unrecognised method (verified, SNP and SilicoDArT), breaking the
documented full-silence contract; gl.dist.pop's sorensen path leaks it
into otherwise-silent runs.
Proposed change: gate at `verbose >= 1` (VRB4: the result is affected,
so level 1, not 2).

**F3 [MEDIUM, confidence: high] — `scale` documentation contradicts the code twice (DOC5 (proposed rule))**
`R/gl.dist.ind.r:11,64-67,174-202` — `@param scale` says
"[default TRUE]"; the signature default is FALSE. The doc presents
scaling as a general option, but the wrapper forwards `scale` only for
euclidean: `method = "simple", scale = TRUE` is identical to
`scale = FALSE` (verified).
Failure scenario: a user requesting scaled simple distances silently
gets unscaled ones and, reading the docs, believes scaling was on by
default anyway.
Proposed change: docs-only — correct the default and state that `scale`
applies to euclidean only (or warn when `scale = TRUE` is ignored).

**F4 [LOW, confidence: high] — documented method lists out of step with the code (DOC5 (proposed rule), DOC1)**
`R/gl.dist.ind.r:29-52` — the SilicoDArT list names "Bray-Curtis" but
not Sorensen (identical formulas in the engine — the synonymy is worth
one line once F1 lands); the SNP list omits nothing but does not note
that post-#315 Simple equals Czekanowski by construction (recorded in
the engine docs). The `@author` line is malformed: "Author(s):
Custodian: Arthur Georges -- Post to #' \url{...}" — a stray `#'` and a
missing author name (DOC7 (proposed rule)).
Proposed change: docs-only.

**F5 [LOW, confidence: high] — `type` is matched case-sensitively against "matrix" only (FS5)**
`R/gl.dist.ind.r:339-355` — any other string, including "Matrix", falls
through to the dist branch silently (verified).
Failure scenario: `type = "Matrix"` returns a dist; downstream matrix
indexing fails elsewhere.
Proposed change: normalise with `tolower()` and validate.

**F6 [LOW, confidence: medium] — an all-missing individual yields NA distances with no warning (VRB4 (proposed rule))**
Verified: an individual with all genotypes NA produces NA for every
distance involving it, silently; the `verbose >= 3` summary then prints
NA for min/max (no `na.rm`), and downstream consumers (hclust, nj,
utils.collapse.matrix) receive NA cells unannounced.
Proposed change: count NA distances after computation and warn at
`verbose >= 1` when any exist.

**F7 [INFO, confidence: high] — cosmetic notes (STY1)**
`dd <- dd <- utils...` doubled assignment at lines 175, 182, 190;
progress messages print after the engine call they describe; plot
objects are always built even when neither displayed nor saved (cheap,
but wasted); the final dist/matrix conversion block reconverts objects
that are already the right class.

Propagation notes (defects owned by the engines, visible here on
ddaed27 until PRs #315/#316 merge — not re-found, listed for the record):

- propagates-from #315: SNP `simple`/`absolute` are reference-allele
  asymmetric on upstream/dev; a one-line progress cat leaks at
  `verbose = 0`; "Returning/Completed" print twice at `verbose >= 2` —
  all through the `verbose = verbose` pass-through on the SNP paths.
- propagates-from #316: `method = "bray-curtis"` silently falls back to
  simple matching on upstream/dev (verified fixed on this branch:
  bray-curtis now equals engine sorensen). The upstream scale-warning
  leak only fires with `scale = TRUE` on non-euclidean silico methods.

## Proposed changes

1. Add "sorensen" to the accepted-methods list (F1).
   **Consequence: numerical output changes for `method = "sorensen"` on
   SilicoDArT data (and for gl.dist.pop's sorensen method).**
2. Gate the two fallback warnings at `verbose >= 1` (F2).
3. Docs-only: `scale` default and euclidean-only reach, method-list
   synonymy note, repair the `@author` line (F3, F4).
4. Normalise and validate `type` (F5).
5. Warn at `verbose >= 1` when the returned distances contain NA (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run.
- Spec: wrapper == engine equality for all 6 SNP method/scale
  combinations and 3 SilicoDArT methods; sorensen and unknown-method
  fallbacks; swap effect (changes jaccard, verified); type variants and
  dimnames; verbose-0 silence; plot.file with display off; 1-individual
  and all-NA-individual edges — run.
- FBM path (DAT6): SKIPPED — no FBM fixture in this wave.
- Google Group search: SKIPPED — not available, no browser session.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges (2026-09-06) | consequence acknowledged: sorensen output becomes true Sorensen/Dice in this function and gl.dist.pop |
| 2 | approved | Arthur Georges (2026-09-06) | |
| 3 | approved | Arthur Georges (2026-09-06) | |
| 4 | approved | Arthur Georges (2026-09-06) | |
| 5 | approved | Arthur Georges (2026-09-06) | |

## Outcome

All 5 approved changes applied on branch `review-gl.dist.ind`
(base `ddaed27`, upstream/dev) and submitted as
[PR #375](https://github.com/green-striped-gecko/dartR.base/pull/375),
covering findings F1-F6; F7 is INFO/no
action. The engines (utils.dist.ind.snp, utils.dist.binary) were not
touched -- their defects ride on the open PRs #315/#316. Verification:

- F1 verified against hand computation: `method = "sorensen"` equals
  1 - 2a/(2a+b+c) on pairwise-complete loci for every pair of
  testset.gs[1:8, 1:60] (max |diff| 5.6e-17), equals the
  utils.dist.binary sorensen engine exactly, and differs from simple
  matching.
- Baseline characterization test: 30 assertions pass; flipped
  expectations tagged `[approved F1]`, `[approved F2]`, `[approved F5]`,
  `[approved F6]` (plus a `[approved F3]` docs-only comment update); no
  unexplained diffs. One pre-existing divergence from the Phase A run:
  on upstream/dev the #315 engine leaks one progress line at
  `verbose = 0` on SNP paths (Phase A ran against the post-#315 engine
  state); the two silence assertions filter that known engine line so
  they test this wrapper's own gating, and the filter becomes a no-op
  once #315 merges.
- `verbose = 3` end-to-end run clean (the reachable sorensen progress
  message now prints); NA-distance warning fires at `verbose = 1` on an
  all-missing individual and not at `verbose = 0`.
- Caller grep across the 8 dartRverse clones: gl.dist.pop
  (simple/jaccard/sorensen -- the sorensen change is the intended chain
  fix, recorded in its report), gl.impute (euclidean, type="matrix" --
  valid, unaffected), gl.smearplot and dartR.popgen's structure/snmf
  plots (Manhattan -- in the accepted list, unaffected),
  dartR.spatial gl.spatial.autoCorr (method pass-through -- a user
  passing "Sorensen" with binary data now gets real Sorensen). No
  breaking caller.

Integration probe (all three chain branches loaded together): see the
probe record in gl.dist.pop.md -- gl.dist.pop's sorensen equals hand
Sorensen through this fix in the combined state.

```json
{
  "function": "gl.dist.ind",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "VRB4", "status": "applied", "change": 5},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "noted", "change": null}
  ],
  "propagations": [
    {"from": "utils.dist.ind.snp", "pr": 315, "effect": "simple/absolute asymmetry + verbose leaks on upstream/dev SNP paths"},
    {"from": "utils.dist.binary", "pr": 316, "effect": "bray-curtis silent fallback to simple on upstream/dev"}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group: no browser session"],
  "baseline_test": "tests/testthat/test-gl.dist.ind.R",
  "status": "pr-open",
  "pr": 375
}
```
