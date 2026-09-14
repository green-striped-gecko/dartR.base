# Review: gl.add.indmetrics (dartR.base)
- Family mode: modify (adds/updates individual metadata; subsets & reorders x)
- Date: 2026-09-14
- Reviewer: Claude (Claude Opus 4.8), dartr-function-review v1.0.0
- Package commit: f02bd34 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.add.indmetrics.r`)
- Datasets: testset.gl, dartR.data extdata (testset_SNPs_2Row.csv + testset_metadata.csv) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.add.indmetrics.R (new file, snapshot captured pre-review; 12 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone (verbosity, flag start/end,
datatype check, history append, explicit return) is present, but the
roxygen block has no `@family` tag, the datatype documentation is wrong, and
a couple of error idioms are non-standard.

**Spec: Rework** — the function crashes on exactly the scenario its own
documentation and warning message describe as supported: a metadata file
whose ids are a superset of (or only partly overlap) the genlight. The core
promise — "add metadata, keeping the individuals that match" — fails
whenever the metadata carries any individual not present in `x`.

## Blast radius

`gl.add.indmetrics` has **no internal callers** anywhere in the dartRverse —
it is a user-facing metadata function (also part of the documented
`gl.read.dart` → `gl.add.indmetrics` workflow).

## Independent verification (spec axis)

Four ordering scenarios were exercised (tests 1-4). The documented example
(250 individuals, metadata of exactly 250 matching ids) works; an
equal-size **reordered** metadata file is correctly aligned to `x` (each
individual's population matches its own metadata row); and a strict
**subset** metadata file (fewer rows, all present in `x`) subsets `x`
correctly. The failure is the **superset / partial-overlap** case (test 2),
traced below.

## Findings

**F1 [HIGH, confidence: high] — crashes when the metadata file contains any individual not in x (DAT2; recurring metadata-desync class)**
`R/gl.add.indmetrics.r:130,138` — after matching, the function subsets `x`
to the matched individuals (`x <- x[ord2, ]`, line 130) but leaves `ind.cov`
at its full row count. Line 138 then does `ind.cov$pop_old <- x@pop`,
assigning a vector of length `nInd(x)` (the matched subset) into a column of
the full-length `ind.cov`.
Failure scenario: a master metadata file covering more individuals than the
current genlight — the routine "apply my sample sheet to this subset"
workflow — errors with "replacement has 6 rows, data has 10" (test 2). This
is the very case the function's own warning at line 88 calls out as
supported ("Maybe this is fine if a subset matches"). Only the exact-match
and strict-subset cases survive; any extra metadata row crashes.
Proposed change: align `ind.cov` to the matched, reordered individuals once,
immediately after `x <- x[ord2, ]` — `ind.cov <- ind.cov[ord, , drop =
FALSE]; ord <- seq_len(nrow(ind.cov))` — so `x@pop` and every subsequent
`ind.cov[ord, ]` share the same length and order. (Verified: `x[ord2, ]` is
in `ind.cov[ord, ]` order, so this subset is the correct alignment.)
**Consequence: the superset/partial-overlap case now succeeds (returns the
matched individuals with their metadata) instead of erroring. Behaviour for
the exact-match and strict-subset cases is unchanged.**

**F2 [MEDIUM, confidence: high] — missing @family and wrong datatype documentation (DOC1, DOC2, DOC5 (proposed rule))**
`R/gl.add.indmetrics.r:1-13` — there is no `@family` tag (every dartR
function carries one; the sibling io/metadata functions use
`@family io` or `@family data manipulation`); `@param x` describes "the
genind object containing the SilocoDArT data" (SilicoDArT is stored in a
genlight with ploidy 1, not a genind; "SilocoDArT" is also a typo); the
`verbose` text is the pre-DOC2 form.
Failure scenario: the manual mis-describes the accepted input; the function
does not appear in any `@family` grouping.
Proposed change: add `@family` (matching the sibling metadata functions),
correct the `@param x` datatype wording and typo, adopt the DOC2 verbose
text.

**F3 [LOW, confidence: high] — non-standard error idiom for the not-unique-ids check (FS5, VRB2)**
`R/gl.add.indmetrics.r:77-82` — the duplicate-id check uses
`cat(error("...")); stop()`, printing the message ungated (even at
`verbose = 0`) and then raising an error with an empty message. The
adjacent no-id check (line 71-72) correctly uses `stop(error(...))`.
Proposed change: `stop(error("Individual names are not unique. ..."))`, for
one gated, informative fatal error.

**F4 [LOW, confidence: high] — confused verbose count expression (STY1)**
`R/gl.add.indmetrics.r:120-121` — `length(ord == nInd(x))` reports the
matched-id count, but the expression evaluates `ord == nInd(x)` (a logical
vector) and takes its length; it returns `length(ord)` by accident, which
happens to be the intended number.
Proposed change: `length(ord)` (the actual match count).

**F5 [INFO, confidence: high] — read.csv(stringsAsFactors = TRUE) (STY3)**
`R/gl.add.indmetrics.r:63-66` — reads all character columns as factors.
Numeric lat/lon stay numeric (verified: `latlon` is stored numeric), so
there is no correctness impact here, but it factorises character columns the
user may want as strings. Modern default is FALSE. Recorded; no change
proposed to avoid altering the ind.metrics column types callers may rely on.

## Proposed changes

1. Align `ind.cov` to the matched/reordered individuals after the subset, so
   the superset/partial-overlap metadata case works instead of crashing
   (F1).
   **Consequence: a metadata file with individuals beyond those in `x` now
   succeeds (matched individuals returned with their metadata); exact-match
   and strict-subset behaviour unchanged.**
2. Documentation: add `@family`, correct the `@param x` datatype wording and
   "SilocoDArT" typo, adopt the DOC2 verbose text (F2).
3. Use `stop(error(...))` for the duplicate-id check, and `length(ord)` for
   the match count (F3, F4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (n/a), PLT (n/a), STY — run
- Spec: four ordering scenarios (exact match, reordered, strict subset,
  superset/partial overlap) traced and exercised; latlon type; no-id error;
  verbose-0 silence; history append — run
- FBM path (DAT6): SKIPPED — the example wraps input in gl.gen2fbm but the
  function operates on metadata / individual subsetting, not the genotype
  matrix
- Caller survey (API3): no internal callers; user-facing metadata function
- GitHub issues / Google Group: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | consequence (superset case succeeds instead of crashing) approved explicitly |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.add.indmetrics",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "f02bd34",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "proposed", "change": 3},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "STY3", "status": "proposed", "change": null}
  ],
  "coverage_skipped": ["DAT6: not exercised on FBM", "GitHub issues not queried", "Google Group not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
