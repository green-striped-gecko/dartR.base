# Review: gl.filter.factorloadings (dartR.base)
- Family mode: modify
- Date: 2026-09-11
- Reviewer: Claude (Claude Fable 5), dartr-function-review v1.0.0
- Package commit: f5e7b72 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.filter.factorloadings.r`)
- Datasets: testset.gl (SNP) — dartR.data 1.2.5; glPca fixture from `gl.pcoa(testset.gl)`
- Baseline: tests/testthat/test-gl.filter.factorloadings.R (new file, snapshot captured pre-review; 16 assertions, all passing)

## Verdicts

**Standards: Needs work** — the plot bundle is correctly decoupled from
display and the function is silent at `verbose = 0`, but the history trail
names only internal calls, the glPca check is a bare `stop()` behind an
ungated `cat`, and `axis` is never validated.

**Spec: Rework** — `retain = TRUE` does the opposite of its documentation:
it returns exactly the same object as `retain = FALSE`, so the parameter has
never worked; and a pca that does not match the genlight is not rejected but
silently recycled into misaligned locus–loading pairs. The filter's central
contract (which loci go, keyed to which loadings) fails on both axes.

## Independent verification (spec axis)

On `gl.pcoa(testset.gl)` (111 polymorphic loci) with `threshold = 0.1`:
31 loci have |loading| >= 0.1 on axis 1. `retain = FALSE` returns the 80
low-loading loci (correct). `retain = TRUE` also returns the identical 80
low-loading loci — not the 31 the documentation promises — while printing
"Retaining 80 of 111 loci with loadings greater than or equal to 0.1".
Locus–loading pairing was checked positionally: `pca$loadings` carries no
rownames on the `gl.pcoa` path, so the function's `cbind(locNames(x), ...)`
positional pairing is the only alignment mechanism and holds on the
canonical example (counts 111 = 111) but is entirely unchecked.

## Findings

**F1 [HIGH, confidence: high] — retain = TRUE is inverted and identical to retain = FALSE (DOC5 (proposed rule); recurring inverted-filter class)**
`R/gl.filter.factorloadings.r:133-136` — the retain branch selects the
high-loading loci into `tmp` and then calls `gl.drop.loc` on them, removing
exactly the loci it was asked to retain. Since the else branch keeps the
complement, both parameter values return the same object.
Failure scenario: any user following the docs ("if TRUE ... holds only the
loci that load high") gets the complement of what they asked for: 80 loci
instead of 31 on the testset fixture, with a progress message asserting the
opposite. Every published use of `retain = TRUE` has filtered the wrong way.
Proposed change: `gl.keep.loc(x, tmp$locus)` in the retain branch (message
already matches the corrected behaviour).

**F2 [HIGH, confidence: high] — a pca that does not match the genlight is silently recycled (DAT2, DAT5)**
`R/gl.filter.factorloadings.r:127-132` — loadings are paired to
`locNames(x)` positionally by `cbind`, with no check that the counts match;
`pca$loadings` has no rownames on the `gl.pcoa` path, so no name-based
recovery is possible. On a length mismatch R warns ("number of rows of
result is not a multiple of vector length") and recycles.
Failure scenario: passing a genlight subset (or any object other than the
one the pca was computed from) produces misaligned locus–loading pairs and
the filter drops essentially random loci, flagged only by a recycling
warning that testthat had to be told to expect. Confirmed with a 150-locus
subset against the full-data pca.
Proposed change: after the internal monomorph removal, stop with a clear
error when `nLoc(x) != nrow(pca$loadings)` (message naming the likely
cause: pca computed on a different object).

**F3 [MEDIUM, confidence: high] — history names internal calls, not this function (FS8)**
`R/gl.filter.factorloadings.r:124,135,139` — the returned object carries
two new history entries, `gl.filter.monomorphs(x = x, verbose = 0)` and
`gl.keep.loc(x = x, loc.list = loclist <- tmp$locus, verbose = 0)`, and
none for `gl.filter.factorloadings` itself.
Failure scenario: `gl.print.history()` shows internal plumbing with an
internal variable name; replaying it does not reproduce the filtering.
Proposed change: restore the input's history and append this function's
`match.call()` as the single entry.

**F4 [MEDIUM, confidence: high] — glPca check: ungated cat + bare stop(), wrong message (FS5, VRB2)**
`R/gl.filter.factorloadings.r:112-118` — `cat(error("To report factor
loadings, require a glPca object\n")); stop()` prints even at `verbose = 0`
outside the condition system and raises an error with an empty message; the
text is copied from the report twin ("To report..."); the check also
overwrites the `datatype` variable from the genlight check.
Failure scenario: `tryCatch(..., error = conditionMessage)` gets an empty
string; scripted callers cannot tell what failed. Confirmed.
Proposed change: `if (!inherits(pca, "glPca")) stop(error("To filter on
factor loadings, a glPca object is required\n"))` — the idiom the reviewed
twin uses.

**F5 [MEDIUM, confidence: high] — axis is never validated (FS5)**
`R/gl.filter.factorloadings.r:127` — `pca$loadings[, axis]` with an axis
beyond the retained axes fails with "subscript out of bounds". The report
twin gained the range check in PR #265; the filter did not.
Failure scenario: `axis = 99` errors uninformatively after the datatype
checks have passed. Confirmed.
Proposed change: the twin's range check (1 <= axis <= ncol(pca$loadings)).

**F6 [LOW, confidence: high] — documented `...` is never forwarded (DOC5 (proposed rule))**
`R/gl.filter.factorloadings.r:186-190` — `@param ...` promises ggsave-style
parameters for saving, but `utils.plot.save` is called without `...`; the
twin forwards it since PR #265.
Failure scenario: `width = 10` is silently ignored.
Proposed change: forward `...` to `utils.plot.save`, as the twin does.

**F7 [LOW, confidence: high] — documentation and style gaps (DOC1, DOC2, STY3)**
`R/gl.filter.factorloadings.r:69` — `@return "The unchanged genlight
object"` is false (the filtered object, returned invisibly); `:3` `@family
matched filters` (house term is "matched filter"); `:29-31` nonstandard
verbose text (DOC2); `:145-146` axis labels render as "... axis PRE-FILTER
1"; `:135,139` `loclist<-tmp$locus` assignment-in-argument instead of
`loc.list = tmp$locus`; `:129-132` the character round-trip through
`cbind`/`as.numeric` where a plain `data.frame(locus, loading)` serves; no
`missing(threshold)` check despite `[required]` (R's default late error is
serviceable but unlabelled).
Failure scenario: misleads readers of the manual; no runtime damage.
Proposed change: correct the doc fields and labels, name the argument,
build the data frame directly, add a threshold-missing check.

## Proposed changes

1. Fix the retain branch: keep the high-loading set (`gl.keep.loc`) when
   `retain = TRUE` (F1).
   **Consequence: numerical output changes for every existing
   `retain = TRUE` caller — the returned object becomes the documented
   high-loading set (31 loci on the fixture) instead of its complement
   (80 loci, identical to retain = FALSE).**
2. Stop with a clear error when the pca does not match the genlight
   (locus-count check after monomorph removal) (F2).
   **Consequence: calls that previously returned silently corrupted
   results now error.**
3. Single history entry naming this function (F3).
4. glPca check via `inherits` + `stop(error(...))` with corrected message
   (F4).
5. Validate `axis` against the retained axes (F5).
6. Forward `...` to `utils.plot.save` (F6).
7. Documentation and style repairs (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (n/a — no Suggests use), PLT,
  STY — run
- Spec: retain semantics vs docs — run (fixture + independent count);
  pca/genlight alignment — run (canonical and mismatched); plot-save path —
  run; silence at verbose 0 — run; input untouched — run
- FBM path (DAT6): SKIPPED — no FBM fixture; the function densifies nothing
  itself (subsetting only)
- GitHub issues: searched "factorloadings" — only pre-campaign PR #90 (doc
  fix); no open issues
- Google Group search: SKIPPED — not queried this session

Out-of-scope note (one line, no action): `gl.report.factorloadings`'s locus
column shows data-frame row numbers, not locus names, because
`gl.pcoa$loadings` carries no rownames on this path — worth a look when
`gl.pcoa` (PR #369) settles.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | consequence (retain = TRUE output changes) approved explicitly |
| 2 | approved | Arthur | consequence (previously corrupted calls now error) approved explicitly |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |

## Outcome

(pending Phase C)

```json
{
  "function": "gl.filter.factorloadings",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "1.0.0",
  "commit": "f5e7b72",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS8", "status": "proposed", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "proposed", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "proposed", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group: not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
