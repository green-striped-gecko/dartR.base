# Review: gl.report.maf (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: ddaed27 (dev, synced with upstream/dev —
  note: dev has advanced, now including merges to PR #293); Branch:
  review-gl.report.maf; Datasets: testset.gl (4-population subset and
  full 30-population set)
- Family mode: report (returns invisible per-population MAF dataframe;
  input untouched verified, as.pop restore verified)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the verbosity contract is broken wholesale: a
  hardcoded `verbose = 3` INSIDE the per-population function overrides
  the caller's setting, so per-population statistics print at every
  verbosity (verified: 77 lines at verbose 0 on a 4-pop subset, 343 on
  the full testset); the overall statistics block and the quantile
  table (`print(df)`) are entirely ungated; the maf.limit/ind.limit
  coercion warnings and the singleton-drop notice are ungated; a
  second handwritten FLAG SCRIPT START block duplicates the "Starting"
  banner (verified: 2 at verbose >= 1) and reads a `build` variable
  that resolves by lexical accident to a stale package-level constant
  in zzz.r; the "3r quartile" typo appears twice.
- **Spec: FAIL (edges)** — the MAF values themselves are verified
  correct (hand computation matches all 235 comparable loci for a
  spot-checked population; the as.pop path works and restores). The
  defects: `@return` documents "An unaltered genlight object" but the
  function returns an invisible per-population MAF dataframe; the
  overall statistics print "No. of loci = nLoc(x)" while the
  statistics are computed on the monomorph-filtered x2 (mismatched
  denominators); the singleton-drop message reports the population
  INDEX ("23") rather than its name (verified); the bad-as.pop error
  directs users to loc.metrics when as.pop is checked against
  ind.metrics.

## Findings

### F1 — Verbosity contract broken (HIGH, confidence: certain)

Rule: VRB (the verbosity-override class). Location: per-pop tmpfun,
overall stats, quantile table, warnings.

Proposed change: delete the hardcoded `verbose = 3` (the captured
caller verbose then gates the per-pop stats at >= 3 as the code
structure intends); gate the overall statistics and quantile table at
verbose >= 1 (the ratified report-results gate); gate the
maf.limit/ind.limit warnings at >= 1 and the singleton notice at
>= 2.

### F2 — Duplicate Starting banner; accidental build reference (MEDIUM, confidence: certain)

Rule: FS/STY. Location: second FLAG SCRIPT START block.

utils.flag.start already emits the banner; the handwritten block
duplicates it and its `build` reference resolves to zzz.r's stale
package constant only by lexical accident. Proposed change: remove
the block.

### F3 — Singleton message content (MEDIUM, confidence: certain)

Rule: STY/spec. The dropped populations are reported by index, not
name, with a grammar slip. Proposed change: report names
(`names(drop_pop)`), fix grammar, gate at verbose >= 2 (with F1).

### F4 — Mismatched counts and wrong-slot error text (LOW, confidence: certain)

Rule: spec/DOC. The overall block reports nLoc(x)/nInd(x) while its
statistics come from the monomorph-filtered x2 — now reports the x2
counts actually used (with the unfiltered count alongside); the
bad-as.pop error names ind.metrics (was loc.metrics).

### F5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC5. `@return` corrected (invisible dataframe of
MAF by locus and population); `@param as.pop` "for the purposes of
deletions" corrected to reporting; verbose canon; "3r quantile" typo
x2; plot.dir/plot.file punctuation.

## Report notes (other functions / not fixed here)

- The per-population plots are computed for every population
  (filter+recalc each) including those later excluded by ind.limit —
  wasted computation, left as-is (structure change out of proportion).
- zzz.r's package-level `build = "v.2023.2"` constant is a footgun
  that masked F2's undefined-variable reference — worth removing when
  zzz.r is reviewed.

## Coverage

`tests/testthat/test-gl.report.maf.R` — 13 assertions: MAF values vs
hand computation, return structure, verbose-0 leak (baseline),
duplicate banner (baseline), as.pop reassign/restore + bad-as.pop
fatal, singleton drop by index (baseline), input untouched, invisible
return. All 13 pass on the pre-fix code.

## Approval

All five findings approved via the approval boxes (2026-08-31).

## Outcome

All five findings applied. Characterization suite: 13/13 pass; flips
map to F1 (silent at verbose 0; singleton notice at >= 2 with names),
F2 (single Starting banner), F3 (names not indices — the singletons in
testset.gl are EmmacNormLeic and EmmacNormSalt). MAF values unchanged
(still match the hand computation; as.pop reassign/restore intact).
End-to-end at verbose 3 clean with per-population statistics
displaying as structured. NEWS entry added.

```json
{"function": "gl.report.maf", "package": "dartR.base", "family_mode": "report",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["VRB"], "loc": "R/gl.report.maf.r gating", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["FS", "STY"], "loc": "R/gl.report.maf.r banner", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["STY", "spec"], "loc": "R/gl.report.maf.r singleton msg", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["spec", "DOC"], "loc": "R/gl.report.maf.r counts/error", "status": "applied"},
  {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC5"], "loc": "R/gl.report.maf.r header", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.report.maf.R", "pr": 294}
```
