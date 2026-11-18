# Review: gl.report.basics (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.report.basics; Datasets: testset.gl, testset.gs,
  constructed 20x30 genlight with an all-NA individual
- Family mode: report (console summary only; returns invisible NULL —
  input untouched verified, no history append, correct)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the all-NA individuals listing prints the
  entire NA-padded array ("NA NA ind3 NA NA ..." — verified; the
  gl.report.allna class, PR #250); the missing-rdepth path emits a raw
  mean.default warning and prints NA; the datatype banner is hardcoded
  to verbose 0 rather than passing verbose through; per-locus for-loops
  and two separate per-pop gl.keep.pop loops where one vectorized pass
  serves; a commented-out `# @import tibble` line (tibble is in
  Imports).
- **Spec: FAIL** — the function has NEVER worked on SilicoDArT data:
  the composition table hard-codes four column names
  (`c("0","1","2","NA")`) onto `table(as.matrix(x), useNA="always")`,
  which has three classes for presence/absence data — fatal "length of
  'dimnames' [2] not equal to array extent" (verified). The same crash
  hits SNP data in which any genotype class is entirely absent
  (verified on a 10-locus subset with no homozygote-alternate scores).
  The SilicoDArT monomorph branch is unreachable in practice. Verified
  correct on standard SNP data: dimensions, monomorph count exactly
  matches gl.filter.monomorphs (144 on testset.gl), all-NA loci count,
  per-pop tables; silent at verbose 0 with invisible NULL.

## Findings

### F1 — Composition table crashes on SilicoDArT and class-absent SNP data (HIGH, confidence: certain)

Rule: spec/DAT dispatch. Location: composition block.

Proposed change: tabulate over explicit factor levels per datatype
(`c("0","1")` for SilicoDArT, `c("0","1","2")` for SNP) with an NA
column, so the table always has the expected shape. Restores the
documented SilicoDArT support; no change to SNP output values.

### F2 — All-NA individuals listing prints the NA-padded array (MEDIUM, confidence: certain)

Rule: STY (the #250 class). Location: all-NA individuals block.

Proposed change: list only the names of all-NA individuals
(`na.omit`), preserving the count line.

### F3 — Average Read Depth unguarded (MEDIUM, confidence: certain)

Rule: DAT/VRB. Location: read-depth block.

Objects without an rdepth locus metric (including testset.gs) trigger
a raw "argument is not numeric or logical" warning and print NA
(verified). Proposed change: print "not available" when the metric is
absent; otherwise `round(mean(rdepth, na.rm = TRUE), 2)`.

### F4 — Style, efficiency, and header (LOW, confidence: certain)

Rules: STY, DOC. Per-locus for-loops for the monomorph and all-NA
counts vectorized; the two per-pop gl.keep.pop loops merged into one
pass computing both per-pop statistics; datatype check passes verbose
through; commented-out `# @import tibble` removed; `@return` states
that results are printed and NULL is returned invisibly.

## Report notes (other functions / not fixed here)

- The monomorph count treats all-NA loci as monomorphic — exactly
  matching gl.filter.monomorphs' behaviour (144 = 144 on testset.gl),
  so it is reported here as consistent, not a defect.

## Coverage

Characterization baseline: `tests/testthat/test-gl.report.basics.R` —
14 assertions: verbose-0 silence/invisible NULL/input untouched, key
SNP statistics (dimensions, monomorph count vs gl.filter.monomorphs,
all-NA loci), SilicoDArT crash (baseline), class-absent SNP crash
(baseline), NA-padded listing (baseline), missing-rdepth warning
(baseline). All 14 pass on the pre-fix code.

## Approval

All four findings approved via the approval boxes (2026-08-31).

## Outcome

All four findings applied. Characterization suite: 17/17 pass; flipped
assertions map to F1 (SilicoDArT and class-absent SNP data now run),
F2 (names-only listing), F3 (no warning, "not available"). SNP key
statistics unchanged (dimensions, monomorph count still 144 matching
gl.filter.monomorphs, all-NA loci 3). End-to-end at verbose 3 clean on
SilicoDArT (0/1/NA composition table) and SNP data. No functional
callers of gl.report.basics in the family. NEWS entry added.

```json
{"function": "gl.report.basics", "package": "dartR.base", "family_mode": "report",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec", "DAT"], "loc": "R/gl.report.basics.r composition", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["STY"], "loc": "R/gl.report.basics.r allNA listing", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["DAT", "VRB"], "loc": "R/gl.report.basics.r rdepth", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["STY", "DOC"], "loc": "R/gl.report.basics.r", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs", "constructed"],
 "baseline_test": "tests/testthat/test-gl.report.basics.R", "pr": 285}
```
