# Review: gl.join (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.join; Datasets: testset.gl, testset.gs, platypus.gl
- Family mode: io/modify (object combination)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — many warnings are ungated (the method
  deprecation notice verified printing at verbose 0; missing-metrics
  and missing-flags warnings likewise ungated); two
  `cat(error()) + stop()` splits; an entire loc.metrics.flags block is
  TRIPLICATED — each join path sets the flags and then a third copy of
  the same block runs again for both paths (double warnings when flags
  are absent); end-of-function summary prints x1's individual count as
  the combined count; header drift.
- **Spec: FAIL** — two verified defects. (1) A join by shared loci
  LOSES ind.metrics entirely: plain rbind returns NULL metadata and
  gl.join never rebuilds it (verified: 7-row + 4-row inputs → NULL) —
  despite the comment claiming "rbind.dartR handles ind.metrics
  merging". (2) A join by shared individuals CRASHES on objects whose
  loc.metrics.flags data.frame lacks OneRatio/PIC columns — which is
  exactly testset.gl-style SNP data — with "replacement has 0 rows"
  (verified; platypus.gl happens to carry OneRatio, which is why the
  documented example works). Also: datatype1/datatype2 are checked
  individually but never compared, so a SNP + SilicoDArT pair with
  matching names would join silently; the description promises the
  history is cleared but it is inherited from x1 and appended
  (verified: 5 entries); @details documents method='sidebyside' /
  'end2end' values that have never existed in the code.

## Findings

### F1 — Loc-join loses ind.metrics (HIGH, confidence: certain)

Rule: DAT2/DAT3 (metadata desync). Proposed: after rbind, rebuild
`x@other$ind.metrics` by row-binding x1's and x2's metrics
(plyr::rbind.fill; plyr is in Imports), with the id column re-synced
after any make.unique.

### F2 — Ind-join flags crash on SNP-style flags (MEDIUM, confidence: certain)

Rule: DAT4. Proposed: multiply only flag fields present in both
objects' flags (iterate over the intersection of the known flag
names), instead of unconditionally assigning OneRatio/PIC products
that may be integer(0).

### F3 — Triplicated flags block; ungated warnings; stop splits (MEDIUM, confidence: certain)

Rule: STY/VRB. Remove the third flags block (per-path blocks remain —
behaviour unchanged, double warnings gone); gate all cat(warn())
notices at verbose >= 2; `cat(error()) + stop()` → `stop(error(...))`;
fix the end-of-function summary (wrong count, duplicated lines).

### F4 — Datatype mismatch not rejected; docs drift (LOW, confidence: certain)

Rules: DAT/DOC. Fatal error when datatype1 != datatype2. Docs: the
history is carried from x1 and appended, not cleared; @details
corrected to describe the actual auto-detection (and the real legacy
method values); "objectlacks" typos; stray "#'" in details; @return
one-liner expanded.

## Coverage

`tests/testthat/test-gl.join.R` — 17 assertions (loc-join and ind-join
happy paths, NULL-metadata baseline, flags-crash baseline, mismatch
fatal, legacy method baseline, duplicate-name handling, history call).
All pass pre-fix.

## Approval

All four findings approved via the approval boxes (2026-08-31).

## Addendum (discovered during apply)

Messages passed `substitute(x1)`/`substitute(x2)` directly to cat(),
which prints garbage at verbose 2 and crashes ("argument 2 (type
'language') cannot be handled by 'cat'") at verbose >= 3 whenever the
arguments are expressions rather than bare names (verified with
testset.gl[1:7, ]). Below-HIGH (display path); fixed within the
approved messaging scope by deparsing the argument expressions once.

## Outcome

All four findings applied, plus the addendum. Characterization suite:
20/20 pass; flipped assertions map to F1 (ind.metrics rebuilt, 11 rows,
id synced), F2 (testset-style ind-join no longer crashes), F3
(deprecation warning gated), F4 (mixed datatypes fatal). End-to-end at
verbose 3 clean on both join paths with correctly deparsed argument
names; after an ind-join the combined flags are exactly the 11 fields
shared by SNP-style inputs. NEWS entry added.

```json
{"function": "gl.join", "package": "dartR.base", "family_mode": "io",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["DAT2", "DAT3"], "loc": "R/gl.join.r loc path", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["DAT4"], "loc": "R/gl.join.r flags", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["STY", "VRB"], "loc": "R/gl.join.r flags/warnings", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["DAT", "DOC"], "loc": "R/gl.join.r", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs", "platypus.gl"],
 "baseline_test": "tests/testthat/test-gl.join.R", "pr": null}
```
