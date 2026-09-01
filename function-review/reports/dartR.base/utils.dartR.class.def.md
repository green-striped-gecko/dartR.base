# Review: utils.dartR.class.def (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.dartR.class.def.
- Datasets: testset.gl subsets, gen-backed and FBM-backed (via
  gl.gen2fbm); adegenet glSum/glMean as references.
- Family mode: S4 class infrastructure (the dartR class, its FBM
  extension, show/subset/cbind/rbind/glSum/glMean/NA.posi/initialize).
- Checks skipped: Google Group not searched (not available: no
  browser session).

## Verdicts

- **Standards: FAIL** — a dead helper (.fbmsub_copy) whose only call
  site is commented out and whose free variable `code` resolves to
  something invalid ("'code' must be of length 256" when invoked);
  the show method carries a `length(temp > 1)` luck idiom and a
  "lables" typo; cbind/rbind @return document "A genlight object"
  (they return dartR); an informal glSum doc ("not sure what
  happens ;-)"); getFromNamespace reaches into adegenet INTERNALS
  (.seppop_internal, .get_pop_inds) at load time.
- **Spec: FAIL (three defects; core verified sound)** — the verified
  core: XOR validity (@gen vs @fbm), subset/rbind/cbind round-trip
  EXACTLY on gen-backed objects; glSum/glMean match adegenet exactly;
  and the FBM layer matches the gen-backed reference exactly
  (conversion, subsetting, glSum/glMean, NA positions). The defects:
  (1) the show method's loc.metrics detail line is gated on
  IND.metrics being present — a copy-paste in the guard (verified:
  detail vanishes when ind.metrics is absent). (2) Character j
  subsetting with any unmatched locus name crashes with the cryptic
  "Cannot subset a SNPbin with mixed subscripts" (nomatch = 0 indices
  reach the SNPbin subsetter). (3) Negative indices on FBM-backed
  objects crash in big_copy (gen-backed objects handle them).

## Findings

### G1 — show-method copy-paste and idioms (MEDIUM, confidence: certain)

The loc.metrics detail guard tests ind.metrics; corrected to
loc.metrics. `length(temp > 1)` corrected to `length(temp) > 0`;
"lables" typo fixed.

### G2 — Unmatched locus names crash cryptically (MEDIUM, confidence: certain)

Character j is validated after matching: unmatched names produce an
informative stop naming the missing loci.

### G3 — Negative indices on FBM objects (MEDIUM, confidence: certain)

i/j numeric indices are normalized through seq_len()[i] before the
FBM branch, so negative indexing behaves as it does for gen-backed
objects.

### G4 — Dead .fbmsub_copy removed (MEDIUM, confidence: certain)

Superseded by bigstatsr::big_copy at the only (commented) call site;
carries a broken free-variable reference. Removed with its commented
call.

### G5 — Docs (LOW, confidence: certain)

cbind/rbind @return corrected to dartR; glSum's useC doc made
factual; minor tidy.

## Report notes (other functions / not fixed here) — GOVERNANCE

- **The duplicate-class message** ("Found more than one class 'dartR'
  in cache") originates at this file's setClass("dartR"): the retired
  monolithic dartR package defines an identical class name, so R's
  methods cache warns whenever both packages are loaded. No code fix
  here is correct: the remedy is ecosystem-level (users should not
  attach the retired dartR alongside dartR.base; longer term, the
  retired package's class definition goes away). Recorded for the
  team.
- getFromNamespace on adegenet internals is version-fragile
  (adegenet refactors would break load); noted, not changed.
- glSum/glMean are plain shadow functions: adegenet-INTERNAL
  S4-dispatched calls bypass them, so FBM-backed objects passed into
  adegenet internals would compute over the empty @gen. Inherent
  limitation of the shadow approach; FBM remains experimental.
- The SNPbin duplicate-index corruption (adegenet upstream) remains
  relevant to `[` with repeated indices.
- rbind's SNPbin fallback drops @other metadata (documented as lazy
  in its own docs); the FBM path preserves it. Inconsistent but
  documented.
- glSum's FBM path reads the FBM row-by-row (O(nInd) FBM reads);
  a column-block implementation would be faster. Performance note.

## Coverage

`tests/testthat/test-utils.dartR.class.def.R` — 17 assertions:
subset/rbind/cbind round-trips (anchor), glSum/glMean vs adegenet
(anchor), show copy-paste (baseline), character-j crash (baseline),
FBM-vs-gen equality suite (anchor), FBM negative-index crash
(baseline). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Suite: 17/17 pass; flips map to G1 (loc.metrics
detail shown without ind.metrics), G2 (informative stop naming the
unknown locus), G3 (FBM negative indexing works, nInd 9 on the
fixture). FULL package sweep after the subset-method change: 53
files, 847 assertions, only the nine pre-existing dev failures -
zero new. Subset-heavy callers (gl.filter.callrate, gl.keep.pop)
clean. PR #328.

```json
{"function": "utils.dartR.class.def", "package": "dartR.base", "family_mode": "infrastructure",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "G1", "severity": "MEDIUM", "rules": ["spec", "STY"], "loc": "R/utils.dartR.class.def.r show", "status": "applied"},
  {"id": "G2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.dartR.class.def.r subset char j", "status": "applied"},
  {"id": "G3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.dartR.class.def.r subset fbm", "status": "applied"},
  {"id": "G4", "severity": "MEDIUM", "rules": ["STY"], "loc": "R/utils.dartR.class.def.r .fbmsub_copy", "status": "applied"},
  {"id": "G5", "severity": "LOW", "rules": ["DOC"], "loc": "R/utils.dartR.class.def.r docs", "status": "applied"}],
 "datasets": ["testset.gl", "constructed FBM"],
 "baseline_test": "tests/testthat/test-utils.dartR.class.def.R",
 "pr": 328}
```
