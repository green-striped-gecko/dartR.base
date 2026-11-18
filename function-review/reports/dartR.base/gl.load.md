# Review: gl.load (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.load; Datasets: testset.gl (round trip via gl.save)
- Family mode: io (RDS load wrapper)
- Checks skipped: fbm = TRUE conversion path not exercised (no FBM
  session configured); Google Group not searched (not available: no
  browser session).

## Verdicts

- **Standards: FAIL** — the "Loaded object" message and the
  fbm-conversion message are ungated (print at verbose 0);
  gl.compliance.check is called with no verbose argument, so its full
  chatter (21 lines total) prints at verbose 0 (verified); the
  dartR-conversion notice is gated `verbose > 2` against the house
  `>= 2`; @return after @export; no @examples.
- **Spec: FAIL** — the `compliance` parameter, documented
  "[default FALSE]", is inert: gl.compliance.check runs
  unconditionally regardless of the setting (verified: TRUE and FALSE
  produce byte-identical output and objects). The member directed at
  nomination that the default behaviour be no compliance check — i.e.
  make the documented default real. Also: @param fbm documents
  "[default TRUE]" but the signature default is FALSE; @param file
  says "file to receive data" (it is the file to read); the
  description says the object is loaded "from the current workspace"
  (it is read from file); a non-genlight RDS produces the cryptic
  error "no applicable method for `@` applied to an object of class
  'dartR'" via the blind class assignment; a missing file surfaces a
  raw connection error. Round-trip integrity verified (genotypes and
  dimensions identical through gl.save/gl.load); invisible return
  correct.

## Findings

### F1 — compliance parameter inert (MEDIUM; DIRECTED, confidence: certain)

Rule: spec axis (inert documented parameter). Location: compliance
call.

Directed by the member at nomination: "make the default for compliance
check to be no compliance check." Applied as: the check runs only when
`compliance = TRUE`, with verbose passed through. Consequence: objects
loaded with the default are no longer silently passed through
gl.compliance.check (which can modify flags/metrics); users wanting
the old behaviour pass compliance = TRUE.

### F2 — Verbose-0 leaks (MEDIUM, confidence: certain)

Rule: VRB. The "Loaded object" and fbm-conversion messages gated at
verbose >= 2; the compliance call receives the caller's verbose; the
dartR-conversion notice gate aligned to >= 2.

### F3 — File and object validation (MEDIUM, confidence: certain)

Rule: spec/DAT. A missing file now gives a clear fatal error before
readRDS; an RDS that does not contain a genlight object gives a clear
fatal error instead of the cryptic "@"-method failure.

### F4 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC. @param fbm corrected to [default FALSE]; @param file
corrected (file to read); description corrected (loads from file);
@return before @export with the invisible return noted; an @examples
block added (tempdir save/load round trip).

## Report notes (other functions / not fixed here)

- gl.save: its "Saved object" / "Load again" messages are ungated
  (print at verbose 0), the fbm-conversion message likewise, and the
  text says "RDA file" for an RDS file — for gl.save's own review.

## Coverage

`tests/testthat/test-gl.load.R` — 12 assertions: round-trip integrity
(genotypes, dimensions, pops), invisible return, verbose-0 leak
baseline, inert-compliance baseline (TRUE vs FALSE identical),
compliance-banner-always baseline, missing-file and non-genlight
errors. All 12 pass on the pre-fix code.

## Approval

F1 directed by the member at nomination ("make the default for
compliance check to be no compliance check"); F2, F3, F4 approved via
the approval boxes (2026-08-31).

## Outcome

All four findings applied. Characterization suite: 12/12 pass; flips
map to F1+F2 (silent at verbose 0; compliance banner only when
compliance = TRUE), and F3 (clear "not found" / "genlight" errors).
Round-trip integrity unchanged. End-to-end at verbose 3 clean with and
without compliance = TRUE. NEWS entry added.

```json
{"function": "gl.load", "package": "dartR.base", "family_mode": "io",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.load.r compliance", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.load.r messages", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.load.r validation", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["DOC1", "DOC"], "loc": "R/gl.load.r header", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.load.R", "pr": 289}
```
