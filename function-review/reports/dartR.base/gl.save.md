# Review: gl.save (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.save; Datasets: testset.gl
- Family mode: io (RDS save wrapper)
- Checks skipped: FBM input path not exercised (no FBM session
  configured); Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the "Saved object" and "Load again" messages
  are ungated (2 lines at verbose 0, verified), as is the
  FBM-conversion message; the message calls the RDS file "a compressed
  RDA file"; @return after @export; description says the object is
  saved "to the current workspace" (it is written to file).
- **Spec: FAIL (minor)** — @return promises "The input object", but
  the returned object is modified: the class-attribute stamping
  (`attributes(class(x)) <- list(package = "dartR.base")`) is applied
  before the return (verified: not identical to the input), and for
  FBM inputs the returned object would additionally be the gen-format
  conversion made for saving. The saved file itself is verified
  correct (loadable; genotypes and dimensions identical).

## Findings

### F1 — Verbose-0 message leaks; RDA wording (MEDIUM, confidence: certain)

Rule: VRB. The two save messages and the FBM-conversion message gated
at verbose >= 2; "RDA" corrected to RDS.

### F2 — Returned object is not the input (LOW-MEDIUM, confidence: certain)

Rule: spec/API (@return contract). The input is held at entry and
returned unchanged; the class-attribute stamping (and any FBM→gen
conversion) applies only to the copy that is saved. Mild behaviour
change: pipelines relying on the stamped/converted return value would
see the true input instead — none found in the family (no functional
callers).

### F3 — Path validation and header (LOW, confidence: certain)

Rule: spec/DOC. A nonexistent target directory now gives a clear fatal
error before saveRDS (was a raw connection error); description
corrected ("to file"); @return before @export with the invisible
return noted.

## Coverage

`tests/testthat/test-gl.save.R` — 7 assertions: file written and
loadable with identical genotypes, invisible return, verbose-0 leak
baseline, returned-object-modified baseline, bad-path error. All 7
pass on the pre-fix code.

## Approval

All three findings approved via the approval boxes (2026-08-31).

## Outcome

All three findings applied. Characterization suite: 7/7 pass; flips
map to F1 (silent at verbose 0), F2 (returned object identical to the
input), F3 (clear directory error). Round trip verified: the saved
copy retains the stamped dartR class and loads via gl.load with
identical genotypes. NEWS entry added.

```json
{"function": "gl.save", "package": "dartR.base", "family_mode": "io",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.save.r messages", "status": "applied"},
  {"id": "F2", "severity": "LOW", "rules": ["spec", "API"], "loc": "R/gl.save.r return", "status": "applied"},
  {"id": "F3", "severity": "LOW", "rules": ["spec", "DOC"], "loc": "R/gl.save.r", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.save.R", "pr": null}
```
