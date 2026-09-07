# Review: gl.filter.locmetric (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.locmetric
- Datasets: testset.gl with a dummy metric `test = 1:nLoc` (as in the
  function's own example), NA-injected and pop-less variants
- Family mode: modify
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — visible return, cryptic failure on an invalid
  `keep` value, missing is.numeric validation (the report sibling has
  it), irrelevant preambles (a full monomorphs scan for a warning; a
  pop-assignment block that stamps 'pop1' onto pop-less objects), header
  breaks DOC1/DOC2/DOC7 and lacks @seealso to the matched report.
- **Spec: FAIL** — `keep = "outside"` has never worked: the condition is
  `field <= lower AND field >= upper`, impossible whenever lower < upper,
  so the subset is always empty and every call crashes ("Subsetting
  resulted in zero loci"). The documented mode is unreachable. 'within'
  mode is numerically exact (56 = independent recomputation, boundaries
  inclusive, metadata in sync, one history entry).

## Findings

### F1 — keep='outside' always crashes; the mode has never worked (HIGH, confidence: certain)

Rule: spec axis (documented mode unreachable). Location:
`R/gl.filter.locmetric.r:136-144`.

`which(field <= lower & field >= upper)` requires a value simultaneously
at-or-below lower AND at-or-above upper — impossible for lower < upper.
The index is always empty, `x[, integer(0)]` throws "Subsetting resulted
in zero loci". Verified on the function's own example data.

Proposed change: `which(field < lower | field > upper)` — the exact
complement of 'within' (which keeps [lower, upper] inclusive), so every
locus with a non-NA metric lands in exactly one of the two modes.
Escalation-gate class: a mode that always crashed becomes functional.

### F2 — Invalid keep value crashes cryptically; no numeric validation (MEDIUM, confidence: certain)

Rule: spec axis (input validation). Location:
`R/gl.filter.locmetric.r:66-71,119-121,127`.

- `keep = "bogus"` → neither branch assigns `x2` → "object 'x2' not
  found" (verified) — meaningless to the user. Proposed: warn (gated at
  verbose >= 1) and coerce to 'within', matching the
  gl.filter.secondaries method-coercion pattern.
- A non-numeric metric (e.g. AlleleID) produces factor-comparison
  warnings then the zero-loci crash (verified). Proposed: add the report
  sibling's `is.numeric` check as a fatal error.

### F3 — Irrelevant preambles: monomorphs scan and pop stamping (MEDIUM, confidence: certain)

Rule: spec axis / STY (function does more than its name); performance.
Location: `R/gl.filter.locmetric.r:95-115`.

- `tmp <- gl.filter.monomorphs(x, verbose = 0)` runs a full O(nLoc×nInd)
  monomorph scan on every call solely to warn (at verbose >= 2) that
  monomorphic loci exist — irrelevant to metric filtering and a real
  cost on large datasets.
- The pop-assignment block stamps pop = 'pop1' onto every individual of
  a pop-less object (verified) — a silent data modification unrelated to
  the filter's job; no other reviewed filter does this and nothing in
  this function uses pop.

Proposed change: delete both blocks. Escalation-gate class: the
monomorph warning disappears, and pop-less objects come back with their
pop state unchanged instead of 'pop1'.

### F4 — NA-metric loci removed but not itemised (LOW, confidence: certain)

Rule: DAT NA lens; convention from PRs #256-#262. Location:
`R/gl.filter.locmetric.r:127-144`.

`which()` already drops NA comparisons, so NA-metric loci are removed
with sync intact (verified: 46 vs 56 with 10 NAs) — the campaign's
NA policy, already satisfied. Proposed change: itemise the NA count in
the verbose >= 3 summary only (no behaviour change).

### F5 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent). Location:
`R/gl.filter.locmetric.r:162`. Proposed change: `return(invisible(x2))`.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.filter.locmetric.r:1-64`.

- `@return` after `@export` (DOC1); DOC2 verbose wording stale.
- `@author` "Luis Mijangos -- Post to …" — add Author(s)/Custodian
  labels (as done for the report sibling in PR #263).
- `@param keep` documented as "[within]" → "[default 'within']"; wording
  garbled ("keep loci within of upper and lower thresholds").
- No `@seealso` to gl.report.locmetric (noted in PR #263's review).
- Document the boundary semantics: 'within' keeps [lower, upper]
  inclusive; 'outside' keeps the strict complement.

Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- The "bug in adegenet" guard (nLoc vs metrics rows) is a legitimate
  integrity check and is retained.
- gl.filter.monomorphs itself is fine (reviewed, PR #251) — the issue in
  F3 is calling it here at all.

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.locmetric.R` —
15 assertions: 'within' counts vs independent recomputation + inclusive
boundary + sync, history +1, input untouched, metric-absent error, NA
removal with sync, 'outside' always-crash (baseline), cryptic invalid-keep
crash (baseline), non-numeric crash-with-warnings (baseline), pop1
stamping (baseline), visible return (baseline). All 15 pass on the
pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all six findings
**approved**, including the escalation-gate consequences (F1: 'outside'
becomes functional; F3: monomorph warning and pop1 stamping removed).

## Outcome

Applied F1–F6. Verification: all 19 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding
('outside' retains the exact complement of 'within' — 46 + 199 + 10 NA =
255, a perfect partition [F1]; invalid keep coerces with a gated warning
and non-numeric metrics get a clear fatal error [F2]; pop-less objects
keep their pop state and no monomorph scan runs [F3]; NA removals
itemised [F4]; invisible return [F5]). End-to-end at verbose = 3 in both
modes with NA-injected metrics. Caller grep across dartR.base + 7
siblings: only the report sibling's docs reference this function — no
code caller. dartr2shiny: not present in the workspace. NEWS entry
added.

```json
{
  "function": "gl.filter.locmetric",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.filter.locmetric.r:136-144", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.locmetric.r:66-71,119-127", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["spec", "STY"], "loc": "R/gl.filter.locmetric.r:95-115", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["DAT"], "loc": "R/gl.filter.locmetric.r:127-144", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.locmetric.r:162", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.locmetric.r:1-64", "status": "proposed"}
  ],
  "datasets": ["testset.gl (dummy metric)", "NA-injected variant", "pop-less variant"],
  "baseline_test": "tests/testthat/test-gl.filter.locmetric.R",
  "pr": 264
}
```
