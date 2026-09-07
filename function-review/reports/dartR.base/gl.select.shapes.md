# Review: gl.select.shapes (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.select.shapes
- Datasets: testset.gl (for the x= path)
- Family mode: graphics utility (interactive; visible vector return is
  the product and judged correct)
- Checks skipped: visual appearance not assessed (headless); Google
  Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — the plot draws unconditionally (no plot.display
  parameter, no verbose gate — its sibling gl.select.colors gates at
  verbose >= 1 && plot.display); the datatype banner prints at
  verbose = 0 when x is supplied; the genlight argument `x` is shadowed
  by the plot x-coordinates mid-function; "Requires shapes" typo;
  DOC1/DOC2/DOC7 drift.
- **Spec: FAIL** — the range validation is a parenthesis slip:
  `min(select < 0 | max(select > 25))` fires only when EVERY element is
  negative, so partially-negative selects (verified: c(-1, 5)) slip
  through to pch; and the x= argument is entirely non-functional — nPop
  is computed into `nshapes` and then discarded, so a select of the
  wrong length passes silently and a NULL select returns all 26 shapes
  regardless of nPop (both verified). Valid inputs behave correctly
  (0:25 default, given select returned as-is, all-negative and >25
  caught).

## Findings

### F1 — Range validation parenthesis slip (MEDIUM, confidence: certain)

Rule: spec axis (the `length(v > 0)` misplaced-parenthesis class).
Location: `R/gl.select.shapes.r` (select validation).

`if (min(select < 0 | max(select > 25)))` evaluates a logical vector's
min: it is 1 (error fires) only when all elements are negative or any
exceeds 25. Verified: select = c(-1, 5) passes and the negative pch
reaches plot(). Proposed change:
`if (min(select) < 0 | max(select) > 25)`.

### F2 — x= argument is non-functional (MEDIUM, confidence: certain)

Rule: spec axis (documented behaviour absent). Location: x= block.

nshapes is set to nPop(x) and never used; the "must agree in number"
message is informational only. Verified: select of length 2 with 30
populations passes silently; NULL select with x returns 26 shapes.
Proposed change: mirror gl.select.colors — with select supplied, a
length mismatch with nPop(x) is a fatal error; with select NULL, the
first nPop(x) shapes (0..nPop-1) are displayed and returned.
Escalation-gate class: x= calls acquire real validation and
nPop-driven defaults.

### F3 — Plot unconditional; no plot.display parameter (MEDIUM, confidence: certain)

Rule: VRB/PLT (house convention). Location: plotting block.

plot() runs at every verbosity with no way to suppress it. Proposed
change: add `plot.display = TRUE` and gate the display with
`verbose >= 1 && plot.display` (exactly the sibling's pattern). API
addition, backward compatible.

### F4 — Datatype banner at verbose 0; variable shadowing; idioms (LOW, confidence: certain)

`utils.check.datatype(x)` without verbose prints the banner at
verbose = 0 (verified); the genlight `x` is overwritten by plot
x-coordinates (works only because nPop is read first — renamed);
`c(seq(1:26) - 1)` → `0:25`; "Requires shapes" → "Required shapes".

### F5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. `@return` after `@export`; DOC2 wording;
Custodian only; the new plot.display parameter documented.

## Report notes (other functions / not fixed here)

- The visible vector return is the interactive product — correct as-is
  (consistent with gl.select.colors, PR #274).

## Coverage

Characterization baseline: `tests/testthat/test-gl.select.shapes.R` —
11 assertions: default 0:25 + silence, valid select passthrough,
all-negative and >25 caught, partial-negative slip (baseline), x=
non-functionality + banner leak (baseline), input untouched. All 11
pass on the pre-fix code.

## Approval

All five findings approved via the approval boxes (2026-08-30):
F1 approved; F2 approved as recommended (fatal on length mismatch,
nPop-driven default); F3 approved as recommended (plot.display added,
chart gated at verbose >= 1); F4 approved; F5 approved.

## Addendum (discovered during apply)

The nPop-driven default for F2 cannot apply when nPop(x) > 26 — only 26
distinct pch shapes exist (testset.gl itself has 30 populations). A NULL
select with more than 26 populations is now a fatal error directing the
user to specify select manually (repeats allowed). Below-HIGH, recorded
here per the flow rule; covered by the baseline test.

## Outcome

All five findings applied. Characterization suite: 16/16 pass; every
diff from baseline maps to an approved finding (partial-negative select
now fatal [F1]; x= length mismatch fatal, NULL select with x returns
one shape per population, >26 populations fatal [F2 + addendum]; chart
suppressible and off at verbose 0 [F3]; banner silent at verbose 0
[F4]). End-to-end at verbose 3 clean on default, select+x, and
NULL-select+x paths. Caller grep all-clear: no functional callers of
gl.select.shapes in dartR.base, the 7 sibling packages, or dartRstartup
— documentation mentions only (gl.pcoa.plot, gl.select.colors,
dartR.captive::gl.grm.network). NEWS entry added.

```json
{
  "function": "gl.select.shapes",
  "package": "dartR.base",
  "family_mode": "graphics-utility",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.select.shapes.r validation", "status": "applied"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.select.shapes.r x= block", "status": "applied"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["VRB", "PLT"], "loc": "R/gl.select.shapes.r plotting", "status": "applied"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB2", "STY"], "loc": "R/gl.select.shapes.r", "status": "applied"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.select.shapes.r header", "status": "applied"}
  ],
  "datasets": ["testset.gl"],
  "baseline_test": "tests/testthat/test-gl.select.shapes.R",
  "pr": 275
}
```
