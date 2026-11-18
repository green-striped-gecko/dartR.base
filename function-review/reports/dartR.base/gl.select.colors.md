# Review: gl.select.colors (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.select.colors
- Datasets: testset.gl (for the x= path); palette libraries probed
  directly
- Family mode: graphics utility (interactive; visible vector return is
  the product and judged correct)
- Checks skipped: visual appearance of palettes not assessed (headless);
  Google Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — missing-package checks return -1; the datatype
  banner prints at verbose = 0 when x is supplied
  (utils.check.datatype called without verbose); header typos
  ("hue_pl", "select bu number"), DOC2/DOC7 drift, 12-line dead
  "Testing scripts" comment block.
- **Spec: FAIL** — an invalid `library` value silently returns base R's
  colors() FUNCTION (the internal variable is never assigned and
  lexical scoping finds grDevices::colors); requested colour counts are
  not honoured for brewer (2 requested → 3 delivered; 12 requested from
  Blues → 9 delivered, with only raw brewer.pal warnings);
  out-of-bounds `select` indices yield NA colours; and baseR
  palette='heat' passes validation but dispatches to rainbow with a
  confusing warning. The internal-caller contract used across the
  reviewed report/filter functions (brewer, Blues, select = c(7,5),
  verbose = 0 → "#2171B5" "#6BAED6", silent) was verified intact and is
  protected by the baseline.

## Findings

### F1 — Invalid library returns the colors() function (HIGH, confidence: certain)

Location: `R/gl.select.colors.r` (library dispatch).

`library = "bogus"` matches no branch, `colors` is never assigned, and
`return(colors)` resolves to grDevices::colors — a function object is
returned as the "colour vector" (verified). Downstream plot code
receives garbage with no error at the source. Proposed change: validate
`library` up front against c('scales','brewer','gr.palette','gr.hcl',
'baseR'); unknown values coerce to the default (scales/hue_pal) with a
warning gated at verbose >= 1 — the same coercion pattern the function
already uses for unknown palettes.

### F2 — Requested colour count not honoured (brewer) (MEDIUM, confidence: certain)

Location: brewer branch.

brewer.pal cannot return fewer than 3 or more than the palette maximum:
ncolors = 2 delivers 3, ncolors = 12 from Blues delivers 9 (verified),
signalled only by raw R warnings. Proposed change: for requests below
3, draw brewer.pal(3, palette) and trim to the requested count (the
request is then honoured exactly); for requests above the palette
maximum, deliver the maximum with a clear warning gated at
verbose >= 1 stating requested vs delivered.

### F3 — select indices unvalidated: NA colours (MEDIUM, confidence: certain)

Location: select application.

`select = c(1, 15)` against 9 colours returns an NA entry (verified) —
invisible until a plot fails downstream. Proposed change: fatal error
when any select index is outside 1..ncolors.

### F4 — baseR 'heat' validated but not dispatched (MEDIUM, confidence: certain)

Location: baseR branch.

'heat' is in the validation list but has no dispatch branch; it falls
to the else, returning rainbow with the warning "nominated palette not
in Base R" (verified). Proposed change: add the heat.colors branch —
palette = 'heat' returns heat.colors(ncolors).

### F5 — Datatype banner at verbose 0; return(-1) package checks (LOW, confidence: certain)

`utils.check.datatype(x)` is called without verbose, so "Processing
genlight object…" prints even at verbose = 0 when x is supplied
(verified); the two package checks `cat(error); return(-1)` from a
function documented to return a colour vector. Proposed change: pass
`verbose = verbose`; `stop(error(...))`.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Typos "scales::hue_pl" and "select bu number";
DOC2 wording; Custodian only; `@return` after the dead 12-line "Testing
scripts" comment block (removed); the verbose >= 1 gate on the palette
display is house-consistent and documented.

## Report notes (other functions / not fixed here)

- The internal-caller contract (brewer/Blues/select 7,5/verbose 0) is
  asserted first in the baseline and unaffected by every proposed fix.
- gl.colors (the sibling fixed-palette helper with the verbose-0 banner
  leak noted in the hwe review) is a separate function, still to be
  reviewed.

## Coverage

Characterization baseline: `tests/testthat/test-gl.select.colors.R` —
13 assertions: internal-caller contract, default path, x= driving
ncolors + select-length validation, function-return bug (baseline),
brewer count mismatches (baseline), NA select (baseline), heat→rainbow
(baseline), datatype banner leak (baseline). All 13 pass on the pre-fix
code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all six findings
**approved** (F1 as coercion-to-default, the recommended option, over
the fatal-error alternative).

## Outcome

Applied F1–F6. Verification: all 17 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (unknown
library coerces to scales/hue_pal, silent at verbose 0 [F1]; 2 requested
= 2 delivered, above-max warns gated and returns the maximum [F2];
out-of-bounds select is a clear fatal error [F3]; palette='heat' returns
heat.colors with no confusing warning [F4]; x= path silent at verbose 0
and missing packages stop [F5]). Caller grep: 29 caller files across
the family in four call patterns, all verified compatible — the
package-wide default-colour contract (brewer, Blues, select = c(7,5),
verbose = 0 → "#2171B5" "#6BAED6", silent) is the first assertion in
the baseline and unchanged. dartr2shiny: not present in the workspace.
NEWS entry added.

```json
{
  "function": "gl.select.colors",
  "package": "dartR.base",
  "family_mode": "graphics-utility",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.select.colors.r library dispatch", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.select.colors.r brewer branch", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.select.colors.r select application", "status": "proposed"},
    {"id": "F4", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.select.colors.r baseR branch", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB2", "spec"], "loc": "R/gl.select.colors.r checks", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.select.colors.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl"],
  "baseline_test": "tests/testthat/test-gl.select.colors.R",
  "pr": 274
}
```
