# Review: gl.colors (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.colors
- Datasets: none required (pure palette getter; no genlight input)
- Family mode: graphics utility (read-only getter; no history, no
  datatype dispatch)
- Checks skipped: visual appearance not assessed (headless); Google
  Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — both error exits use the split
  `cat(error(...)) + stop(-1)` idiom, so the condition message an
  upstream tryCatch sees is literally "-1" while the real message leaks
  to stdout even inside try(silent = TRUE); the return is invisible,
  contradicting the function's own first example; header drift (missing
  @details, "@param verbose --", a `\item  4"` missing-quote typo, DOC2
  wording).
- **Spec: FAIL** — the description documents a `"pal"` category that
  the implementation rejects (verified: fatal error), while the
  `"structure"` type (35 colors, used by dartR.popgen::gl.plot.snmf) is
  implemented but entirely undocumented. Separately, the function is
  used as a lazy-evaluated default argument in 22 signature sites
  across 6 family packages, and each evaluation prints the 3-line
  banner regardless of the enclosing caller's verbose (verified; two
  sites already hand-patch `verbose = 0` around it).

## Findings

### F1 — cat + stop(-1) error mechanics (MEDIUM, confidence: certain)

Rule: STY/error conventions (the `cat(error())+stop()` split class).
Location: `R/gl.colors.r` (both validation exits).

`tryCatch(..., error = function(e) conditionMessage(e))` receives "-1";
the informative text goes to stdout where try(silent = TRUE) cannot
suppress it (verified). Proposed change: `stop(error("No valid color
option specified. Check ?gl.colors\n"))` at both exits.

### F2 — Banner leak at every default-argument evaluation (MEDIUM, confidence: certain)

Rule: VRB (verbose-0 leak class, lazy-default variant). Location:
signature/verbosity handling.

22 live signature sites across dartR.base, dartR.captive, dartR.popgen,
dartR.sim, and dartR.spatial use `gl.colors(...)` as a parameter
default. R evaluates those inside the callee, so gl.colors runs with
its own default verbosity (2) and prints "Starting gl.colors /
Selected color type / Completed:" no matter what verbose the user gave
the caller (verified). gl.report.fstat and gl.report.ld.map already
work around it with explicit `verbose = 0`. Proposed change: default
`verbose = 0` — a silent getter; banners only when explicitly
requested. Escalation-gate class: changes the parameter default —
interactive `gl.colors(2)` calls stop printing the banner, and a
session `gl.set.verbosity()` no longer applies unless verbose is
passed.

### F3 — "pal" documented but rejected; "structure" implemented but undocumented (LOW, docs-only, confidence: certain)

Rule: spec axis / DOC. Location: header.

`gl.colors("pal")` is a fatal error although the description lists it
as a category; `gl.colors("structure")` returns a 35-color palette
absent from the docs (and is relied on by dartR.popgen::gl.plot.snmf).
Proposed change: docs-only — drop the "pal" bullet, list "div", "dis",
"con", "vir" as the palette types directly, and document "structure".

### F4 — Invisible return contradicts the documented examples (LOW, confidence: certain)

Rule: spec axis / house product-utility convention. Location: return.

The first documented example is a bare `gl.colors(2)` — which displays
nothing, because the return is `invisible(cols)` (verified). The
sibling product-getters (gl.select.colors, gl.select.shapes) return
their vector visibly. Proposed change: `return(cols)`. User-visible
behavior change (bare calls start printing the palette).

### F5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7-adjacent. Missing @details (details prose sits
inside @description); `@param verbose --` stray dash and non-canon
wording; `\item  4"` missing opening quote; @return "returns colors as
a vector" is wrong for the four palette types, which return a function;
`cols <- NA` dead initialisation.

## Report notes (other functions / not fixed here)

- The 22 default-argument caller sites remain as they are; with F2
  rejected, each still leaks the banner unless its default passes
  `verbose = 0`. The two hand-patched sites (gl.report.fstat:322,
  gl.report.ld.map:325,333) show the correct pattern; the remaining 20
  can adopt it as each function comes up for its own review.
- dartR.popgen::gl.plot.snmf depends on the undocumented "structure"
  type — documenting it (F3) protects that dependency.

## Coverage

Characterization baseline: `tests/testthat/test-gl.colors.R` — 21
assertions: all four fixed palettes, the four palette functions'
outputs at n=4, structure length/first color, invalid-type errors
(bogus, 7, 1, "pal"), the "-1" condition message (baseline), the
3-line default-evaluation banner leak (baseline), invisible return
(baseline). All 21 pass on the pre-fix code.

## Approval

Decisions via the approval boxes (2026-08-30): F1 approved; F2
REJECTED — the verbose default stays 2; the default-argument banner
leak is retained deliberately and callers must pass verbose = 0; F3
approved; F4 approved (visible return); F5 approved.

## Outcome

F1, F3, F4, F5 applied; F2 not applied per rejection. Characterization
suite: 22/22 pass; the two flipped assertions map to F1 (condition
message now "No valid color option...", nothing to stdout at verbose 0)
and F4 (visible return); the banner-leak assertion is retained
unchanged as deliberate behavior. End-to-end at verbose 3 clean on
fixed, structure, and palette types. Caller impact of F4: none — the
visible return only affects bare top-level calls; the 22
default-argument sites never auto-print. NEWS entry added.

```json
{
  "function": "gl.colors",
  "package": "dartR.base",
  "family_mode": "graphics-utility",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["STY"], "loc": "R/gl.colors.r validation exits", "status": "applied"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB", "API2"], "loc": "R/gl.colors.r signature", "status": "rejected"},
    {"id": "F3", "severity": "LOW", "rules": ["DOC"], "loc": "R/gl.colors.r header", "status": "applied"},
    {"id": "F4", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.colors.r return", "status": "applied"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2"], "loc": "R/gl.colors.r header", "status": "applied"}
  ],
  "datasets": [],
  "baseline_test": "tests/testthat/test-gl.colors.R",
  "pr": 276
}
```
