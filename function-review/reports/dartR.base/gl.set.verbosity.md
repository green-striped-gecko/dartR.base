# Review: gl.set.verbosity (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: ddaed27 (dev, synced with upstream/dev);
  Branch: review-gl.set.verbosity; Datasets: none required
  (environment setter)
- Family mode: utility (global-option setter; partner of
  gl.check.verbosity)
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — validation happens after utils.flag.start has
  already consumed the raw value (so an out-of-range value triggers
  gl.check.verbosity's coercion warning on the banner path before the
  setter's own logic runs); no @family tag (its partner carries
  "environment"); the example misleadingly assigns the invisible NULL
  (`gl <- gl.set.verbosity(...)`).
- **Spec: FAIL** — invalid values are a silent no-op that CLAIMS
  success: `gl.set.verbosity(7)` leaves the option untouched and then
  prints "Global verbosity set to: 7" (verified); a character value
  rides string comparisons through the echo; `gl.set.verbosity(NULL)`
  crashes with "argument is of length zero" (verified — the `&` chain
  yields logical(0) in the if). The valid path is verified correct
  (option set, picked up by gl.check.verbosity, silent at 0); @return
  promises "verbosity value" while the code returns invisible(NULL).

## Findings

### F1 — Invalid values: false success, garbage echo, NULL crash (MEDIUM, confidence: certain)

Rule: spec axis. Location: validation and messaging.

Proposed change: validate upfront (before flag.start): NULL,
non-numeric, or out-of-range values warn and coerce to the default 2
(the gl.check.verbosity coercion pattern), the option is then actually
set, and the confirmation message reports the value actually set.

### F2 — Return value vs contract (LOW, confidence: certain)

Rule: DOC5/API. @return promises the value; invisible(NULL) is
returned. Proposed change: return the value actually set, invisibly,
and document it.

### F3 — Header and ordering (LOW, confidence: certain)

Rule: DOC/FS. Validation moved before flag.start; @family environment
added; the example no longer assigns; verbose param prose canon.

## Report notes (other functions / not fixed here)

- gl.check.verbosity (Bernd): a vector argument (e.g. c(2,3)) reaches
  `if` with a length-2 condition — an error in current R; and its
  coercion warning is ungated (defensible for an explicit user
  argument). For its own review.

## Coverage

`tests/testthat/test-gl.set.verbosity.R` — 9 assertions: valid set +
pickup + silent-at-0, NULL-return baseline, false-success baseline
(option untouched, "set to: 7" printed), NULL crash baseline. All 9
pass on the pre-fix code. Global option restored by every test.

## Approval

All three findings approved via the approval boxes (2026-08-31).

## Outcome

All three findings applied. Suite: 12/12 pass; flips map to F1
(invalid values including NULL warn, coerce to 2, and genuinely set
it; the confirmation reports the value actually set) and F2 (the set
value returned invisibly). Valid-path behaviour unchanged (option set
and picked up by gl.check.verbosity; silent at 0). NEWS entry added.

```json
{"function": "gl.set.verbosity", "package": "dartR.base", "family_mode": "utility",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.set.verbosity.r validation", "status": "applied"},
  {"id": "F2", "severity": "LOW", "rules": ["DOC5", "API"], "loc": "R/gl.set.verbosity.r return", "status": "applied"},
  {"id": "F3", "severity": "LOW", "rules": ["DOC", "FS"], "loc": "R/gl.set.verbosity.r header", "status": "applied"}],
 "datasets": [],
 "baseline_test": "tests/testthat/test-gl.set.verbosity.R", "pr": 297}
```
