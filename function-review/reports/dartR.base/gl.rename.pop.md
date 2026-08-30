# Review: gl.rename.pop (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.rename.pop
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT), pop-less variant
- Family mode: modify (data manipulation)
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — visible return, redundant `is(x, "genlight")`
  check after the datatype check, unstyled progress message, header
  breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — renaming a nonexistent population is a silent no-op
  that still appends a history entry claiming the rename happened;
  renaming TO an existing name silently merges the two populations (an
  artefact of `levels<-` merge semantics; verified: 30 pops → 29, sizes
  summing); a pop-less object crashes with "attempt to set an attribute
  on NULL". The normal rename itself is exact (verified on both
  datatypes; input untouched; membership preserved).

## Findings

### F1 — Nonexistent population: silent no-op with a false history entry (MEDIUM, confidence: certain)

Rule: spec axis. Location: `R/gl.rename.pop.r:66`.

`levels(pop(x))[levels(pop(x)) == old] <- new` matches nothing when `old`
is not a level — the object is returned unchanged, no warning, and the
history records a rename that never happened (verified). Proposed
change: fatal error listing the populations present.

### F2 — Renaming to an existing name silently merges populations (MEDIUM, confidence: certain)

Rule: spec axis (silent API trap). Location: `R/gl.rename.pop.r:66`.

`levels<-` merges levels when the assigned value already exists —
verified: renaming EmmacRoss to EmmacTweeUki produced 29 populations with
the merged size equal to the sum of both. Undocumented; a typo can
silently collapse two populations. Options offered: document + warn, or
forbid. **Decision: forbid** — fatal error directing users to
gl.recode.pop for deliberate amalgamation.

### F3 — Pop-less object crashes cryptically (LOW, confidence: certain)

Rule: spec axis (input validation). Location: `R/gl.rename.pop.r:66`.

"attempt to set an attribute on NULL" (verified). Proposed change: the
standard pop-presence check with a clear fatal error (as in
gl.recode.pop).

### F4 — Redundant class check; unstyled message (LOW, confidence: certain)

Rule: STY. Location: `R/gl.rename.pop.r:53-58,60`.

`if (!is(x, "genlight")) stop(...)` is dead — utils.check.datatype has
already validated the object. The progress cat is not report()-styled.
Proposed change: remove the dead check; style the message.

### F5 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5. Location: `R/gl.rename.pop.r:77`. Proposed change:
`return(invisible(x))`.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.rename.pop.r:1-28`.

Description paragraphs belong in @details; `@return` after `@export`;
DOC2 wording variant; Custodian only; `@param x` says SNP only though
SilicoDArT works (verified); no @seealso. Proposed change: rewrite to
the ratified template, documenting the new existence/uniqueness rules.

## Report notes (other functions / not fixed here)

- The forbid decision makes gl.rename.pop strictly a 1:1 rename;
  gl.recode.pop (PR #266) is the documented amalgamation route.

## Coverage

Characterization baseline: `tests/testthat/test-gl.rename.pop.R` — 16
assertions pre-fix: exact rename on both datatypes, membership
preserved, history +1, input untouched, missing-argument errors, silent
no-op (baseline), silent merge (baseline), NULL-attribute crash
(baseline), visible return (baseline). All 16 passed pre-fix; the four
baseline assertions were flipped to the approved behaviours (14 pass
post-fix).

## Approval

Decisions recorded 2026-08-30 via approval boxes — all six findings
**approved**. F2 was a policy question (document + warn vs forbid vs keep):
**forbid** chosen — renaming to an existing name is now a fatal error.

## Outcome

Applied F1–F6 with F2 as forbid. Verification: all 14 characterization
assertions pass post-fix; every diff from baseline maps to an approved
finding (nonexistent pop → fatal error naming the valid pops [F1];
rename-to-existing → fatal error directing to gl.recode.pop [F2];
pop-less → clear error [F3]; invisible return [F5]). Normal renames
verified byte-identical membership on both datatypes. Caller grep across
dartR.base + 7 siblings: two tutorial calls renaming to fresh names —
unaffected. dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.rename.pop",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.rename.pop.r:66", "status": "approved"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.rename.pop.r:66", "status": "approved-forbid"},
    {"id": "F3", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.rename.pop.r:66", "status": "approved"},
    {"id": "F4", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.rename.pop.r:53-60", "status": "approved"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.rename.pop.r:77", "status": "approved"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.rename.pop.r:1-28", "status": "approved"}
  ],
  "datasets": ["testset.gl", "testset.gs", "pop-less variant"],
  "baseline_test": "tests/testthat/test-gl.rename.pop.R",
  "pr": null
}
```
