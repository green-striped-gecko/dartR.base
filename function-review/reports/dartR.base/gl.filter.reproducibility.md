# Review: gl.filter.reproducibility (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.reproducibility
- Datasets: testset.gl (SNP, RepAvg), testset.gs (SilicoDArT,
  Reproducibility)
- Family mode: modify
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — internal gl.drop.loc delegation leaks a second
  history entry (FS8), threshold warning prints at verbose = 0 (VRB2),
  visible return, missing `verbose == 0 → plot.display = FALSE` guard,
  indented `@family` tag leaking into the rendered help title, header
  breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — a dataset missing the RepAvg/Reproducibility metric is
  returned UNFILTERED with no error or warning (silent no-op), while the
  function demands AlleleID/CloneID which it never uses; `plot.file` with
  `plot.display = FALSE` crashes; NA-metric loci are silently retained,
  contradicting the NA-removal convention set in the rdepth/taglength
  pairs. On clean data the filter is numerically exact on both datatypes
  with the boundary retained (>= threshold, matching the verbose summary;
  the @return wording says "greater than").

## Findings

### F1 — gl.drop.loc delegation leaks a history entry (MEDIUM, confidence: certain)

Rule: FS8 (one call, one history entry); campaign precedent:
gl.filter.monomorphs (PR #251), gl.filter.allna (PR #252). Location:
`R/gl.filter.reproducibility.r:142,248-250`.

`x2 <- gl.drop.loc(x, loc.list, verbose = 0)` appends gl.drop.loc's own
`match.call()` to the history; this function then appends its own — one
user call yields two new entries, the first exposing an internal
implementation detail (`gl.drop.loc(x = x, loc.list = loc.list, verbose =
0)`). Verified: history length 1 → 3.

Proposed change: restore the pre-call history after the delegation
(`x2@other$history <- x@other$history`) before appending this function's
own entry — the established fix from PR #251.

### F2 — Missing metric is a silent no-op; irrelevant AlleleID/CloneID demanded (MEDIUM, confidence: certain)

Rule: spec axis (function silently fails to do its job). Location:
`R/gl.filter.reproducibility.r:92-100,108-128`.

If RepAvg (SNP) or Reproducibility (SilicoDArT) is absent,
`which(NULL < threshold)` is empty and the dataset is returned unfiltered
— no error, no warning. Verified: 255 loci in, 255 out, silence. The
user believes their data has been quality-filtered when nothing happened.
Meanwhile lines 92-100 stop on missing AlleleID/CloneID — metrics this
function never uses (check copied from the secondaries sibling), so a
dataset with RepAvg but no AlleleID errors needlessly.

Proposed change: check the actual metric per datatype and stop with the
same message as the report sibling; drop the AlleleID/CloneID check.
Escalation-gate class: missing-metric datasets change from silent
pass-through to fatal error; AlleleID-less datasets change from error to
working.

### F3 — plot.file + plot.display=FALSE crashes (MEDIUM, confidence: certain)

Rule: spec axis; PLT3. Location:
`R/gl.filter.reproducibility.r:155-204,208-213`.

Same class as PRs #255-#257: plots built only inside `if (plot.display)`,
save block references `p3` unconditionally. Verified: "object 'p3' not
found". Proposed change: build plots unconditionally; gate only
`print(p3)`.

### F4 — NA-metric loci silently retained (LOW, confidence: certain)

Rule: spec axis (DAT NA lens); consistency with PRs #256/#258. Location:
`R/gl.filter.reproducibility.r:116,127`.

`which(repeatability < threshold)` drops NA comparisons, so NA-RepAvg
loci pass the filter silently (verified: 10/10 retained; no desync — the
delegation keeps sync). The rdepth/taglength filters now remove
NA-metric loci with an itemised count. Proposed change: include NA-metric
loci in the removal list and itemise them in the verbose >= 3 summary,
aligning the convention. Escalation-gate class: NA-bearing datasets lose
their NA-metric loci.

### F5 — Threshold warning ungated (LOW, confidence: certain)

Rule: VRB2. Location: `R/gl.filter.reproducibility.r:82-90`. Verified: 2
lines at verbose = 0. Proposed change: gate at `verbose >= 1`.

### F6 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent). Location:
`R/gl.filter.reproducibility.r:258`. Proposed change:
`return(invisible(x2))`.

### F7 — @family leaks into help title; header conformance; dead code (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7, STY3. Location:
`R/gl.filter.reproducibility.r:1-50,110-131,226-246`.

- Line 4 indented `@family` → rendered title ends "…at a locus @family
  matched filter", no \concept (verified — same as PRs #256/#258).
- DOC2 wording; DOC7 Custodian only; `@return` after `@export`;
  `@return` says "greater than the specified threshold" where behaviour
  (and the verbose summary) is >=.
- Dead code: two commented-out loops, 21-line commented SAVE
  INTERMEDIATES block, redundant `loc.list <- array(NA, nLoc(x))`
  initialisation and post-hoc NA strip.

Proposed change: rewrite header to the ratified template; delete dead
code.

### F8 — Structure: section order (LOW, confidence: certain)

Rule: FS4. Location: `R/gl.filter.reproducibility.r:60-75`. SET WORKING
DIRECTORY and SET COLOURS precede FLAG SCRIPT START. Proposed change:
move FLAG SCRIPT START to directly follow SET VERBOSITY.

## Report notes (other functions / not fixed here)

- gl.drop.loc's own history append is correct behaviour for direct calls;
  the leak is this caller's to fix (per PR #251 precedent).
- No flag reset needed (per-locus metrics remain valid for retained
  loci).

## Coverage

Characterization baseline:
`tests/testthat/test-gl.filter.reproducibility.R` — 14 assertions:
retained counts vs independent recomputation on both datatypes, boundary
retained, metadata sync, input untouched, history +2 (baseline), silent
no-op on missing metric (baseline), NA loci retained (baseline), p3 crash
(baseline), ungated warning (baseline), visible return (baseline). All 14
pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all eight findings
**approved**, including the three escalation-gate consequences (F2:
missing metric silent no-op → fatal error, AlleleID requirement dropped;
F4: NA-metric loci removed; F1/F3 behaviour fixes).

## Outcome

Applied F1–F8. Verification: all 16 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (history
+1 with only this function's own entry [F1]; missing metric errors like
the report sibling [F2]; RDS saved without display [F3]; NA loci removed
and itemised — 22 discarded of which 10 NA on the injected set [F4];
warnings silent at verbose = 0 [F5]; invisible return [F6]). End-to-end
at verbose = 3 on the NA-injected SNP set and testset.gs (240 retained =
independent recomputation). The @family concept now renders; sibling
matched-filter Rd cross-references update mechanically (F7).

Caller grep across dartR.base + 7 siblings: two code callers in
dartR.captive (gl.filter.parent.offspring, gl.report.parent.offspring) —
both pre-check RepAvg presence themselves (the new F2 error cannot fire
for them) and assign the return; gl.report.parent.offspring has the same
verbose = 0 + plot.display = plot.filters interaction as recorded for PR
#256, covered by the same standing decision (guard kept; noted for the
dartR.captive review). gl2fasta example calls with partial-matched `t=1`
— docs only. dartr2shiny: not present in the workspace. NEWS entry
added.

```json
{
  "function": "gl.filter.reproducibility",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["FS8"], "loc": "R/gl.filter.reproducibility.r:142,248-250", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.reproducibility.r:92-100,108-128", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["spec", "PLT3"], "loc": "R/gl.filter.reproducibility.r:155-204,208-213", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.filter.reproducibility.r:116,127", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.filter.reproducibility.r:82-90", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.reproducibility.r:258", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7", "STY3"], "loc": "R/gl.filter.reproducibility.r:1-50,110-131,226-246", "status": "proposed"},
    {"id": "F8", "severity": "LOW", "rules": ["FS4"], "loc": "R/gl.filter.reproducibility.r:60-75", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.filter.reproducibility.R",
  "pr": 260
}
```
