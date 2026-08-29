# Review: gl.filter.overshoot (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.overshoot
- Datasets: testset.gl (SNP — 21 genuine overshoot loci, live fixture),
  crafted no-op variant (SnpPosition all 0), NA-injected variant,
  testset.gs (SilicoDArT rejection)
- Family mode: modify
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — the no-op branch's message prints at verbose = 0
  (VRB2), visible return in both branches, duplicated end-of-script code,
  header breaks DOC1/DOC2.
- **Spec: PASS (core logic)** — removal is numerically exact against
  independent recomputation (21 removed, 234 retained on testset.gl, no
  overshoot loci left, metadata in sync, one history entry when
  filtering). Loci whose overshoot status cannot be assessed (NA
  TrimmedSequence) silently pass the filter, contradicting the NA-removal
  convention set in PRs #256/#258/#260 — the only behavioural gap.

## Findings

### F1 — No-op message ungated (LOW, confidence: certain)

Rule: VRB2. Location: `R/gl.filter.overshoot.r:108`.

"There were no loci…" prints at verbose = 0 (verified: 1 line). The
filter branch is correctly gated. Proposed change: gate at
`verbose >= 1`.

### F2 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent PRs #251-#260). Location:
`R/gl.filter.overshoot.r:105,115`.

Both branches return visibly; the end-of-script block is duplicated in
each branch. Proposed change: single exit with `return(invisible(...))`.

### F3 — NA-assessable loci silently retained (LOW, confidence: certain)

Rule: spec axis (DAT NA lens); consistency with PRs #256/#258/#260.
Location: `R/gl.filter.overshoot.r:80`.

`which(snpos > nchar(trimmed))` drops NA comparisons, so loci with NA
TrimmedSequence (or NA SnpPosition) pass the filter silently (verified:
5/5 retained; no desync). The convention set in the rdepth/taglength/
reproducibility filters is to remove loci that cannot be assessed, with
an itemised count. Proposed change: include NA-assessable loci in the
removal set and itemise them at verbose >= 3. Escalation-gate class:
NA-bearing datasets lose those loci.

### F4 — Listing trailing comma (LOW, confidence: certain)

Rule: STY. Location: `R/gl.filter.overshoot.r:87`.

`cat(paste0(locNames(x)[os], ","))` glues a comma to each name, leaving a
stray trailing comma (verified at verbose = 3) — same slip fixed in the
report sibling (PR #261). Proposed change:
`paste(locNames(x)[os], collapse = ", ")`.

### F5 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2. Location: `R/gl.filter.overshoot.r:1-32`.

- Description paragraphs 2-3 (analysis consequences, ordering advice)
  belong under @details; `@return` sits after `@export` (DOC1).
- DOC2 verbose wording stale.
- `@author` has both but with a semicolon ("Author(s): …; Custodian:") —
  normalised to the campaign format.
- No `@seealso` to the matched report.

Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- The no-op branch appends no history entry — consistent with
  gl.filter.allna (PR #252: history only when the object changes); left
  as-is and asserted in the baseline.
- No flag reset needed (per-locus metrics remain valid).

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.overshoot.R` —
14 assertions: removal counts vs independent recomputation (21/234), no
overshoot left, metadata sync, history +1 on filter / +0 on no-op, input
untouched, SilicoDArT + missing-metric errors, ungated no-op message
(baseline), NA loci retained (baseline), visible return (baseline). All
14 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all five findings
**approved**, including the escalation-gate consequence (F3: NA-assessable
loci removed instead of silently passing).

## Outcome

Applied F1–F5. Verification: all 16 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (no-op
message silent at verbose = 0 [F1]; invisible return via a single exit
[F2]; NA-assessable loci removed and itemised — 26 removed of which 5 NA
on the injected set [F3]; clean listing [F4]). End-to-end at verbose = 3.
No-op runs still append no history (asserted; allna precedent). Caller
grep across dartR.base + 7 siblings: two code callers — gl2bpp and
gl2fasta both assign the return at verbose = 0 and are unaffected
(indeed improved: the ungated no-op message no longer leaks into their
output); remainder are docs/tutorial references. dartr2shiny: not present
in the workspace. NEWS entry added.

```json
{
  "function": "gl.filter.overshoot",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.filter.overshoot.r:108", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.overshoot.r:105,115", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["spec", "DAT"], "loc": "R/gl.filter.overshoot.r:80", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["STY"], "loc": "R/gl.filter.overshoot.r:87", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2"], "loc": "R/gl.filter.overshoot.r:1-32", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.filter.overshoot.R",
  "pr": null
}
```
