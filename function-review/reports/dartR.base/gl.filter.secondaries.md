# Review: gl.filter.secondaries (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.secondaries
- Datasets: platypus.gl (working dataset — 9 multi-SNP clones), testset.gl
  (no-secondaries case), testset.gs (SilicoDArT rejection)
- Family mode: modify
- Checks skipped: FBM (dartR_fbm) code path not exercised (not available:
  no FBM-converted fixture in this run); Google Group not searched (not
  available: no browser session — GitHub issues in
  green-striped-gecko/dartR.base showed no open complaint naming this
  function).

## Verdicts

- **Standards: FAIL** — invalid-method warning prints at verbose = 0 (VRB2),
  visible return (family precedent is invisible), unreachable SilicoDArT
  block, unused imports, header breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — `method = "best"` does not select the best SNP: the sort
  runs on the full AlleleID string, which is unique per locus, so the
  RepAvg/AvgPIC keys never engage and selection is effectively alphabetical.
  Verified: 5 of platypus.gl's 9 multi-SNP clones keep the wrong locus,
  including two unambiguous RepAvg losses (0.95 kept over 1.0; 0.9898 kept
  over 1.0). Additionally the output locus order is scrambled for both
  methods, even when nothing is removed.

## Findings

### F1 — method='best' selects alphabetically, not by RepAvg/AvgPIC (HIGH, confidence: certain)

Rule: spec axis (documented selection criterion not implemented).
Location: `R/gl.filter.secondaries.r:95-101`.

`order(AlleleID, -RepAvg, -AvgPIC)` sorts on AlleleID first. AlleleID
includes the SNP descriptor (e.g. `45066085|F|0-23:T>G-23:T>G`) and is
unique per locus, so ties never occur and the `-RepAvg, -AvgPIC` keys are
dead code. Within a clone, the retained SNP is whichever AlleleID string
sorts first — its tag position, nothing to do with quality. Verified on
platypus.gl: 5/9 multi-SNP clones keep a locus that is not the
best-RepAvg/AvgPIC one; clone 45066085 keeps RepAvg 0.95 over 1.0, clone
45067382 keeps 0.9898 over 1.0.

Proposed change: rank within clone (the extracted CloneID `b`), not within
AlleleID: `ord <- order(b, -RepAvg, -AvgPIC)`, keep the first occurrence of
each clone in that ranking, then subset by the kept indices. Behaviour
change: the set of retained loci changes for method='best' on any dataset
with secondaries — escalation-gate class (numerical output changes
downstream).

### F2 — Output locus order is scrambled, even for a no-op filter (MEDIUM, confidence: certain)

Rule: spec axis / DAT (least-surprise for a filter). Location:
`R/gl.filter.secondaries.r:100-110`.

Both branches reorder the whole object (by AlleleID for 'best', by a random
permutation for 'random') and never restore the input order. Verified:
testset.gl (zero secondaries) returns all 255 loci in random order — a
filter that removes nothing still shuffles the loci, breaking positional
comparisons with the input and making default runs non-reproducible in
locus order. Also the shuffle idiom `sample(1:(n + 10000), size = n)` is an
obfuscated `sample(n)`.

Proposed change: compute the kept locus indices per clone (per F1 for
'best', random pick for 'random'), then subset once with `sort(keep)` so
the output preserves the input locus order. Escalation-gate class (output
locus order changes for both methods; a no-op filter becomes a true no-op
apart from the history entry).

### F3 — Invalid-method warning ungated and method never reassigned (LOW, confidence: certain)

Rule: VRB2. Location: `R/gl.filter.secondaries.r:54-56`.

`cat(warn("  Warning: method must be specified, set to 'random'\n"))`
prints at verbose = 0. The message claims method is "set to 'random'" but
the variable is never reassigned — behaviour is random only because the
else-branch catches everything non-'best'. Verified: 1 line at verbose = 0.

Proposed change: gate at `verbose >= 1` and actually assign
`method <- "random"`.

### F4 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent: gl.filter.allna, gl.filter.monomorphs
both changed to invisible in PRs #251/#252). Location:
`R/gl.filter.secondaries.r:145`.

`return(x3)` auto-prints the genlight summary on a bare call. Proposed
change: `return(invisible(x3))`.

### F5 — Unreachable SilicoDArT block (LOW, confidence: certain)

Rule: STY3. Location: `R/gl.filter.secondaries.r:58-67`.

`utils.check.datatype(x, accept = c("genlight", "SNP"))` rejects
SilicoDArT before the `if (datatype == "SilicoDArT")` Reproducibility check
can run (verified: testset.gs errors in the datatype check). The docs
correctly say the filter is not implemented for presence/absence data.
Proposed change: delete the dead block.

### F6 — Header: tag order, verbose wording, author, unused imports, typo (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.filter.secondaries.r:1-35`.

- Tag order: `@return` sits after `@export`; ratified order is description,
  details, param, return, author, examples, export.
- `@param verbose` uses the old wording — replace with the ratified canon.
- `@author` has Custodian only — add "Author(s): Arthur Georges." (DOC7).
- `@importFrom stats dpois` and `@import patchwork` are unused (no Poisson,
  no plots in this function — copied from the report sibling). Remove.
- Description typo: "on based on repeatability".
- Example uses platypus.gl without `require("dartR.data")`.
- Add a note that method='random' depends on the RNG state (set.seed for
  reproducibility).

Proposed change: rewrite the header to the ratified template (docs-only).

## Report notes (other functions / not fixed here)

- gl.report.secondaries (PR #253) is the matched report; its counts agree
  with this filter's removal set (9 secondaries on platypus.gl).
- History is appended even when no loci are removed (testset.gl) —
  consistent with gl.filter.monomorphs; left as-is.
- The `verbose > 2` / `> 1` / `> 0` comparison style (vs the house `>=`)
  is equivalent in effect; left as-is.

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.secondaries.R` —
20 assertions: one-locus-per-clone + all clones represented + ploidy/nInd
(random, seeded), genotype-metadata sync, current 'best' selection for the
two unambiguous clones (baseline), scrambled no-op order on testset.gl
(baseline), verbose-0 silence for valid method + 1-line warning for invalid
(baseline), history +1, input untouched, visible return (baseline),
SilicoDArT rejected. All 20 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all six findings
**approved**, including both escalation-gate consequences:

- F1: retained-locus set changes for method='best' (true best-RepAvg/AvgPIC
  SNP kept; downstream results differ) — approved.
- F2: output locus order preserved for both methods (no-op filter returns
  loci unchanged) — approved.
- F3, F4, F5, F6 — approved.

## Outcome

Applied F1–F6. Verification: all 30 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding ('best' now
keeps the max-RepAvg SNP for every multi-SNP clone with AvgPIC tie-break
[F1]; input locus order preserved, testset.gl no-op returns identical
locNames [F2]; invalid-method warning silent at verbose = 0 and method
coerced [F3]; invisible return [F4]). End-to-end at verbose = 3 for both
methods and the no-op case; filter retains 991 loci on platypus.gl,
agreeing with gl.report.secondaries' n.total.tags = 991. Caller grep
across dartR.base + 7 siblings: doc/example references only
(gl.report.secondaries example, dartRstartup tutorial); no code caller
depends on the selection or ordering. dartr2shiny: not present in the
workspace. NEWS entry added. Note: with method='random' the retained set
for a given seed differs from pre-fix runs (the shuffle idiom changed from
`sample(1:(n+10000), n)` to `sample(n)`) — random selection was never
seed-stable across versions and falls under the approved F2 rework.

```json
{
  "function": "gl.filter.secondaries",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.filter.secondaries.r:95-101", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.filter.secondaries.r:100-110", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.filter.secondaries.r:54-56", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.secondaries.r:145", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["STY3"], "loc": "R/gl.filter.secondaries.r:58-67", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.secondaries.r:1-35", "status": "proposed"}
  ],
  "datasets": ["platypus.gl", "testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.filter.secondaries.R",
  "pr": null
}
```
