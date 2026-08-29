# Review: gl.filter.taglength (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.taglength
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT — has TrimmedSequence)
- Family mode: modify
- Checks skipped: FBM path not exercised (not available: no FBM fixture);
  Google Group not searched (not available: no browser session — GitHub
  issues showed no open complaint naming this function).

## Verdicts

- **Standards: FAIL** — threshold warnings print at verbose = 0 (VRB2), one
  of them names the wrong parameter ('verbose' instead of 'lower'), visible
  return, indented `@family` tag leaking into the rendered help title,
  header breaks DOC1/DOC2/DOC7.
- **Spec: FAIL** — loci with NA tag length silently desync genotypes from
  the locus metadata (185 loci vs 195 metric rows with 10 NAs injected) —
  the same corruption class as gl.filter.rdepth (PR #256) — and the
  progress message says loci *between* the thresholds are removed when
  those are the loci retained. On clean data the filter is numerically
  correct on both datatypes with inclusive boundaries, and the
  threshold-swap safeguard works.

## Findings

### F1 — NA tag lengths silently desync genotypes and metadata (HIGH, confidence: certain)

Rule: DAT2/DAT3. Location: `R/gl.filter.taglength.r:117-121`.

`index <- (nchar.tags >= lower & nchar.tags <= upper)` is NA wherever
TrimmedSequence is NA; the genlight subset drops those loci but the
data.frame subset keeps an all-NA row for each. Verified: 10 NAs → 185
loci vs 195 metric rows. Identical class to gl.filter.rdepth F1 (PR #256).

Proposed change:
`index <- !is.na(nchar.tags) & nchar.tags >= lower & nchar.tags <= upper`,
with NA-length removals itemised in the verbose >= 3 summary.
Escalation-gate class: on affected datasets the output changes from a
corrupt object to a clean one with NA-length loci removed.

### F2 — Progress message states the opposite of the behaviour (LOW, confidence: certain)

Rule: spec axis. Location: `R/gl.filter.taglength.r:108-116`.

"Removing loci with taglength between 60 and 69" — those are the loci
retained. Proposed change: "Removing loci with tag length < lower or >
upper".

### F3 — Threshold warnings ungated and mislabelled (LOW, confidence: certain)

Rule: VRB2. Location: `R/gl.filter.taglength.r:67-95`.

The swap warning and both range warnings print at verbose = 0 (verified:
2-4 lines). The lower-range warning text reads "Parameter 'verbose' must
be an integer between 0 and 250, set to 20" — it is checking `lower`, not
verbose (copy-paste slip). Proposed change: gate all three at
`verbose >= 1`; correct the parameter name in the message.

### F4 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent PRs #251/#252/#254/#256). Location:
`R/gl.filter.taglength.r:154`. Proposed change: `return(invisible(x2))`.

### F5 — @family tag leaks into the help title; header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.filter.taglength.r:1-37`.

- Line 4 `#'  @family matched filter` is indented after the title
  continuation — rendered `?gl.filter.taglength` title ends with "…tag
  length @family matched filter" and the Rd has no \concept (verified,
  same defect as gl.filter.rdepth, PR #256).
- DOC2 verbose wording; DOC7 author missing "Author(s):"; `@return` after
  `@export`; no `@seealso` to the matched report; description is one
  sentence with no statement of what the function does; no @details.
- Stray code comment "# Remove SNP loci with rdepth < threshold" (line
  107) copied from the rdepth sibling — corrected in passing.

Proposed change: rewrite the header to the ratified template (docs-only;
fixes the visible title corruption); fix the stray comment.

## Report notes (other functions / not fixed here)

- No flag reset needed: removing loci leaves per-locus metrics valid
  (consistent with gl.filter.secondaries/rdepth).
- The threshold-swap safeguard (upper < lower → swap) works correctly and
  is retained as-is apart from the verbosity gate.

## Coverage

Characterization baseline: `tests/testthat/test-gl.filter.taglength.R` —
15 assertions: retained counts vs independent recomputation on both
datatypes, inclusive boundaries, metadata sync on clean data, swap
safeguard, history +1, input untouched, TrimmedSequence-absent error, NA
desync 185/195 (baseline), ungated warning (baseline), visible return
(baseline), wrong message wording (baseline). All 15 pass on the pre-fix
code.

## Approval

Decisions recorded 2026-08-29 via approval boxes — all five findings
**approved**, including the escalation-gate consequence (F1: NA-length
loci removed cleanly instead of corrupting the object).

## Outcome

Applied F1–F5. Verification: all 17 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (NA case:
185 loci = 185 metric rows, no NA rows, itemised in the verbose >= 3
summary [F1]; message "tag length < 60 or > 69" [F2]; warnings silent at
verbose = 0 with corrected parameter name [F3]; invisible return [F4]).
End-to-end at verbose = 3 on the NA-injected SNP set and testset.gs
(196 = independent recomputation). The @family concept now renders and
the sibling matched-filter Rd files gain the cross-reference line
(mechanical roxygen consequence of F5). Caller grep across dartR.base +
7 siblings: one code caller — gl.dist.phylo assigns the return at
verbose = 0, unaffected by all changes; remainder are docs/tutorial
references. dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.filter.taglength",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["DAT2", "DAT3"], "loc": "R/gl.filter.taglength.r:117-121", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.filter.taglength.r:108-116", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.filter.taglength.r:67-95", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.taglength.r:154", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.taglength.r:1-37", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.filter.taglength.R",
  "pr": null
}
```
