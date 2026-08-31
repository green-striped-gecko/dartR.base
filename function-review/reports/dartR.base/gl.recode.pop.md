# Review: gl.recode.pop (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.recode.pop
- Datasets: testset.gl with dartR.data's testset_pop_recode.csv (plain
  recode) and a crafted variant flagging EmmacRoss + EmmacTweeUki as
  Delete (20 individuals)
- Family mode: modify (data manipulation)
- Checks skipped: SilicoDArT path not exercised end-to-end (not
  available: no pop recode fixture for testset.gs; the SNP-only gating of
  mono.rm/recalc was verified by code inspection); FBM path not exercised
  (not available: no FBM fixture); Google Group not searched (not
  available: no browser session — GitHub issues showed no open complaint
  naming this function).

## Verdicts

- **Standards: FAIL** — the object's flag state depends on verbosity
  (utils.reset.flags inside the verbose >= 2 gate — the gl.drop.ind bug
  class), the Delete path appends two history entries, visible return,
  header breaks DOC1/DOC2/DOC7 and @seealso links to itself.
- **Spec: FAIL** — the verbose = 3 deletions listing prints the wrong
  individuals and the wrong count: it indexes indNames() by a recycled
  length-30 logical (one entry per recode row), so it listed 16
  arbitrary individuals ("A total of 16 dropped") when 20 were deleted,
  none of the listed names being actual deletions; and the gate is
  `verbose == 3` exactly, so verbose = 5 gets no listing. The recode
  mapping and the Delete mechanics themselves are exact (verified against
  independent recomputation; input untouched; error paths work).

## Findings

### F1 — Deletions listing: wrong individuals, wrong count, wrong gate (MEDIUM, confidence: certain)

Rule: spec axis (recycled-index class). Location:
`R/gl.recode.pop.r:156-165`.

`indNames(x)[tolower(recode.table[, 2]) == "delete"]` indexes 250
individual names by a length-30 logical (one per recode-table row) —
recycling selects arbitrary individuals. Verified: 16 wrong names listed,
"A total of 16 individuals dropped" printed, actual deletions = 20 and
absent from the listing. Also `verbose == 3` should be `>= 3`.

Proposed change: compute the deletions from the recoded pops —
`indNames(x)[as.character(pop(x)) %in% c("Delete", "delete")]` — and gate
at `verbose >= 3`.

### F2 — Flag state depends on verbosity (MEDIUM, confidence: certain)

Rule: VRB (object content must not depend on verbosity; the
gl.drop.ind/keep.ind class). Location: `R/gl.recode.pop.r:197-202`.

`utils.reset.flags(x, verbose = 0)` sits INSIDE `if (verbose >= 2)`. A
pure renaming run (no deletions — locus metrics unaffected) invalidates
every locus-metric flag at verbose >= 2 but leaves them valid at verbose
0-1 (verified: CallRate flag FALSE vs TRUE). In the Delete path the reset
is redundant — the gl.drop.pop delegation already resets flags (verified:
FALSE at every verbosity).

Proposed change: remove the utils.reset.flags call; keep the gated
warning. Pure renaming then preserves valid flags at every verbosity, and
deletion runs keep the reset performed by gl.drop.pop.

### F3 — History delegation leak on the Delete path (MEDIUM, confidence: certain)

Rule: FS8; campaign precedent PRs #251/#260. Location:
`R/gl.recode.pop.r:166-169,223-225`.

gl.drop.pop (and gl.filter.monomorphs / gl.recalc.metrics when invoked)
append their own history entries before this function appends its own —
verified: +2 entries on a Delete run. Proposed change: snapshot the
history before the delegations and restore it before the single own
append.

### F4 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent). Location: `R/gl.recode.pop.r:234`.
Proposed change: `return(invisible(x))`.

### F5 — Unguarded monomorphs-flag test (LOW, confidence: certain)

Rule: DEP guard (the utils.recalc.callrate isTRUE class). Location:
`R/gl.recode.pop.r:181`.

`if (x@other$loc.metrics.flags$monomorphs == FALSE)` errors with
"argument is of length zero" if the flag is absent (crafted objects).
Proposed change: `if (isFALSE(x@other$loc.metrics.flags$monomorphs))`.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. Location: `R/gl.recode.pop.r:1-62`.

- `@details` after `@param`; `@return` after `@export` (DOC1); DOC2
  verbose wording variant.
- `@author` Custodian only (DOC7).
- `@seealso` links to gl.recode.pop itself; corrected to gl.recode.ind
  and gl.drop.pop alongside gl.filter.monomorphs.

Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- The recode loop is O(nInd × ntr) nested exact-match — correct, and
  fast enough at these scales; left as-is.
- The recode-table validation requires the unique population sets to
  match in count exactly (extra unmatched rows in the table are fatal) —
  documented-adjacent behaviour, left as-is.

## Coverage

Characterization baseline: `tests/testthat/test-gl.recode.pop.R` — 16
assertions: mapping vs independent recomputation, Delete mechanics
(counts, pops gone, right individuals gone), input untouched,
short-table error, history +1/+2 (baseline), wrong listing + wrong count
+ absent-at-5 (baseline), verbosity-dependent flags (baseline), visible
return (baseline). All 16 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all six findings
**approved**, including both escalation-gate consequences (F1: corrected
listing; F2: flags no longer invalidated by pure renaming at
verbose >= 2).

## Outcome

Applied F1–F6. Verification: all 19 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (listing
names the 20 true deletions and appears at verbose >= 3 including 5
[F1]; plain-recode flags stay valid at every verbosity while deletion
runs still get gl.drop.pop's reset [F2]; exactly one history entry on
both paths, bearing this function's own call [F3]; invisible return
[F4]). Caller grep across dartR.base + 7 siblings: docs cross-references
and one tutorial call that assigns the return — no incompatibility.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.recode.pop",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.recode.pop.r:156-165", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.recode.pop.r:197-202", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["FS8"], "loc": "R/gl.recode.pop.r:166-169,223-225", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.recode.pop.r:234", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["DEP"], "loc": "R/gl.recode.pop.r:181", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.recode.pop.r:1-62", "status": "proposed"}
  ],
  "datasets": ["testset.gl + testset_pop_recode.csv", "crafted Delete variant"],
  "baseline_test": "tests/testthat/test-gl.recode.pop.R",
  "pr": 266
}
```
