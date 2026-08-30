# Review: gl.recode.ind (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.recode.ind
- Datasets: testset.gl with dartR.data's testset_ind_recode.csv (250
  rows, one Delete: UC_00126c) and a crafted pure-rename variant
- Family mode: modify (data manipulation)
- Checks skipped: SilicoDArT path not exercised end-to-end (not
  available: no ind recode fixture for testset.gs); FBM path not
  exercised (not available: no FBM fixture); Google Group not searched
  (not available: no browser session — GitHub issues showed no open
  complaint naming this function).

## Verdicts

- **Standards: FAIL** — near-clone of gl.recode.pop (PR #266) with the
  same defect family: flag state depends on verbosity, the Delete path
  appends two history entries, visible return, DOC1/DOC2/DOC7 drift; plus
  a warning leak at verbose = 0 and a results summary gated at >= 2 where
  the pop sibling uses >= 3.
- **Spec: FAIL** — the verbose = 3 deletions listing prints the literal
  word "Delete" instead of the deleted individual's original identifier
  (it lists post-rename names, and deleted individuals have already been
  renamed to Delete); the `verbose == 3` gate hides the listing at
  verbose = 5; and even verbose = 0 emits a line (gl.drop.ind's "Listed
  individual delete not present — ignored" warning, because the
  delegation always passes both spellings). The recode mapping and the
  deletion itself are exact (verified; input untouched; error paths
  work); `@return` claims "A genlight or genind object" — genind is
  never returned.

## Findings

### F1 — Deletions listing useless; verbose-0 warning leak; wrong gate (MEDIUM, confidence: certain)

Rule: spec axis. Location: `R/gl.recode.ind.r:139-160`.

- The listing indexes post-rename `indNames(x)` by recode-table row
  order; the deleted individual has been renamed to "Delete", so the
  listing prints " Delete " instead of UC_00126c (verified). The
  informative content is the ORIGINAL identifiers, available directly
  from the recode table.
- `verbose == 3` hides the listing at verbose = 5 (verified).
- `gl.drop.ind(x, ind.list = c("Delete", "delete"), verbose = 0)` warns
  "Listed individual delete not present — ignored" even at verbose = 0
  (verified: 1 line leaks) because only one spelling is ever present.

Proposed change: list
`recode.table[tolower(recode.table[, 2]) == "delete", 1]` (original
names) at `verbose >= 3`; pass only the spellings actually present:
`intersect(c("Delete", "delete"), indNames(x))`.

### F2 — Flag state depends on verbosity (MEDIUM, confidence: certain)

Rule: VRB (object content must not depend on verbosity). Location:
`R/gl.recode.ind.r:187-192`.

Identical to gl.recode.pop F2 (PR #266): `utils.reset.flags` inside
`if (verbose >= 2)`. Pure renaming invalidates flags at verbose >= 2 but
not below (verified: CallRate FALSE vs TRUE); the Delete path is already
reset by the gl.drop.ind delegation at every verbosity (verified).
Proposed change: remove the misplaced call; keep the gated warning.

### F3 — History delegation leak on the Delete path (MEDIUM, confidence: certain)

Rule: FS8; campaign precedent PRs #251/#260/#266. Location:
`R/gl.recode.ind.r:156-159,213-215`.

Verified: +2 entries on a Delete run (gl.drop.ind's own call plus this
function's). Proposed change: snapshot the history before the
delegations and restore it before the single own append.

### F4 — Results summary gated at verbose >= 2 (LOW, confidence: certain)

Rule: VRB3 (results summary belongs at >= 3; the pop sibling uses >= 3).
Location: `R/gl.recode.ind.r:197`. Proposed change: gate at
`verbose >= 3`.

### F5 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5. Location: `R/gl.recode.ind.r:224`. Proposed change:
`return(invisible(x))`.

### F6 — Unguarded monomorphs-flag test (LOW, confidence: certain)

Rule: DEP guard. Location: `R/gl.recode.ind.r:171`. Same as
gl.recode.pop F5. Proposed change: `isFALSE(...)`.

### F7 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC5, DOC7. Location: `R/gl.recode.ind.r:1-57`.

- `@return` "A genlight or genind object" — genind is never returned;
  corrected.
- `@details` after `@param`; `@return` after `@export` (DOC1); DOC2
  wording variant; `@author` Custodian only (DOC7).

Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- gl.drop.ind's "not present — ignored" warning prints at verbose = 0 —
  a defect in that (already-reviewed) function's warning gate; this
  review sidesteps it by passing only present spellings, but the callee
  gate deserves a follow-up when gl.drop.ind is next touched.
- The O(nInd × ntr) recode loop is correct; left as-is (as in PR #266).

## Coverage

Characterization baseline: `tests/testthat/test-gl.recode.ind.R` — 15
assertions: mapping vs independent recomputation, fixture deletion
(UC_00126c), input untouched, short-table error, verbose-0 warning leak
(baseline), history +2 (baseline), "Delete" listing + absent-at-5
(baseline), verbosity-dependent flags (baseline), summary at verbose 2
(baseline), visible return (baseline). All 15 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all seven findings
**approved**, including the escalation-gate consequences (F1: corrected
listing + silent verbose 0; F2: flags no longer invalidated by pure
renaming; F4: summary moves from verbose >= 2 to >= 3).

## Outcome

Applied F1–F7. Verification: all 17 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (listing
names UC_00126c and appears at verbose >= 3 including 5; verbose = 0
fully silent — only present spellings passed to gl.drop.ind [F1];
pure-rename flags valid at every verbosity, deletion runs still reset
via the delegation [F2]; exactly one history entry bearing this
function's own call [F3]; summary at >= 3 [F4]; invisible return [F5]).
Caller grep across dartR.base + 7 siblings: one doc mention
(gl.make.recode.ind) — no code caller. dartr2shiny: not present in the
workspace. NEWS entry added. The gl.drop.ind verbose-0 warning gate
remains flagged as a callee follow-up (report note).

```json
{
  "function": "gl.recode.ind",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.recode.ind.r:139-160", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["VRB"], "loc": "R/gl.recode.ind.r:187-192", "status": "proposed"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["FS8"], "loc": "R/gl.recode.ind.r:156-159,213-215", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB3"], "loc": "R/gl.recode.ind.r:197", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.recode.ind.r:224", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DEP"], "loc": "R/gl.recode.ind.r:171", "status": "proposed"},
    {"id": "F7", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC5", "DOC7"], "loc": "R/gl.recode.ind.r:1-57", "status": "proposed"}
  ],
  "datasets": ["testset.gl + testset_ind_recode.csv", "crafted pure-rename variant"],
  "baseline_test": "tests/testthat/test-gl.recode.ind.R",
  "pr": null
}
```
