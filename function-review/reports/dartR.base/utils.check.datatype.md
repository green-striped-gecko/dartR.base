# Review: utils.check.datatype (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: ddaed27 (dev, synced with upstream/dev);
  Branch: review-utils.check.datatype; Datasets: testset.gl,
  testset.gs, constructed abnormal-ploidy and all-NA-individual
  objects
- Family mode: analysis/utility (the package's central datatype
  dispatcher — called on entry by nearly every exported function;
  blast radius assessed via a family-wide census of accept= usage)
- Checks skipped: FBM-backed objects not exercised (no FBM session);
  Google Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — the all-NA screen runs a FULL
  gl.filter.allna pass on every call at verbose >= 2 (measured ~53 ms
  per call on the small testset — a tax on every function entry in
  the package); the unknown-class warning prints at verbose 0
  (verified); dead code (the unreachable misspecification stop — the
  ploidy-null case is fataled earlier, so the `is.null(ploidy)==F`
  branch always captures the remainder; the unused `type` variable
  and commented remnants in the fd branch); `== F` idiom.
- **Spec: FAIL** — (1) the all-NA INDIVIDUALS warning carries the
  LOCI wording ("loci that are scored NA across all individuals" for
  an object whose defect is an all-NA individual — verified). (2)
  Mixed or abnormal ploidy (e.g. c(1,2,...) or uniform 4) is silently
  classified "SNP" with no notice at any verbosity (verified) — note
  this catch-all is load-bearing for the family's polyploid paths, so
  it must warn, not fail. (3) accept = "genlight" (or "dartR") alone
  is fatal for every genlight object, because datatype is never
  literally "genlight" (verified) — latent: the family census shows
  all 10 such callers pair it with SNP/SilicoDArT. (4) @return
  documents 'mat' and 'glPCA' where the code returns "matrix" and
  "glPca"; the details reference gl.check/fd.check/dist.check/
  mat.check parameters from a defunct interface.

## Findings

### F1 — All-NA screen: full filter per call + wrong message text (MEDIUM-HIGH, confidence: certain)

Rule: spec/perf. Location: all-NA block.

Proposed change: replace the gl.filter.allna call with a direct
single-pass check (colSums/rowSums over is.na of the matrix) —
strictly less work, same warnings; and correct the individuals
message to "individuals that are scored NA across all loci".
Consequence: the individuals-case message text changes (right for the
first time); the loci-case text and all verbosity behaviour are
unchanged.

### F2 — Abnormal ploidy silently SNP (MEDIUM, confidence: certain)

Rule: DAT dispatch. Location: ploidy branch.

Proposed change: retain the SNP classification (the polyploid support
paths depend on it) but emit a gated (verbose >= 2) warning that
ploidy is not uniformly 2 and the data are being treated as SNP;
remove the unreachable stop.

### F3 — accept = "genlight"/"dartR" semantics (LOW-MEDIUM, confidence: certain)

Rule: spec/API. Location: accept check.

Proposed change: an accept vector containing "genlight" or "dartR"
now matches both genotype datatypes (SNP and SilicoDArT) — the
intuitive semantics. Backward compatible: every existing caller
already lists the specific datatypes alongside.

### F4 — Unknown-class warning ungated (LOW, confidence: certain)

Rule: VRB. Gated at verbose >= 1 (the accept-check fatal that follows
names the class anyway).

### F5 — Header and style (LOW, docs-only, confidence: certain)

Rules: DOC, STY. @return aligned to the actual strings ("matrix",
"glPca"); the defunct gl.check/fd.check note removed; accept default
documented correctly (bracket typo, dartR included); dead `type` and
commented remnants removed; `== F` idiom.

## Report notes (other functions / not fixed here)

- The per-call all-NA screen, even in its cheaper form, remains an
  O(individuals x loci) scan on every function entry at verbose >= 2.
  A flags-based approach (checking loc.metrics.flags$allna) would
  eliminate it entirely but changes semantics when flags are stale —
  a design question for the coordinators, recorded not actioned.

## Coverage

`tests/testthat/test-utils.check.datatype.R` — 16 assertions:
dispatch across six classes, invisibility, accept rejection,
abnormal-ploidy silence (baseline), wrong individuals wording
(baseline), loci warning presence and verbose-0 silence,
unknown-class leak (baseline), accept-genlight fatal (baseline). All
16 pass on the pre-fix code.

## Approval

All findings (F1-F5) approved via the approval boxes (2026-08-31).

## Outcome

All findings applied. F3 implemented with a governs-clause: 'genlight'
or 'dartR' in accept admits both genotype datatypes ONLY when no
specific genotype datatype is listed — so the family's
c("genlight","SNP") callers retain their SNP-only restriction
(verified by assertion). Own suite: 19/19; flips map to F1 (correct
individuals wording), F2 (ploidy notice at >= 2, silent at 0), F3
(accept="genlight" returns the datatype; specific listing governs),
F4 (no unknown-class leak at verbose 0). Blast radius verified by a
full-tree regression sweep of all 53 committed test files: the same
nine failures occur with and without this change — all nine are
PRE-EXISTING on upstream dev (ddaed27) and are recorded below. NEWS
entry added.

## Pre-existing dev failures (recorded for the merge session, not
caused or touched by this review)

test-gl.filter.hwe.R:29, test-gl.fixed.diff.R:37/40,
test-gl.report.allelerich.R:21-24, test-gl.report.basics.R:22,
test-gl.report.hwe.R:83 — failing on clean upstream dev. Symptoms
suggest merge-resolution interactions among the recently merged
review PRs (e.g. gl.report.hwe returns the old 245-row
verbosity-dependent tested set where the merged test expects the
fixed 196) and possibly the dartR.data 1.2.4 dependency (PR #10 over
there). Flagged for Bernd/the coordinators when the PR backlog is
merged.

## Addendum 2026-09-06 (v2 review)

The second-pass review (`utils.check.datatype-v2.md`, this directory on the
follow-up branch) established that on the reviewed state (ddaed27,
upstream/dev) the courtesy all-NA screen could kill the caller's run
outright: at `verbose >= 2` (the default), `gl.filter.allna` cannot return
a zero-locus object, so any object in which every locus is scored NA died
in the preamble with "Subsetting resulted in zero loci.", and a zero-locus
object died with "subscript out of bounds" — before the calling function,
or its own `accept` gate, ever saw the object. The same call at
`verbose = 0/1` proceeded: a verbosity flag changed control flow. This is
the residual lead recorded by the PR #367 (gl.compliance.check) review,
now confirmed and root-caused (v2 review, its finding F1, HIGH).

This PR's F1 change (the direct `colSums`/`rowSums` scan) removes the
filter call and with it every crash path. The fix was not recorded here at
the time — the v1 review replaced the call for cost and message-wording
reasons without observing the crash it also removed. The behaviour is now
pinned against regression by the edge baseline added with this addendum
(`tests/testthat/test-utils.check.datatype-edges.R`, 37 assertions,
passing on this branch), which covers the all-all-NA and zero-locus
survivals alongside the degenerate-input, accept-gate, FBM and read-only
contracts.

```json
{"function": "utils.check.datatype", "package": "dartR.base", "family_mode": "utility",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "MEDIUM-HIGH", "rules": ["spec", "perf"], "loc": "R/utils.check.datatype.r allna", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["DAT"], "loc": "R/utils.check.datatype.r ploidy", "status": "applied"},
  {"id": "F3", "severity": "LOW-MEDIUM", "rules": ["spec", "API"], "loc": "R/utils.check.datatype.r accept", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["VRB"], "loc": "R/utils.check.datatype.r unknown class", "status": "applied"},
  {"id": "F5", "severity": "LOW", "rules": ["DOC", "STY"], "loc": "R/utils.check.datatype.r header", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs", "constructed"],
 "baseline_test": "tests/testthat/test-utils.check.datatype.R", "pr": 296}
```
