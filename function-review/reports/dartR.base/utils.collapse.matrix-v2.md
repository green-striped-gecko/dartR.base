# Review: utils.collapse.matrix — v2, chain context (dartR.base)

Follow-up to the infrastructure-wave review (utils.collapse.matrix.md,
PR #325, all findings applied on this branch). This v2 pass was run as
part of the population-distance chain wave because gl.dist.pop's
SilicoDArT methods depend on this collapse; it verifies the
chain-relevant behaviour the first review did not pin, and records how
the divergence between upstream/dev and the PR state affects the chain.
Naming precedent: utils.check.datatype-v2.md.

## Provenance

- Model: Claude Fable 5 (claude-fable-5, Claude Code) via dartr-dev agent;
  Skill: dartr-function-review v2.0.0; Base: upstream/dev at ddaed27
  with `git diff upstream/dev -- R/utils.collapse.matrix.r` NON-empty —
  the local file is the applied PR #325 state (commit 73f6795).
  Handling: the loaded code is the PR state; the upstream state was
  emulated by hand for the diagonal comparison. Off-diagonal values are
  identical between the two states (verified on the fixture), and the
  dist route drops the diagonal entirely, so every chain result through
  gl.dist.pop is identical on ddaed27 and on this branch.
- Datasets: constructed 6-individual/3-population fixture (pops of 3, 2,
  and 1 individuals; non-alphabetical factor levels; one NA pair; one
  all-NA block), testset.gl/testset.gs via the gl.dist.pop chain.
- Family mode: analysis/infrastructure utility.
- Baseline: tests/testthat/test-utils.collapse.matrix.chain.R (7 tests,
  all pass), companion to the tracked infrastructure-wave file
  test-utils.collapse.matrix.R.
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

**Standards: Ready** — at the PR #325 state the stops, gates, and docs
findings are applied; nothing new on this axis.
**Spec: Ready** — the collapse statistic is confirmed and label-safe;
the two new findings are guard-quality items on an internal-only
utility.

The collapse statistic, stated for the record (gl.dist.pop's SilicoDArT
methods inherit it): a between-population cell is the mean over ALL
cross-population individual pairs (n_i x n_j cells), `na.rm = TRUE`; the
within-population diagonal (matrix route only — `as.dist` drops it) is
the mean over distinct pairs (PR #325 state; upstream includes the zero
self-distances, deflating it); a single-individual population collapses
to 0; an all-NA block propagates as NaN. All verified by hand on the
fixture. Population indexing is name-based, so non-alphabetical factor
levels come out correctly labelled (verified — the exact hazard that
bites gl.dist.pop's dcast path, its F1, is absent here).

## Findings

**F1 [LOW, confidence: high] — the name guard is one-directional (FS5)**
`R/utils.collapse.matrix.r:77-79` — the guard checks that every matrix
name exists in the genlight (`rownames(mat) %in% indNames(x)`) but not
that every individual of the genlight exists in the matrix.
Failure scenario: a D computed on a subset of individuals passes the
guard and dies later with "subscript out of bounds" (verified) instead
of a clear message.
Proposed change: also require `indNames(x) %in% rownames(mat)` (or
intersect explicitly and say what was dropped).

**F2 [INFO, confidence: high] — unnamed matrix input fails with a raw dimnames error (FS5)**
An unnamed square matrix passes the square check and dies at the `%in%`
guard with "no 'dimnames' attribute for array" (verified). Internal-only
utility, so INFO: a `is.null(dimnames(mat))` check with a plain message
would do.

## Proposed changes

1. Two-directional name guard with a clear message (F1, F2 — one
   change).

## Coverage

- Chain spec: between-pop statistic vs hand means; PR-vs-upstream
  diagonal emulation (off-diagonals identical; dist route diagonal-free);
  label safety under non-alphabetical levels; NA pair and all-NA block
  propagation; singleton population; dist-in/dist-out and
  matrix-in/matrix-out; missing-individual and unnamed-matrix guards;
  verbose-0 silence — run.
- Standards walk: not repeated — covered by the infrastructure-wave
  review and its applied PR #325.
- Google Group search: SKIPPED — not available, no browser session.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur Georges (2026-09-06) | the two guard notes (F1 LOW, F2 INFO) as one change |

## Outcome

Change 1 applied as an addendum commit on the open PR #325 branch
(`review-utils.collapse.matrix`), per the routing note above: the name
guard is now two-directional and a dimnames-presence check precedes it,
each with a clear fatal message. Verification:

- Chain companion baseline (this file's
  test-utils.collapse.matrix.chain.R, added to the branch in the same
  commit): all tests pass; the missing-individual expectation flipped
  from "subscript out of bounds" to the clear message
  (`[approved F1]`), and a new unnamed-matrix expectation pins the
  dimnames message (`[approved F2]`).
- The infrastructure-wave baseline (test-utils.collapse.matrix.R,
  already on the branch) still passes unchanged -- valid inputs are
  unaffected.
- The gl.dist.pop chain is unaffected: its D always arrives from
  gl.dist.ind carrying the full individual set with names.
- PR #325's body updated with an addendum note (gh pr edit).

```json
{
  "function": "utils.collapse.matrix",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27+73f6795",
  "report_version": 2,
  "verdict_standards": "ready",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "INFO", "confidence": "high", "rule": "FS5", "status": "applied", "change": 1}
  ],
  "coverage_skipped": ["standards walk: covered by v1/PR #325", "Google Group: no browser session"],
  "baseline_test": "tests/testthat/test-utils.collapse.matrix.chain.R",
  "status": "pr-open",
  "pr": 325
}
```
