# Review: rbind.dartR (dartR.base)

Follow-up to the `utils.dartR.class.def` review (PR #328), which recorded
the lost metadata as documented behaviour ("rbind is a bit lazy") and did
not change it. Defects reported by the dartR.sim session while reviewing
`gl.sim.offspring`; fix requested by Luis.

- Family mode: modify (combines individuals of several objects)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 1dd6a3a (`origin/dev`)
- Datasets: platypus.gl, testset.gs
- Baseline: tests/testthat/test-rbind.dartR.R (11 expectations), captured
  before review; all pass

## Verdicts

**Standards: Needs work** — the two paths (FBM and in-memory) build their
output differently; argument handling has a dead branch and drops
arguments silently.
**Spec: Rework** — on the in-memory path, which is the default for all
non-FBM data, genotypes are misaligned when locus order differs, and all
`@other` metadata is lost; the FBM path keeps metadata for the first
object's individuals only.

What works well: same-order, same-coding inputs give exact genotypes on
both paths, and the FBM path already aligns loci by name.

## Findings

**F1 [HIGH, confidence: high] — genotypes misaligned (DAT2)**
`R/utils.dartR.class.def.r:868` — the in-memory path concatenates each
object's `@gen` in that object's own locus order, then labels the result
with the first object's locus names. The column map built at :722-740 is
used only by the FBM path.
Failure scenario: `rbind(x[1:5, ], x[6:10, sample(nLoc(x))])` on
platypus.gl gives 1722 of 4636 genotypes (37%) of individuals 6-10 at the
wrong locus, with no error or warning.
Proposed change: reorder every object to the first object's locus order
(`objs[[k]][, colmap[[k]]]`) before concatenating.

**F2 [HIGH, confidence: high] — metadata lost or desynchronised (DAT2, DAT4, FS8)**
In-memory path (:865-881): `@other` comes back empty — `loc.metrics`,
`ind.metrics`, `loc.metrics.flags`, `latlon` and `history` are all lost;
`other_out` (:790) is computed but never used. FBM path (:856): `@other`
is the first object's, so `rbind` of 5 + 5 individuals gives 10
individuals with 5 `ind.metrics` rows. No history entry on either path;
locus-metric flags stay TRUE although the set of individuals changed.
Failure scenario: `gl.sim.offspring` output combined with its parents
(the `gl.sim.relatedness` pattern) loses `ind.metrics`; any per-individual
join after an FBM `rbind` reads the wrong rows.
Proposed change: both paths build `@other` the same way —
`loc.metrics` from the first object (in reference order); `ind.metrics`
row-bound across objects on the union of columns, missing columns filled
with NA (as `gl.join` does); `latlon` row-bound when every object has it;
locus-metric flags reset to FALSE (the individuals changed, so call rate,
allele frequencies and similar metrics are stale); the first object's
history plus one entry for this call; other `@other` elements from the
first object.

**F3 [MEDIUM, confidence: high] — allele coding not checked (DAT1, DAT2)**
Objects whose `loc.all` differs at the same locus (reference and
alternate swapped) are combined silently, and the result takes the first
object's alleles, so genotype 2 means different alleles in different rows.
Failure scenario: two DArT reports of the same loci with opposite allele
orientation are combined; allele frequencies of the merged data are wrong
for every swapped locus.
Proposed change: after aligning loci, stop with an error that names the
number of loci whose `loc.all` differs.

**F4 [MEDIUM, confidence: high] — SNP and SilicoDArT combined (DAT1)**
An SNP object (ploidy 2) and a SilicoDArT object (ploidy 1) with the same
locus names combine silently into one object with ploidy 2 and 1.
Failure scenario: SilicoDArT presence (1) is then read as a heterozygote
by every SNP function.
Proposed change: stop with an error when the inputs mix ploidy.

**F5 [LOW, confidence: high] — argument handling and docs (API2, DOC5)**
`:704-706` — non-genlight arguments are dropped silently
(`rbind(x, "a")` returns `x`); the `is.list(objs[[1]])` branch can never
run, because `objs` holds only genlight objects; `rbind(list(a, b))`
returns a base matrix because S3 dispatch sees a list. The roxygen says
`@param ... list of dartR objects` and that `@other` is lost.
Failure scenario: a user passes a list and gets a matrix; a typo in an
argument is ignored.
Proposed change: remove the dead branch; stop when an argument other than
`NULL` is not a genlight object (`NULL` stays allowed: `gl.impute` starts
from `rbind(NULL, pop)`); document `do.call(rbind, list_of_objects)` and
the new metadata behaviour.

## Proposed changes

1. Align loci by name on the in-memory path (F1). **Consequence: genotypes
   of inputs with different locus order change (they are wrong today).**
2. Carry metadata on both paths: `loc.metrics`, row-bound `ind.metrics`
   and `latlon`, flags reset, history entry (F2). **Consequence: the
   returned object now has `@other` metadata (it was empty on the
   in-memory path); on the FBM path `ind.metrics` covers all individuals
   instead of the first object's.**
3. Stop when `loc.all` differs between objects (F3). **Consequence: calls
   that now combine opposite allele codings stop with an error.**
4. Stop when SNP and SilicoDArT objects are mixed (F4). **Consequence:
   calls that now combine them stop with an error.**
5. Argument handling and docs (F5). **Consequence: a non-genlight,
   non-NULL argument stops with an error instead of being ignored.**

Callers checked: dartR.base `gl.join` (overwrites `ind.metrics`,
`loc.metrics` and flags after `rbind`), `gl.impute` (`rbind(NULL, pop)`,
then restores metadata from the input), `gl.report.pa`,
`gl2faststructure`; dartR.captive `gl.sim.relatedness`; dartR.sim
`gl.sim.emigration`. Their tests are to be re-run in Phase C.

## Coverage

- Standards walk: FS, DOC, DAT, API — run (S3 method: FS2-FS9 entry
  structure does not apply).
- Spec: both paths on platypus.gl; FBM path via `gl.gen2fbm` — run.
- Locus-order, allele-coding, mixed-ploidy and list-input cases — run.
- Google Group: SKIPPED — not available: no browser session.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (genotypes change) stated in the question |
| 2 | approved | Luis | consequence (metadata now carried) stated in the question |
| 3 | approved | Luis | now errors |
| 4 | approved | Luis | now errors |
| 5 | approved | Luis | now errors |

## Outcome

All five changes applied in `R/utils.dartR.class.def.r`; metadata built by
one helper (`.rbind_other`) used by both paths, flags and history by
`.rbind_finish`. The history entry is `rbind.dartR(n.objects = k)`:
`match.call()` under `do.call(rbind, ...)` embeds the whole objects
(6.3 MB for two 2x3 objects), so a compact call is stored instead.

Snapshot result: the baseline failed in exactly four places, each mapped to
an approved change — empty `@other` (2), misaligned genotypes (1), silent
allele-coding mismatch (3), silent SNP + SilicoDArT (4). No other diff.

Tests: test-rbind.dartR.R 24/24 (in-memory and FBM paths give identical
genotypes and 10 `ind.metrics` rows; locus-order case equals the input
exactly). Callers: test-gl.join.R 25/25, test-gl.report.pa.R 13/13,
test-gl2faststructure.R 17/17, test-utils.dartR.class.def.R 17/17,
test-gl.impute.R 74/75 — the failure (a `verbose = 0` silence check) fails
identically on unmodified origin/dev and is the `gl.alf` leak fixed by
PR #420. `devtools::document()` run. Full R CMD check left to CI.

PR: #422

```json
{
  "function": "rbind.dartR",
  "package": "dartR.base",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "1dd6a3a",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "API2", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["Google Group: no browser session"],
  "status": "pr-open",
  "pr": 422
}
```
