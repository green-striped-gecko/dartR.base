# Review: gl2demerelate (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2demerelate.R (snapshot captured pre-review; 18 assertions passing)

## Verdicts

**Standards: Needs work** — the body is clean and the FS skeleton correct, but the datatype check admits SilicoDArT into a diploid two-allele recode, and the roxygen is missing required tags outright.

**Spec: Ready** — the output matches the Demerelate input contract on SNP data: `Sample-ID` / `Population` then two adjacent allele columns per locus, coding verified against the genotype matrix (`0 -> 1/1`, `1 -> 1/2`, `2 -> 2/2`, NA preserved), locus names sanitized of `-`, `|`, `/`.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted into diploid allele recode (DAT7)**
`R/gl2demerelate.r:35` — permissive `utils.check.datatype` default; the recode at lines 51-59 assumes dosages 0/1/2.
Failure scenario: `gl2demerelate(testset.gs, ...)` runs without error (verified 2026-09-05) and converts presence/absence scores into fake genotypes — every present tag (1) becomes a "heterozygote" 1/2 — feeding Demerelate meaningless relatedness input.
Proposed change: `accept = "SNP"` in the datatype call.

**F2 [MEDIUM, confidence: high] — required roxygen tags absent (DOC1)**
`R/gl2demerelate.r:1-19` — no `@description` and no `@details` at all; the header jumps from `@family` to `@param`; `@return` sits last after `@export`; the DOC2 verbose clause is the outdated wording and lacks its final period.
Failure scenario: the man page renders with the title doubling as the description and no statement of the coding scheme — the one thing a Demerelate user needs to know.
Proposed change: add `@description` (one paragraph) and a `@details` documenting the 1/2 allele coding and NA-as-missing convention; house tag order; DOC2 text.

**F3 [LOW, confidence: high] — dead assignment before the preamble (FS2, STY1)**
`R/gl2demerelate.r:23` — `x_temp <- x` executes before SET VERBOSITY and is overwritten unconditionally at line 61.
Failure scenario: none at runtime; it misleads the reader into hunting for a use.
Proposed change: delete the line.

**F4 [LOW, confidence: high] — author block names a custodian only (DOC7) (proposed rule)**
`R/gl2demerelate.r:11-12` — "Custodian: Luis Mijangos" with no `Author(s):` line.
Proposed change: add the `Author(s):` part per DOC7.

## Proposed changes

1. Add `accept = "SNP"` to the datatype check (F1).
2. Roxygen pass: add `@description`/`@details` (documenting the allele coding), house order, DOC2 verbose text, DOC7 author structure (F2, F4).
3. Delete the dead `x_temp <- x` assignment (F3).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP (no external packages used), PLT (n/a), STY — run
- Spec: full coding check against `as.matrix(testset.gl)` — run (0/1/2/NA all verified);
  `_1`/`_2` column adjacency after the alphabetical sort — verified on a 30-locus subset;
  name sanitation (`-`, `|` to `_`; `/` removed) — verified
- Missing data: NA is preserved (the commented-out zero-recode at line 63 is inert);
  Demerelate reads NA as missing — checked, nothing found
- SilicoDArT admission probe on testset.gs — run (admitted; F1)
- verbose = 0 silence via `capture.output` — run (0 lines)
- Genome-position lens: n/a — no chromosome/position fields in the output
- FBM path (DAT6): SKIPPED — no FBM fixture; double `as.matrix(x)` densification noted
  (two full copies held simultaneously) but not flagged absent an FBM-scale fixture
- Demerelate package round-trip: SKIPPED — Demerelate was removed from CRAN; format
  checked against its archived input documentation

## Approval (Phase B)

All findings at every severity approved 2026-09-05 by Arthur Georges via the
formal approval boxes, as part of the blanket class approval for the io
converter batch: converted outputs change where they were wrong, and the DAT7
class fix (fatal `accept = "SNP"` gate wherever SilicoDArT is silently
admitted) applies.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | DAT7 fatal gate (class approval) |
| 2 | Approved | Arthur Georges | |
| 3 | Approved | Arthur Georges | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2demerelate` (base `upstream/dev` =
ddaed27). All three changes applied. Baseline characterization test updated
for the single approved behaviour change (F1 SilicoDArT rejection); 18
assertions pass -- the SNP dataframe (shape, coding, name sanitation,
adjacency, NA handling) is byte-identical to the pre-review snapshot.
End-to-end run on `testset.gl` at `verbose = 3` verified (274 x 1512
dataframe). Sibling-clone grep across the 8 dartR-verse packages: no callers
of `gl2demerelate` -- all clear. NEWS.md entry added.

```json
{
  "function": "gl2demerelate",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 2}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Demerelate round-trip: package archived on CRAN"],
  "status": "pr-open",
  "pr": null
}
```
