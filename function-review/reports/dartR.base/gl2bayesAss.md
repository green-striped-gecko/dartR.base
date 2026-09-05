# Review: gl2bayesAss (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2bayesAss.R (snapshot captured pre-review; 14 assertions passing)

## Verdicts

**Standards: Needs work** — the FS skeleton (verbosity, wd check, flag start/end, datatype check) is in place and `data.table` is a full NAMESPACE import needing no guard, but the datatype check admits SilicoDArT into diploid-only math and the roxygen deviates from the house template.

**Spec: Needs work** — on SNP data the output matches the BA3/BA3-SNPs row format exactly (verified against the genotype matrix: `indiv pop locus a1 a2`, missing as `0 0`, het as `1 2`); the `ploidy` parameter, however, is a guard dressed as an option, and presence/absence data passes through unrejected.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted into diploid-only recode (DAT7)**
`R/gl2bayesAss.r:55` — `utils.check.datatype(x, verbose = verbose)` uses the permissive default; the recode table at line 65 is defined for diploid dosages 0/1/2 only.
Failure scenario: `gl2bayesAss(testset.gs, ...)` runs without error (verified 2026-09-05) and writes presence/absence scores as fake diploid genotypes (absent tag `0 -> 1 1` homozygote, present tag `1 -> 1 2` heterozygote), producing a plausible-looking but meaningless BA3 input. The `ploidy != 2` check at line 58 does not catch it because `ploidy` is a user argument, not the object's ploidy.
Proposed change: `accept = "SNP"` in the datatype call.

**F2 [MEDIUM, confidence: high] — `ploidy` parameter is inert but documented as settable (DOC5, proposed rule)**
`R/gl2bayesAss.r:38,58` — the roxygen says "Set the ploidy [defaults 2]" yet the only accepted value is 2; any other value stops. The parameter configures nothing.
Failure scenario: a user with haploid or polyploid data reads the docs, passes `ploidy = 1`, and gets a refusal for a parameter the signature invited them to set.
Proposed change: either document it honestly as a confirmation guard ("must be 2; the function only supports diploid data") or remove the parameter (removal is API-affecting — API2 requires a loud error for old callers).

**F3 [LOW, confidence: high] — plain `stop()` without the error helper (VRB2)**
`R/gl2bayesAss.r:58` — `stop("This function only caters for ploidy=2")` bypasses `error()` and the "Fatal Error:" prefix used package-wide.
Failure scenario: the message prints uncoloured and off-pattern; scripted matching on "Fatal Error" misses it.
Proposed change: `stop(error("Fatal Error: this function only caters for ploidy = 2\n"))`.

**F4 [LOW, confidence: high] — roxygen order and completeness (DOC1)**
`R/gl2bayesAss.r:1-36` — `@return` sits last, after `@export`; there is no `@details`; `@param ploidy` uses "[defaults 2]" rather than "[default 2]".
Failure scenario: none at runtime; the file diverges from the ratified header order.
Proposed change: reorder to the house sequence and add a short `@details` (even one line on the BA3 row format).

**F5 [LOW, confidence: high] — nonstandard verbose default clause (DOC2)**
`R/gl2bayesAss.r:15-17` — "[default 2 or as specified using gl.set.verbosity]" instead of the ratified NULL-cascade wording.
Proposed change: adopt the DOC2 standard text.

**F6 [LOW, confidence: high] — author block names a custodian only (DOC7) (proposed rule)**
`R/gl2bayesAss.r:19-20` — "Custodian: Carlo Pacioni" with no `Author(s):` line.
Proposed change: `Author(s): Carlo Pacioni. Custodian: Carlo Pacioni -- Post to \url{...}`.

## Proposed changes

1. Add `accept = "SNP"` to the `utils.check.datatype` call (F1).
2. Re-document `ploidy` as a diploid-confirmation guard, or remove it with a loud unused-argument error (F2). **Consequence if removed: callers passing `ploidy=` error explicitly.**
3. Replace the plain `stop()` with the `stop(error("Fatal Error: ..."))` idiom (F3).
4. Roxygen pass: house tag order, add `@details`, fix "[defaults 2]", DOC2 verbose text, DOC7 author/custodian structure (F4, F5, F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: output format vs BA3/BA3-SNPs row layout on testset.gl — run; full recode
  mapping (0/1/2/NA -> allele pairs) verified against `as.matrix(x)` in the baseline test
- SilicoDArT admission probe on testset.gs — run (admitted; F1)
- verbose = 0 silence via `capture.output` — run (0 lines)
- Genome-position lens: n/a — the function writes no chromosome/position fields
- FBM path (DAT6): SKIPPED — no FBM fixture; `as.matrix(x)` densifies by design in this family
- BA3 executable round-trip: SKIPPED — BayesAss binary not available; format checked against Mussmann et al. (2019) immanc layout

## Approval (Phase B)

All findings at every severity approved 2026-09-05 by Arthur Georges via the
formal approval boxes, as part of the blanket class approval for the io
converter batch: converted outputs change where they were wrong, and the DAT7
class fix (fatal `accept = "SNP"` gate wherever SilicoDArT is silently
admitted) applies.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | DAT7 fatal gate (class approval) |
| 2 | Approved | Arthur Georges | Documented as diploid-confirmation guard; parameter retained |
| 3 | Approved | Arthur Georges | |
| 4 | Approved | Arthur Georges | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2bayesAss` (base `upstream/dev` =
ddaed27). All four changes applied; for change 2 the documentation option was
taken (`ploidy` retained, documented as a confirmation guard) — no API
removal. Baseline characterization test updated for the two approved
behaviour changes (F1 SilicoDArT rejection, F3 error-message format); all 13
assertions pass. End-to-end run on `testset.gl` at `verbose = 3` verified
(206870 rows; SNP recode unchanged). Sibling-clone grep across the 8
dartR-verse packages: no callers of `gl2bayesAss` — all clear. NEWS.md entry
added.

```json
{
  "function": "gl2bayesAss",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "applied", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 4}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "BA3 executable round-trip: binary unavailable"],
  "status": "pr-open",
  "pr": null
}
```
