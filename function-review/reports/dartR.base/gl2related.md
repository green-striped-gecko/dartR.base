# Review: gl2related (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT)
- Baseline: tests/testthat/test-gl2related.R (snapshot captured pre-review; 13 assertions, all passing against current code)

## Verdict

**Standards: Ready** — the cleanest function in this batch: verbosity gating is correct (`verbose = 0` verified fully silent), the visible data-frame return matches `@return`, and the vectorised recoding is tidy; only roxygen wording items remain.
**Spec: Needs work** — the two-column allele coding is verified correct (het 1/3, hom-ref 1/1, hom-alt 3/3, missing 0/0), but the exported file wraps individual names in quotes, and SilicoDArT data is admitted into diploid coding.

## Findings

**F1 [MEDIUM, confidence: medium] — individual names quoted in the output file (DOC5)**
`R/gl2related.r:92-98` — `write.table` runs with its default `quote = TRUE`, so the first column is written as `"AA010915"` (verified in the file). The description says the file "can be used and loaded into coancestry".
Failure scenario: COANCESTRY (the standalone Fortran program, the file's documented consumer) reads the quote characters as part of the individual identifier or rejects the line; round-tripping through R's `read.table` masks the defect, which is why it survives. Confidence is medium because COANCESTRY's exact parse behaviour was not tested here.
Proposed change: add `quote = FALSE` to the `write.table` call.

**F2 [MEDIUM, confidence: high] — SilicoDArT admitted into diploid allele coding (DAT7)**
`R/gl2related.r:70` — `utils.check.datatype(x, verbose = verbose)` keeps the permissive default while the recoding assumes 0/1/2 dosage. Verified: testset.gs converts, with every presence score (1) emitted as a heterozygote pair `1 3`.
Failure scenario: a SilicoDArT user feeds COANCESTRY/related a file in which presence/absence is silently reinterpreted as heterozygosity — relatedness estimates are meaningless, with no error at any point. The roxygen says "SNP data" only.
Proposed change: `accept = "SNP"` (the recurring DAT7 defect class).

**F3 [LOW, confidence: high] — output coding undocumented (DOC5)**
`R/gl2related.r:5-13,49` — neither `@description` nor `@return` states the file/data-frame coding (alleles as 1 and 3, missing as 0, two columns per locus, tab-separated, no header). A user validating the export against COANCESTRY's manual has to reverse-engineer it from the source.
Proposed change: add a `@details` paragraph stating the coding.

**F4 [LOW, confidence: high] — roxygen wording drift (DOC2, DOC7 (proposed rule), DOC1)**
`R/gl2related.r:24-29` — the `verbose` default clause ("[default 2, unless specified using gl.set.verbosity]") is not the ratified DOC2 wording; the author block reads "Bernd Gruber (bugs? Post to ...)" with neither an `Author(s):` nor a `Custodian:` label (DOC7 — default fix: `Author(s): Bernd Gruber. Custodian: Bernd Gruber`); `@return` sits after `@export` (DOC1 old order — flag the file).

## Proposed changes

1. Add `quote = FALSE` to the `write.table` call (F1). **Consequence: file output changes for every caller (quotes disappear from the name column).**
2. Restrict datatype with `accept = "SNP"` (F2). **Consequence: SilicoDArT input errors instead of converting.**
3. Document the allele coding in `@details` (F3).
4. Align verbose/author/return roxygen per DOC2/DOC7/DOC1; run `devtools::document()` (F4, DOC4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plots), STY — run
- Spec: converter executed on a testset.gl subset; returned data frame and written file verified cell-by-cell against `as.matrix` for all four genotype states; tab separation, header absence, and row count checked — run
- `save = FALSE` and `verbose = 0` silence: run — clean
- SilicoDArT admission: run — converts silently (F2)
- Load of the exported file into COANCESTRY or `related::coancestry`: SKIPPED — neither program installed; `related` is R-Forge-only. F1 confidence capped accordingly
- FBM path (DAT6): SKIPPED — no FBM fixture; converter densifies via `as.matrix` by design

## Approval (Phase B)

All findings at every severity approved by Arthur Georges, 2026-09-05, via
the formal approval boxes -- a blanket class approval covering this batch,
explicitly acknowledging that converted outputs change where they were
wrong, and ratifying the DAT7 fatal-gate policy (applied here as change 2).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05; file output changes for every caller |
| 2 | Approved | Arthur Georges | 2026-09-05; DAT7 policy: fatal accept = "SNP" gate |
| 3 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 4 | Approved | Arthur Georges | 2026-09-05, blanket class approval |

## Outcome (Phase C)

All four changes applied on branch `review-gl2related` (base
`upstream/dev` = ddaed27), 2026-09-05. Verification: baseline
characterization test updated -- every diff maps to an approved finding
(F1 unquoted name column, F2 SilicoDArT rejection) -- 13 assertions
pass, including the cell-by-cell allele-coding checks unchanged.
End-to-end on a testset.gl subset at verbose = 3: clean; the written
file starts with the bare individual name. Sibling grep across the 8
clones: no package-code callers of `gl2related` -- all clear.
PR: https://github.com/green-striped-gecko/dartR.base/pull/PRNUM

```json
{
  "function": "gl2related",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "ready",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "applied", "change": 4}
  ],
  "coverage_skipped": ["COANCESTRY/related round-trip: programs not installed", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": null
}
```
