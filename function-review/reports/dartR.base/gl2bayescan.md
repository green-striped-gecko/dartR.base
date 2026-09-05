# Review: gl2bayescan (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2bayescan.R (snapshot captured pre-review; 14 assertions passing)

## Verdicts

**Standards: Needs work** — preamble and file handling follow the house pattern (`gl.check.wd`, `file.path`, flag start/end), but SilicoDArT slips through the datatype check and the `sink()` redirection has no error protection.

**Spec: Needs work** — the GESTE codominant layout and the allele counts are correct on SNP data (verified by hand for population 1, locus 1: `1 22 2 22 0` matches `2*nobs`, alt count, ref count from the genotype matrix), but 815 of 23405 data rows on testset.gl are zero-sample rows (`j 0 2 0 0`) for populations with no genotyped individual at a locus, and BayeScan's tolerance of zero gene copies is not established.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted into diploid gene-copy math (DAT7)**
`R/gl2bayescan.r:51` — permissive `utils.check.datatype` default; the writer computes gene copies as `2 * nobs` (lines 93-96), which presumes ploidy 2.
Failure scenario: `gl2bayescan(testset.gs, ...)` runs without error (verified 2026-09-05) and writes presence/absence tag counts doubled into a codominant BayeScan file; the docs offer BayeScan's dominant-marker path for such data, which this output is not.
Proposed change: `accept = "SNP"` in the datatype call.

**F2 [MEDIUM, confidence: medium] — zero-sample rows written for all-missing population x locus cells (io-family missing-data trap; no catalogue rule fits — flagged as a principle)**
`R/gl2bayescan.r:91-98` — when a locus has no genotyped individual in a population, the row `<j> 0 2 0 0` is written. On testset.gl this happens 815 times (verified).
Failure scenario: BayeScan reads a population sample of size 0 for those loci; behaviour ranges from undefined estimates to a crash depending on version. The user gets no warning that their file contains empty samples.
Proposed change: warn at `verbose >= 1` with the count of zero-sample cells (VRB4 spirit), and document that heavy filtering should precede export; optionally offer to drop affected loci.

**F3 [LOW, confidence: high] — `sink()` without `on.exit` protection (no catalogue rule fits — flagged as a principle)**
`R/gl2bayescan.r:76-102` — the console is redirected between two bare `sink()` calls.
Failure scenario: any error inside the loop (e.g. a write failure on a full disk) leaves the session's console redirected; every later message vanishes into the file.
Proposed change: `sink(outfilespec); on.exit(sink(), add = TRUE)` — the closing `sink()` may stay.

**F4 [LOW, confidence: high] — visible NULL return (FS10)**
`R/gl2bayescan.r:116` — `return(NULL)` is visible; an unassigned call prints `NULL` at the console (verified with `withVisible`). Sibling `gl2eigenstrat` returns `invisible(NULL)`.
Proposed change: `return(invisible(NULL))`.

**F5 [LOW, confidence: high] — roxygen order and completeness (DOC1, DOC2)**
`R/gl2bayescan.r:1-31` — no `@details`; `@return` after `@export`; verbose clause is the outdated "[default 2 or as specified...]" wording; "BAyescan" typo in the description.
Proposed change: roxygen pass to the house template.

**F6 [LOW, confidence: high] — author block names a custodian only (DOC7) (proposed rule)**
`R/gl2bayescan.r:18-19` — "Custodian: Luis Mijangos" with no `Author(s):` line.
Proposed change: add the `Author(s):` part per DOC7.

## Proposed changes

1. Add `accept = "SNP"` to the datatype check (F1).
2. Count zero-sample population x locus cells and warn at `verbose >= 1`; document the recommendation to filter on call rate before export (F2).
3. Add `on.exit(sink(), add = TRUE)` after opening the sink (F3).
4. Return `invisible(NULL)` (F4).
5. Roxygen pass: `@details`, tag order, DOC2 verbose text, typo, DOC7 author structure (F5, F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: GESTE layout vs Foll & Gaggiotti codominant format on testset.gl — run;
  allele counts hand-verified for pop 1 locus 1; `sum` from `gl.allele.freq(..., by = 'popxloc')`
  confirmed to be the raw alternate-allele count (unaffected by `percent = TRUE`)
- Locus-order consistency across populations — verified: `gl.allele.freq` orders by
  `loc_order` then `popn`; the stable re-sort on `popn` preserves locus order within pops
- SilicoDArT admission probe on testset.gs — run (admitted; F1)
- verbose = 0 silence via `capture.output` — run (0 lines)
- Genome-position lens: n/a — the format carries no chromosome/position fields
- FBM path (DAT6): SKIPPED — no FBM fixture
- BayeScan executable round-trip: SKIPPED — binary unavailable; zero-sample behaviour (F2) therefore unconfirmed

## Approval (Phase B)

All findings at every severity approved 2026-09-05 by Arthur Georges via the
formal approval boxes, as part of the blanket class approval for the io
converter batch: converted outputs change where they were wrong, and the DAT7
class fix (fatal `accept = "SNP"` gate wherever SilicoDArT is silently
admitted) applies.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | DAT7 fatal gate (class approval) |
| 2 | Approved | Arthur Georges | Warning at verbose >= 1 + docs; rows still written |
| 3 | Approved | Arthur Georges | Depth-guarded on.exit |
| 4 | Approved | Arthur Georges | |
| 5 | Approved | Arthur Georges | |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2bayescan` (base `upstream/dev` =
ddaed27). All five changes applied. For change 3 the guard records the sink
depth before opening and pops back to it in `on.exit()`, so an error during
writing restores the console without disturbing any enclosing diversion
(e.g. `capture.output`); verified empirically -- `sink.number()` is 0 after a
forced write failure. Baseline characterization test updated for the three
approved behaviour changes (F1 SilicoDArT rejection, F2 warning, F4
invisible return); 14 assertions pass. End-to-end run on `testset.gl` at
`verbose = 3` verified: identical GESTE records, and the F2 warning reports
the 815 zero-sample cells. Sibling-clone grep across the 8 dartR-verse
packages: no callers of `gl2bayescan` -- all clear. NEWS.md entry added.

```json
{
  "function": "gl2bayescan",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "medium", "rule": "principle:io-missing-data", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "principle:sink-on.exit", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "BayeScan executable round-trip: binary unavailable"],
  "status": "pr-open",
  "pr": null
}
```
