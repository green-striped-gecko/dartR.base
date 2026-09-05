# Review: gl2hiphop (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2hiphop.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the body is short and clean, but SilicoDArT data
is admitted where the recode only makes sense for SNP dosages.
**Spec: Ready** — for SNP input the output matches the hiphop package's
expected layout exactly (verified): a data frame with individuals as rows
(row names = individual IDs, which hiphop keys on), loci as columns, and the
hiphop coding 0 = homozygote, 1 = other homozygote, 2 = heterozygote — i.e.
dartR 1 -> 2 and 2 -> 1 — with NA preserved cell-for-cell.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted (DAT7)**
`R/gl2hiphop.r:51` — `utils.check.datatype(x, verbose = verbose)` uses the
default `accept`, so presence/absence data passes a recode written for
0/1/2 dosages.
Failure scenario: verified — `gl2hiphop(testset.gs)` returns a data frame of
0s and 2s: every presence call (1) is relabelled 2, "heterozygote" in hiphop
coding. Fed to `hiphop::hothiphop()` this scores parentage mismatches on
meaningless pseudo-genotypes, silently. The roxygen says "the genlight
object containing the SNP data"; the code does not enforce it. This is the
recurring defect class already logged for `gl2structure`, `gl2vcf` and
others.
Proposed change: `accept = "SNP"`.

**F2 [LOW, confidence: high] — roxygen gaps (DOC1)**
`R/gl2hiphop.r:1-37` — no `@details` tag (the reference-vs-hiphop coding
swap deserves one line); `@return` sits after `@export` (pre-2026-08-27
order).
Failure scenario: none at runtime; a user comparing dartR and hiphop
encodings has to read the source to learn that 1 and 2 are swapped.
Proposed change: add `@details` stating the coding map
(0 -> 0, 1 -> 2, 2 -> 1, NA -> NA), move `@return` before `@author`.

**F3 [LOW, confidence: high] — author block lacks Custodian line (DOC7) (proposed rule)**
`R/gl2hiphop.r:19` — `@author` names Luis Mijangos with no `Custodian:`
label.
Proposed change: add `Custodian: Luis Mijangos` per the DOC7 default.

**F4 [INFO, confidence: high] — character round-trip recode and a dplyr dependency for as.numeric (STY3)**
`R/gl2hiphop.r:55-61` — the 1 <-> 2 swap coerces the numeric matrix to
character via a `"het"` sentinel, then converts back with
`dplyr::mutate_all(as.numeric)`.
Failure scenario: none observed — values, dimnames and NAs all survive
(verified). Noted because a numeric one-liner
(`m[m == 1] <- 3; m[m == 2] <- 1; m[m == 3] <- 2` or `c(0, 2, 1)[m + 1]`)
avoids the type churn and the dplyr import.
Proposed change: optional simplification only.

## Proposed changes

1. Restrict datatype with `accept = "SNP"` (F1).
   **Consequence: SilicoDArT callers that previously received 0/2
   pseudo-genotypes now receive an error.**
2. Roxygen: add `@details` with the coding map, reorder `@return`, add the
   Custodian line (F2, F3).
3. Optional: replace the character-sentinel recode with a numeric recode and
   drop the `mutate_all` import (F4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no file output so
  FS7 not applicable; no plot bundle; dplyr is in Depends so DEP1 not
  applicable; no history append needed — the return is a data frame).
- Spec: recode verified value-for-value against `as.matrix(testset.gl)`
  (0 -> 0, 1 -> 2, 2 -> 1, NA -> NA); row/column name survival verified
  (hiphop matches individuals by row name); hiphop's documented coding
  checked against the package's published input description — run.
- `verbose = 0` silence: verified (zero captured lines) — run.
- End-to-end run through `hiphop::hothiphop()`: SKIPPED — hiphop not
  installed on the review machine; layout checked against its documented
  input contract instead.
- FBM path (DAT6): SKIPPED — no FBM fixture. `as.matrix(x[,])` densifies;
  same team-level consideration as the other converters.
- Note (not a finding): the 3 all-NA loci of testset.gl are retained as
  all-NA columns; hiphop treats NA as unscored, so this is harmless.

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 2 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 3 | Approved | Arthur Georges | 2026-09-05, via approval boxes |

All findings at every severity approved 2026-09-05 (blanket class
approval, explicitly acknowledging that converted outputs change where
they were wrong, and including the DAT7 fatal `accept = "SNP"` gate —
applied here as change 1).

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2hiphop` (base `upstream/dev`,
ddaed27). PR: green-striped-gecko/dartR.base#354.

- F1: `accept = "SNP"` — SilicoDArT input now errors instead of
  returning 0/2 pseudo-genotypes.
- F2/F3: `@details` added stating the coding map (0 -> 0, 1 -> 2,
  2 -> 1, NA preserved); `@return` moved to the house position;
  `Author(s)`/`Custodian` block; re-documented.
- F4: the character-sentinel recode replaced with the numeric lookup
  `c(0, 2, 1)[m + 1]`; the `mutate_all` import dropped from NAMESPACE.
  The `%>%` import is retained in this file although the function no
  longer uses it: `gl.map.interactive`, `gl.pcoa.plot` and
  `utils.heatmap` use the pipe at run time and gl2hiphop's roxygen was
  the only `@importFrom` declaring it — removing it would break them.
  Follow-up: move the declaration to a package-level file.

Verification: all 11 baseline assertions pass, including the
value-for-value recode check against `as.matrix(testset.gl)` and
`verbose = 0` silence; the single updated expectation is the SilicoDArT
rejection, per F1. End-to-end `gl2hiphop(testset.gl, verbose = 3)`
returns the 274 x 755 data frame with values 0/1/2 and NA. Sibling
caller grep across the 8 clones: no callers of `gl2hiphop` — all-clear.
NEWS entry added.

Note: the NAMESPACE change is scoped by hand to the `mutate_all`
removal only; a full `devtools::document()` also wants to add
`importFrom(grDevices, heat.colors)` and drop `importFrom(plyr, count)`
— pre-existing drift from other files' headers (DOC4), left for their
own reviews.

```json
{
  "function": "gl2hiphop",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 2},
    {"id": "F4", "severity": "INFO", "confidence": "high", "rule": "STY3", "status": "applied", "change": 3}
  ],
  "coverage_skipped": ["hiphop end-to-end: hiphop not installed", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": 354
}
```
