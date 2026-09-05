# Review: gl2structure (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2structure.R (snapshot captured pre-review; 16 assertions passing)

## Verdicts

**Standards: Needs work** — sound preamble and parameter validation, but the
datatype gate admits SilicoDArT (a named DAT7 offender) and the file write
appends instead of overwriting.
**Spec: Needs work** — the two-row-per-individual layout, 1/2 allele codes
and -9 missing code all verified correct against the genotype matrix, but
two realistic paths (re-runs, duplicated locus names) silently produce files
STRUCTURE will misparse.

What works well: genotype coding verified cell-by-cell — dosage 0 -> 1/1,
1 -> 1/2, 2 -> 2/2, NA -> -9/-9 — and the `ind.names`/`add.columns` length
checks fail fast with clear messages.

## Findings

**F1 [HIGH, confidence: high] — unconditional append corrupts the output on re-run (no catalogue rule: file-overwrite contract)**
`R/gl2structure.r:155-163` — `write.table(..., append = TRUE)` never
truncates. With `export.marker.names = TRUE` the preceding `cat(...,
file = outfilespec)` (148-152) happens to truncate first, but with
`export.marker.names = FALSE` nothing does.
Failure scenario: verified — two successive
`gl2structure(x, export.marker.names = FALSE)` calls to the same path
produce 32 data rows for 8 individuals; STRUCTURE reads every individual
twice and the run is silently wrong.
Proposed change: truncate first (e.g. `append = export.marker.names`, or an
explicit `file.create()`/single `write.table` with `append = FALSE` when no
header is written).

**F2 [MEDIUM, confidence: high] — SilicoDArT admitted, exported as fabricated diploid genotypes (DAT7)**
`R/gl2structure.r:64` — `utils.check.datatype(x, verbose = verbose)` keeps
the default `accept`; this function is one of the two named offenders in
the DAT7 rule.
Failure scenario: verified — `gl2structure(testset.gs[1:4, 1:10])` writes a
syntactically valid STRUCTURE file in which presence/absence scores appear
as diploid genotypes (0 -> 1/1, 1 -> 1/2); STRUCTURE runs happily on noise.
Proposed change: `accept = "SNP"`.

**F3 [MEDIUM, confidence: high] — duplicated locus names silently collapse genotype columns (DAT2 analogue: output columns must track loci 1:1)**
`R/gl2structure.r:140-145` — genotype columns are assigned by locus name
(`StructTab[[dimnames(genmat)[[2]][i]]] <- ...`), so a repeated name
overwrites the earlier column.
Failure scenario: verified — with 2 of 15 locus names duplicated, the
header row lists 15 markers but data rows carry 14 genotype columns, and
the first duplicate's genotypes are replaced by the second's; STRUCTURE
misaligns every column after the collision.
Proposed change: assign by position (`StructTab[[ncol(StructTab) + 1]]`)
and set names afterwards, or stop with a clear message on duplicated
`locNames(x)`.

**F4 [LOW, confidence: high] — NULL return is visible (FS10, VRB5)**
`R/gl2structure.r:177` — `return(NULL)` prints `NULL` on every unassigned
call, including at `verbose = 0` (verified via `withVisible`).
Proposed change: `invisible(NULL)`.

**F5 [LOW, confidence: high] — roxygen drift (DOC1, DOC5; DOC2, DOC7 proposed rules)**
`R/gl2structure.r:5-39` — `@description`/`@param x` claim the object must
contain "location data, lat longs", which the function never touches;
`@param ind.names` says "defaults to ind.names(x)" (the accessor is
`indNames`); no `@details` tag; the `verbose` default clause is the
outdated wording (DOC2); the `@author` block lacks the labelled
`Author(s): ... Custodian: ...` structure (DOC7).
Failure scenario: users pre-filter for coordinate completeness they do not
need; help text misnames an accessor.
Proposed change: drop the lat/long claim, fix the accessor name, adopt the
DOC2 clause and DOC7 author labels, add a minimal `@details`; run
`devtools::document()` (DOC4).

## Proposed changes

1. Truncate the output file before writing when no marker-name header is
   emitted (F1).
2. Restrict the datatype gate with `accept = "SNP"` (F2).
   **Consequence: SilicoDArT callers now error instead of receiving a
   fabricated diploid file.**
3. Build genotype columns positionally or reject duplicated locus names
   loudly (F3).
4. Replace `return(NULL)` with `invisible(NULL)` (F4).
5. Roxygen conformance pass + `devtools::document()` (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: layout/coding/missing-code verified cell-by-cell against
  `as.matrix` on testset.gl subsets; re-run append path; SilicoDArT
  admission; duplicate-name path; ind.names/add.columns validation — run
  empirically
- Parsing the file with STRUCTURE itself: SKIPPED — no STRUCTURE binary on
  the review machine; layout checked against the documented two-row format
- `ploidy != 2` path: SKIPPED — no packaged non-diploid genlight fixture
- FBM path (DAT6): SKIPPED — no FBM fixture for converters

## Approval (Phase B)

Blanket approval by Arthur Georges, 2026-09-05, via the formal approval
boxes: all findings at every severity approved, explicitly acknowledging
that converted outputs change where they were wrong; the DAT7 class fix
(fatal `accept = "SNP"` gate wherever SilicoDArT is silently admitted) is
approved campaign-wide and applies here as change 2.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05 |
| 2 | Approved | Arthur Georges | 2026-09-05; DAT7 fatal gate, SilicoDArT callers now error |
| 3 | Approved | Arthur Georges | 2026-09-05; output-change consequence acknowledged |
| 4 | Approved | Arthur Georges | 2026-09-05 |
| 5 | Approved | Arthur Georges | 2026-09-05 |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2structure` (base upstream/dev,
ddaed27). All five changes applied:

1. F1 — `write.table(append = export.marker.names)`: the data write
   appends only when the marker-name header has just truncated the file;
   a headerless re-run now overwrites instead of doubling the file.
2. F2 — `utils.check.datatype(x, accept = "SNP", ...)`: SilicoDArT is a
   fatal error instead of a fabricated diploid export.
3. F3 — genotype columns assigned positionally
   (`StructTab[[ncol(StructTab) + 1]]`, names set afterwards), so
   duplicated locus names no longer collapse columns.
4. F4 — `invisible(NULL)`.
5. F5 — roxygen: lat/long claim dropped, `indNames(x)` accessor named
   correctly, DOC2 verbose clause, DOC7 author labels, minimal
   `@details` added, house tag order; `devtools::document()` run
   (gl2structure.Rd only).

Verification: baseline characterization test passes (17 assertions; the
flipped expectations map to F1 — re-run no longer doubles the file, F2 —
SilicoDArT now fatal, F3 — 15 genotype columns retained under duplicate
names, F4 — invisible return). End-to-end run on testset.gl at
`verbose = 3` clean: 549 lines = 1 header + 2 x 274 individuals.
Sibling grep across the 8 clones: the only external reference,
`dartR.captive/R/gl2colony.r:212`, calls a same-file local
`gl2structure` defined at line 431 which shadows dartR.base's export —
no caller breaks. NEWS.md entry added.

```json
{
  "function": "gl2structure",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "none:file-overwrite-contract", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 5}
  ],
  "coverage_skipped": ["STRUCTURE parse test: no STRUCTURE binary", "ploidy != 2: no fixture", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": null
}
```
