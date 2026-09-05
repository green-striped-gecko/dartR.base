# Review: gl2treemix (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2treemix.R (snapshot captured pre-review; 130 assertions passing)

## Verdicts

**Standards: Needs work** — clean structure and preamble, but the datatype
gate admits SilicoDArT and the gz write is not sink-safe.
**Spec: Ready** — the output was verified cell-by-cell against manually
computed allele counts across the full test matrix: gzipped, header of
population names whose order matches the columns, one `ref,alt` row per
locus, `0,0` for an all-missing population-by-locus cell. For SNP data the
file is exactly what the treemix manual asks for.

What works well: the fragile-looking index walk over `gl.allele.freq`'s
popxloc table (`k:(k + nPop - 1)`) is in fact correct — the table is sorted
locus-major with populations in factor-level order, and a full-matrix
recomputation matched every cell.

## Findings

**F1 [MEDIUM, confidence: high] — SilicoDArT admitted; each haploid individual counted as two allele copies (DAT7, DAT1)**
`R/gl2treemix.r:54,59` — `utils.check.datatype` keeps the default `accept`,
and `freq$ref <- freq$nobs * 2 - freq$sum` hard-codes diploidy.
Failure scenario: verified — `gl2treemix(testset.gs[1:6, 1:5])` writes a
file where the allele-copy total per cell is `2 * n_scored` for ploidy-1
data; treemix's drift model sees twice the real number of sampled
chromosomes and understates drift throughout.
Proposed change: `accept = "SNP"` (or, if presence/absence support is
intended, compute `ref <- nobs - sum` under a documented haploid branch —
custodian call; the doubling comes from `gl.allele.freq`'s deliberate
"biallelic, no heterozygotes" recoding, so intent is ambiguous).

**F2 [MEDIUM, confidence: medium] — unguarded sink on an inline gz connection (no catalogue rule; gl2fasta F3 precedent)**
`R/gl2treemix.r:78-90` — `sink(gzfile(outfilespec))` has no
`on.exit(sink())`, and the connection created inline is never explicitly
closed.
Failure scenario: an error or interrupt inside the write loop leaves the
session's console diverted into a half-written gz file; even on success the
gz connection is only reclaimed by the garbage collector ("closing unused
connection" warnings in long sessions). Success path verified clean;
failure path established from source.
Proposed change: open the connection into a variable,
`on.exit({sink(); close(con)}, add = TRUE)`.

**F3 [LOW, confidence: high] — record count reports individuals, but rows are loci (VRB1, DOC5 proposed)**
`R/gl2treemix.r:91-95` — the `verbose > 2` summary prints "Records written
... : nInd(x)".
Failure scenario: verified — a 3-individual, 5-locus object writes 6 lines
and reports "3".
Proposed change: report `nLoc(x)` (and the population count).

**F4 [LOW, confidence: high] — NULL return is visible (FS10, VRB5)**
`R/gl2treemix.r:104` — `return(NULL)` prints `NULL` on every unassigned
call, including at `verbose = 0` (verified via `withVisible`).
Proposed change: `invisible(NULL)`.

**F5 [LOW, confidence: high] — roxygen drift (DOC1; DOC2, DOC7 proposed rules)**
`R/gl2treemix.r:17-22` — the `verbose` `@param` is split by a stray blank
line (the default clause renders as a separate paragraph) and uses the
outdated wording (DOC2); the `@author` block names only a custodian
(DOC7); no `@details` tag (DOC1).
Proposed change: rejoin and modernise the `verbose` text, add
`Author(s):`, add a minimal `@details`; run `devtools::document()` (DOC4).

Note (not a finding): `gl.allele.freq(x, percent = TRUE, ...)` requests
percentages, but the `frequency` column is immediately discarded — the
argument is dead and can be dropped whenever the file is next touched.

## Proposed changes

1. Restrict the datatype gate with `accept = "SNP"`, or implement and
   document an explicit haploid branch (F1).
   **Consequence: SilicoDArT callers now error (or get halved counts)
   instead of silently doubled allele copies.**
2. Guard the sink and close the gz connection deterministically (F2).
3. Correct the verbose > 2 record count to loci (F3).
4. Replace `return(NULL)` with `invisible(NULL)` (F4).
5. Roxygen conformance pass + `devtools::document()` (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: full-matrix cell-by-cell recomputation of ref/alt counts on
  testset.gl (8 x 15, 7 populations, incl. an all-NA popxloc cell); gz
  integrity via `gzfile`; header/column order; SilicoDArT admission — run
  empirically
- Feeding the file to treemix itself: SKIPPED — no treemix binary on the
  review machine; format checked against the manual's specification
- FBM path (DAT6): SKIPPED — no FBM fixture for converters

## Approval (Phase B)

Blanket approval by Arthur Georges, 2026-09-05, via the formal approval
boxes: all findings at every severity approved, explicitly acknowledging
that converted outputs change where they were wrong; the DAT7 class fix
(fatal `accept = "SNP"` gate wherever SilicoDArT is silently admitted) is
approved campaign-wide and applies here as change 1 (the fatal gate, not
the haploid-branch alternative).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05; DAT7 fatal gate chosen, SilicoDArT callers now error |
| 2 | Approved | Arthur Georges | 2026-09-05 |
| 3 | Approved | Arthur Georges | 2026-09-05 |
| 4 | Approved | Arthur Georges | 2026-09-05 |
| 5 | Approved | Arthur Georges | 2026-09-05 |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2treemix` (base upstream/dev,
ddaed27). All five changes applied:

1. F1 — `utils.check.datatype(x, accept = "SNP", ...)`: SilicoDArT is a
   fatal error instead of silently doubled allele copies.
2. F2 — gz connection opened into a variable; depth-recorded `on.exit`
   unwinds the sink and closes the connection (tryCatch-guarded) on
   failure; the success path closes both explicitly.
3. F3 — `verbose > 2` summary reports loci and populations, not
   individuals.
4. F4 — `invisible(NULL)`.
5. F5 — roxygen: `verbose` param rejoined with the DOC2 default clause,
   DOC7 author labels, minimal `@details` added, house tag order;
   `devtools::document()` run (gl2treemix.Rd only).

Verification: baseline characterization test passes (128 assertions; the
count drops from 130 because the SilicoDArT admission block collapsed to
a single rejection check under F1, and the visibility flip maps to F4).
The 125-assertion cell-by-cell ref/alt recomputation on testset.gl is
unchanged and passes — SNP output identical. End-to-end run on
testset.gl at `verbose = 3` clean: 756 lines = header + 755 loci for 31
populations; no sink or connection left open. Sibling grep across the 8
clones: no callers of gl2treemix — all clear. NEWS.md entry added.

```json
{
  "function": "gl2treemix",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "medium", "rule": "none:sink-safety", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 5}
  ],
  "coverage_skipped": ["treemix run test: no treemix binary", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": null
}
```
