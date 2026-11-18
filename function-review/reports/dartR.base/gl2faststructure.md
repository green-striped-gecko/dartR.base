# Review: gl2faststructure (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2faststructure.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the preamble conforms, but the permissive datatype
default admits presence/absence data, the roxygen block has no `@examples`,
and the visible `NULL` return leaks output at `verbose = 0`.
**Spec: Ready** — on SNP data the file matches the fastStructure str layout
(two rows per individual, six leading columns, 1/2 allele codes, `-9`
missing); recoding verified locus-by-locus against the 0/1/2 dosage.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted and written as pseudo-diploid (DAT7)**
`R/gl2faststructure.r:49` — `utils.check.datatype(x, verbose = verbose)` uses
the permissive default, so presence/absence data passes. The converter then
duplicates each 0/1 row and recodes presence (1) as a heterozygote pair
`1 2`, emitting a well-formed diploid SNP file from data that has no diploid
meaning.
Failure scenario: `gl2faststructure(testset.gs)` runs without error and
writes a file fastStructure will happily analyse; results are genetically
meaningless. Confirmed empirically.
Proposed change: pass `accept = "SNP"` to `utils.check.datatype()`.

**F2 [MEDIUM, confidence: high] — roxygen block has no @examples or @details (DOC1)**
`R/gl2faststructure.r:1-27` — the header lacks `@examples` entirely (DOC1
requires it; every sibling converter reviewed carries one) and has no
`@details`.
Failure scenario: no runnable example is checked by CRAN or available from
`?gl2faststructure`; format-specific caveats (six ignored columns, `-9`
missing code) have no home.
Proposed change: add a `\donttest{}` example on `testset.gl` with
`outpath = tempdir()`, and a short `@details` describing the emitted layout.

**F3 [LOW, confidence: high] — visible NULL return prints at verbose = 0 (FS10, VRB5)**
`R/gl2faststructure.r:102` — `return(NULL)` is visible, so an unassigned call
prints `NULL` to the console even at `verbose = 0`.
Failure scenario: `gl2faststructure(x, verbose = 0)` at the prompt emits
`NULL`, violating the fully-silent contract.
Proposed change: `return(invisible(NULL))`.

**F4 [LOW, confidence: high] — unused file connection; file reopened per individual (STY3)**
`R/gl2faststructure.r:61-79` — `zz <- file(outfilespec, "w")` is opened but
never written to; its only effect is truncating the file. Each `write(...,
file = outfilespec, append = TRUE)` inside the loop reopens the file by
path, once per individual.
Failure scenario: no wrong output, but O(nInd) file opens slow large exports
and the dangling connection idiom obscures intent.
Proposed change: write through the open connection (`file = zz`) and drop
`append`, or drop the connection and truncate explicitly before the loop.

**F5 [LOW, confidence: high] — verbose text and author block off-canon (DOC2, DOC7) (proposed rule)**
`R/gl2faststructure.r:19-23` — the `verbose` param text ends
"[default 2 or as specified using gl.set.verbosity]" instead of the canonical
default clause (DOC2), and `@author` names Bernd Gruber without the
`Author(s):` / `Custodian:` labels (DOC7, proposed).
Failure scenario: doc drift across the family; custodianship unstated.
Proposed change: adopt the DOC2 canonical clause; relabel as
`Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to ...`.

Note (no finding): the six leading columns hold the individual's row index
rather than `indNames(x)`. fastStructure ignores these columns, so output is
usable, but carrying real individual labels in column 1 would make the file
self-documenting. Worth folding into any touch-up of the writer.

## Proposed changes

1. Restrict datatype with `accept = "SNP"` (F1).
2. Add `@examples` (\donttest, tempdir) and `@details`; run
   `devtools::document()` (F2).
3. Return `invisible(NULL)` (F3).
4. Tidy the writer: single open connection, optionally emit `indNames(x)` in
   column 1 (F4, note).
5. Canonical DOC2 verbose text and DOC7 author labels (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: behaviour vs roxygen and fastStructure str layout on testset.gl — run
  (locus-by-locus recoding check in baseline test)
- SilicoDArT admission probe on testset.gs — run
- verbose = 0 silence probe — run
- Edge cases: single individual (works), all-NA locus (coded -9) — run
- FBM path (DAT6): SKIPPED — no FBM fixture exercised for this function
- Genome-position lens: N/A — the format carries no coordinate fields

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 2 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 3 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 4 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 5 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |

All findings approved 2026-09-05 via the formal approval boxes: blanket
class approval at every severity, explicitly including the consequence
that converted outputs change where they were wrong, and the DAT7 fatal
`accept = "SNP"` gate wherever SilicoDArT was silently admitted.

## Outcome (Phase C)

All five changes applied on branch `review-gl2faststructure` (base
`upstream/dev`, ddaed27). Baseline characterization test updated for two
approved diffs: SilicoDArT rejection (F1) and invisible NULL return (F3);
17/17 assertions pass. End-to-end run on `testset.gl` at `verbose = 3`
wrote 548 rows (2 x 274 individuals) to tempdir. Sibling-caller grep
across the 8 dartRverse clones: no code callers of `gl2faststructure`
(one doc mention in `gl2structure.r`) — all clear. NEWS entry added.
PR: green-striped-gecko/dartR.base#333.

```json
{
  "function": "gl2faststructure",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "applied", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "genome-position lens: format has no coordinate fields"],
  "status": "phase-c",
  "pr": 333
}
```
