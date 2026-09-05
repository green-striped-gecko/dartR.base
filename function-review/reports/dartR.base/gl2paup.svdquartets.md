# Review: gl2paup.svdquartets (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl (SNP), testset.gs (SilicoDArT)
- Baseline: tests/testthat/test-gl2paup.svdquartets.R (snapshot captured pre-review; 22 assertions, all passing against current code)

## Verdict

**Standards: Needs work** — the SNP path is sound (the ambiguity lookup table is verified correct: A/T→W, A/G→R, A/C→M, T/G→K, T/C→Y, G/C→S, and method 1 splits ref/alt lines consistently), but `verbose = 0` leaks internal `gl.sort` output, and the monomorph check breaks on plain genlight objects.
**Spec: Needs work** — the SilicoDArT branch, which the roxygen explicitly supports, writes 0/1 characters under a `format datatype = dna` declaration that PAUP will reject; several PAUP-block filenames ignore `outfile`.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT output declared as DNA (DOC5)**
`R/gl2paup.svdquartets.r:276` — the single `format datatype = dna gap = - ;` line serves both branches. The ploidy-1 branch (lines 226-245) writes 0/1/? presence-absence strings, so the resulting NEXUS declares DNA data whose matrix contains `0` and `1` — symbols outside the DNA alphabet. Verified in the written file.
Failure scenario: the documented tag P/A workflow ("If the data are tag presence/absence, then method=2 is assumed") produces a file PAUP cannot execute; the parsimony converter next door writes `datatype = standard` for the same data.
Proposed change: emit `format datatype = standard gap = - ;` (with an appropriate symbols list) when `datatype == "SilicoDArT"`.

**F2 [MEDIUM, confidence: high] — `verbose = 0` leaks internal output (VRB3, VRB5)**
`R/gl2paup.svdquartets.r:81,131` — the monomorph warning is ungated, and `gl.sort(x, sort.by = "pop")` is called without `verbose = 0`, so the nested function prints its full start/warning/end banner. Verified: six lines captured at `verbose = 0`, including "Starting gl.sort".
Failure scenario: silent pipelines emit multi-line noise; the leak scales with global verbosity settings.
Proposed change: pass `verbose = 0` to `gl.sort` and gate the warning per VRB3. Note the SilicoDArT branch never sorts — the SNP branch sorts by population while the ploidy-1 branch writes individuals in input order, so taxpartitions for unsorted SilicoDArT input can be wrong (see F4 companion note).

**F3 [MEDIUM, confidence: high] — plain genlight objects fail opaquely (DAT5)**
`R/gl2paup.svdquartets.r:80` — `if (!x@other$loc.metrics.flags$monomorphs)` assumes the dartR bookkeeping exists. A genlight built outside dartR dies with "invalid argument type". Verified.
Failure scenario: users converting an externally built genlight (the very audience of a format converter) get an uninterpretable error instead of the DAT5-mandated clear message or a compliance-check pass.
Proposed change: test for the flag's existence first (the pattern already applied to `gl2fasta` in this campaign), or fall back to `gl.filter.monomorphs` counting as `gl2paup.parsimony` does.

**F4 [MEDIUM, confidence: high] — single-population input crashes; SilicoDArT branch skips sorting (DAT5)**
`R/gl2paup.svdquartets.r:251-258` — the taxpartition loop `for (i in 2:length(b))` indexes `b[0]` when only one population exists: "replacement has length zero". Verified. Additionally, the ploidy-1 branch (line 226) uses `pop(x)` in input order without the `gl.sort` applied to SNP data, so a SilicoDArT object whose individuals are not already grouped by population gets taxpartition ranges that do not correspond to the matrix rows.
Failure scenario: one-population export dies; unsorted SilicoDArT input silently mislabels taxa blocks in PAUP.
Proposed change: guard the loop for one population; sort once, before the branch split.

**F5 [MEDIUM, confidence: high] — `method = 0` silently accepted (DOC5)**
`R/gl2paup.svdquartets.r:88` — the guard `method < 0 | method > 2` admits 0, which then behaves as method 2 with no warning. Verified: no "method must be" warning at `method = 0`.
Failure scenario: a typo (`method = 0`) silently produces the one-line format when the two-line format was intended.
Proposed change: `!method %in% c(1, 2)`.

**F6 [LOW, confidence: high] — PAUP block filenames ignore `outfile` (DOC5)**
`R/gl2paup.svdquartets.r:298,307,308` — `log file=svd.txt`, `treeFile=svd.tre` and `savetrees file=svd_boot.tre` are hardcoded. Two conversions run from the same directory overwrite each other's logs and trees even with distinct `outfile` values. Verified in the written file.
Proposed change: derive these names from `outfile` as `gl2paup.parsimony` does from `outfileprefix`.

**F7 [LOW, confidence: high] — taxon-name sanitisation narrower than the sibling converter (DOC5)**
`R/gl2paup.svdquartets.r:98-101` — populations get spaces and parentheses replaced; individual names only get spaces replaced. Parentheses in individual names reach the NEXUS matrix and sets block, where PAUP treats them as punctuation.
Proposed change: apply the parsimony converter's full substitution set to `indNames` too.

**F8 [LOW, confidence: high] — visible NULL return (FS10)**
`R/gl2paup.svdquartets.r:327` — `return(NULL)` prints `NULL` on every call. Proposed change: `invisible(NULL)`.

**F9 [LOW, confidence: high] — roxygen order and wording drift (DOC1, DOC2, DOC7 (proposed rule))**
`R/gl2paup.svdquartets.r:1-44` — `@return` after `@export` (pre-ratification order; per DOC1 flag the file); `verbose` default clause not the ratified DOC2 wording; author block names only a custodian (DOC7 — default fix: add `Author(s): Arthur Georges.`). `nbootstraps` is never validated (contrast the parsimony converter's advisory warnings).

## Proposed changes

1. Branch the `format datatype` line on `datatype`, writing `standard` for SilicoDArT (F1).
2. Call `gl.sort` with `verbose = 0`, move it ahead of the ploidy branch, and gate the monomorph warning (F2, F4 sorting half).
3. Make the monomorph check tolerant of missing `loc.metrics.flags` (F3).
4. Guard taxpartition construction for a single population (F4).
5. Validate `method` with `%in% c(1, 2)` (F5).
6. Derive PAUP log/tree filenames from `outfile` (F6).
7. Extend individual-name sanitisation to parentheses (F7).
8. Return `invisible(NULL)` (F8).
9. Reorder/reword roxygen per DOC1/DOC2/DOC7; run `devtools::document()` (F9, DOC4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plots), STY — run
- Spec: methods 1 and 2 executed on testset.gl subset; NEXUS inspected (dimensions, format, matrix records, ambiguity codes verified base-by-base against `loc.all`, taxpartition); SilicoDArT branch executed on testset.gs subset — run
- Ambiguity lookup table: verified against all six heterozygote pairs — correct
- One-population, `method = 0`, plain-genlight, `verbose = 0` probes: run
- Generated PAUP command block executed in PAUP itself: SKIPPED — no PAUP binary on this machine; syntax reviewed statically
- FBM path (DAT6): SKIPPED — no FBM fixture; converter densifies via `as.matrix`/`data.frame` by design

## Approval (Phase B)

All findings at every severity approved by Arthur Georges, 2026-09-05, via
the formal approval boxes -- a blanket class approval covering this batch,
explicitly acknowledging that converted outputs change where they were
wrong. DAT7 policy note: this function documents SilicoDArT support, so
the approved fix is the report's proposed one (emit datatype = standard
for the silico branch), not a fatal SNP-only gate.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 2 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 3 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 4 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 5 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 6 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 7 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 8 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 9 | Approved | Arthur Georges | 2026-09-05, blanket class approval |

## Outcome (Phase C)

All nine changes applied on branch `review-gl2paup.svdquartets`
(base `upstream/dev` = ddaed27), 2026-09-05. Two notes beyond the letter
of the findings, both recorded here for the custodian:

- Implementing F5's method validation surfaced that the documented
  assumption "If the data are tag presence/absence, then method=2 is
  assumed" was never enforced -- SilicoDArT with `method = 1` crashed on
  the unbuilt `refseq` ("object 'refseq' not found", pre-existing, loud).
  The documented assumption is now implemented: silico input coerces
  method to 2 with a warning at `verbose >= 1`. Classified MEDIUM (loud
  crash on an explicit-argument path, no wrong output), folded into
  change 5.
- F3 is applied as scoped (the monomorph flag check tolerates flag-less
  objects), but a fully plain genlight still fails downstream inside
  `gl.sort`, whose own `ind.metrics` subsetting is not flag-tolerant
  ("invalid subscript type 'list'"). That is `gl.sort`'s defect, outside
  this function's scope; the baseline test documents the relocated
  failure.

Verification: baseline characterization test updated -- every diff maps
to an approved finding (F1 standard format, F2/F8 silence at verbose = 0,
F3 relocated failure, F4 single-population, F5 method coercion, F6
outfile-derived PAUP filenames) -- 22 assertions pass. End-to-end runs at
verbose = 3 on testset.gl (methods 1 and 2) and testset.gs completed
cleanly. Sibling grep across the 8 clones: no package-code callers of
`gl2paup.svdquartets` -- all clear.
PR: https://github.com/green-striped-gecko/dartR.base/pull/PRNUM

```json
{
  "function": "gl2paup.svdquartets",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9}
  ],
  "coverage_skipped": ["PAUP execution of generated command block: no PAUP binary", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": null
}
```
