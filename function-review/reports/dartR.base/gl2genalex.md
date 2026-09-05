# Review: gl2genalex (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl2genalex.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the poppr wrapper is sound, but the dependency
guard returns `-1` instead of stopping, SilicoDArT slips through, and the
function is not silent at `verbose = 0`.
**Spec: Ready** — the written csv matches the GenAlEx codominant layout
(header counts, pop-name row, Ind/Pop/locus columns, A=1 C=2 G=3 T=4 allele
codes, 0 for missing), verified cell-level against the input dosages.

## Findings

**F1 [HIGH, confidence: high] — SilicoDArT admitted and exported as codominant genotypes (DAT7)**
`R/gl2genalex.r:59` — the datatype check keeps the permissive default. A
presence/absence object passes through `gl2gi()` and poppr, producing a
well-formed GenAlEx codominant file whose "genotypes" are artefacts of
coercing 0/1 ploidy-1 data (confirmed on `testset.gs`: plausible-looking
allele pairs per locus).
Failure scenario: a user exports SilicoDArT data and analyses meaningless
codominant genotypes in GenAlEx with no hint anything is wrong.
Proposed change: `accept = "SNP"` in `utils.check.datatype()`.

**F2 [MEDIUM, confidence: high] — dependency guard returns -1 instead of stopping (DEP1)**
`R/gl2genalex.r:62-70` — when poppr is absent the function prints the error
text with `cat(error(...))` and `return(-1)`.
Failure scenario: in a script, the call "succeeds" with value `-1`, no file
is written, no condition is raised, and the pipeline continues; `tryCatch`
cannot intercept it.
Proposed change: use the DEP1 idiom —
`stop(error("Package poppr needed for this function to work. Please install it.\n"))`.

**F3 [MEDIUM, confidence: high] — poppr chatter at verbose = 0 (VRB5)**
`R/gl2genalex.r:78-83` — `poppr::genind2genalex()` prints "Extracting the
table ... Writing the table to <path> ... Done." unconditionally (confirmed
at `verbose = 0`).
Failure scenario: silent pipelines still emit per-call progress text.
Proposed change: pass `quiet = (verbose < 2)` through to
`genind2genalex()`.

**F4 [LOW, confidence: high] — all-NA loci silently dropped (VRB4) (proposed rule)**
`R/gl2genalex.r:75` — `gl.filter.allna(x, verbose = 0)` can remove loci
before export with no message at any verbosity, so the file's locus count
can differ from the input object's.
Failure scenario: a user cross-referencing locus counts between object and
file sees an unexplained discrepancy.
Proposed change: report the number of dropped loci at `verbose >= 1` when it
is non-zero.

**F5 [LOW, confidence: high] — visible NULL return prints at verbose = 0 (FS10, VRB5)**
`R/gl2genalex.r:97` — `return(NULL)` is visible; an unassigned call prints
`NULL`.
Failure scenario: console noise at `verbose = 0`.
Proposed change: `return(invisible(NULL))`.

**F6 [LOW, confidence: medium] — overwrite semantics documented as skip, implemented as error (DOC5) (proposed rule)**
`R/gl2genalex.r:21-22` — the `overwrite` doc says "the file will not be
overwritten"; poppr actually raises an error when the file exists
(confirmed; baseline test captures it). The existing file does survive, but
the call fails rather than quietly declining.
Failure scenario: a scripted export with `overwrite = FALSE` against an
existing file aborts the script instead of skipping.
Proposed change: document the error, or catch it and warn.

## Proposed changes

1. `accept = "SNP"` (F1).
2. Replace the guard with the DEP1 `stop(error(...))` idiom (F2).
   **Consequence: callers relying on the -1 return see an error instead.**
3. Pass `quiet = (verbose < 2)` to `genind2genalex()` (F3).
4. Report non-zero all-NA locus drops at `verbose >= 1` (F4).
5. Return `invisible(NULL)` (F5).
6. Align `overwrite` documentation with the raise-on-exists behaviour (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a), STY — run
- Spec: written csv vs GenAlEx codominant layout on testset.gl — run
  (header rows, pop sizes/names, allele codes, missing code verified
  cell-level in the baseline test)
- SilicoDArT probe — run (F1)
- verbose = 0 silence probe — run (F3, F5)
- overwrite = FALSE probe — run (F6)
- poppr-absent path (F2): reasoned from code only — poppr is installed here,
  so the -1 return was not executed
- FBM path (DAT6): SKIPPED — no FBM fixture exercised for this function
- Genome-position lens: N/A — the format carries no coordinate fields

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 2 | Approved | Arthur Georges (2026-09-05) | Callers relying on -1 return now see an error |
| 3 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 4 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 5 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |
| 6 | Approved | Arthur Georges (2026-09-05) | Blanket approval, all severities |

All findings approved 2026-09-05 via the formal approval boxes: blanket
class approval at every severity, explicitly including the consequence
that converted outputs change where they were wrong, and the DAT7 fatal
`accept = "SNP"` gate wherever SilicoDArT was silently admitted.

## Outcome (Phase C)

All six changes applied on branch `review-gl2genalex` (base
`upstream/dev`, ddaed27): `accept = "SNP"` (F1), DEP1 `stop(error(...))`
guard replacing the -1 return (F2), `quiet = (verbose < 2)` passed to
`poppr::genind2genalex()` (F3), non-zero all-NA locus drops reported at
`verbose >= 1` (F4), `return(invisible(NULL))` (F5), and the `overwrite`
doc aligned with the raise-on-exists behaviour (F6). Baseline
characterization test updated for two approved diffs (F1 rejection;
F3/F5 silence and invisibility); 14/14 assertions pass — the SNP-path
csv snapshot (header counts, pop names, allele codes, first data row)
is byte-identical. End-to-end run on `testset.gl` at `verbose = 3`
wrote 277 lines (3 header + 274 individuals), with the F4 drop warning
reporting this checkout's 3 all-NA loci. The poppr-absent guard path
remains reasoned from code only (poppr is installed here). Sibling-caller
grep across the 8 dartRverse clones: no callers of `gl2genalex` — all
clear. NEWS entry added.

```json
{
  "function": "gl2genalex",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "VRB4", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 6}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "poppr-absent guard path not executed", "genome-position lens: format has no coordinate fields"],
  "status": "phase-c",
  "pr": null
}
```
