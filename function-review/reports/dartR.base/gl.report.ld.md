# Review: gl.report.ld (dartR.base)
- Family mode: report
- Date: 2026-09-10
- Reviewer: Claude (Claude Fable 5), dartr-function-review v1.0.0
- Package commit: f5e7b72 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.report.ld.r`)
- Datasets: platypus.gl (25- and 15-locus subsets), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.report.ld.R (new file, snapshot captured pre-review; 13 assertions, all passing)

## Verdicts

**Standards: Needs work** — dependency guards use `cat` + `return(-1)`
instead of `stop`, the internal `gl2gi` conversion prints at `verbose = 0`,
chunk files are written regardless of the `save` switch, and the roxygen
block is out of date on several points.

**Spec: Rework** — the title, description and return documentation promise
"pairwise population based LD ... between subpopulations", but the function
pools all individuals and never separates populations; the documented genind
input is rejected at the datatype check; and the crash-restart feature —
the function's stated reason for existing in this form — itself crashes at
`verbose < 2`. The per-population promise cannot be met by the finding list
alone: either the computation or the contract has to change.

## Independent verification (spec axis)

On a 25-locus platypus subset (81 individuals, 3 populations) the function
returns one pooled data frame of 300 pairs with `n = 81` — all individuals,
no population column, `seppop` absent from the code. The LD estimates
themselves are sound: spot checks of `R2` against `snpStats::ld`
R-squared on the same pooled data agree to ~1% (different maximum-likelihood
procedures; e.g. 0.0389 vs 0.0394), and `D`, `Dprime`, `r`, `X2`, `p` are
mutually consistent.

## Findings

**F1 [HIGH, confidence: high] — pools all populations despite the documented per-population purpose (DOC5 (proposed rule))**
`R/gl.report.ld.r:2-3,36-41,95-98` — the title says "pairwise population
based Linkage Disequilibrium", the return doc says "pairwise LD across all
loci between subpopulations"; the code converts the whole object with
`gl2gi(x)` and computes LD over all individuals pooled. Pooling across
differentiated populations inflates LD (the two-locus Wahlund effect), so
the numbers are not the within-population LD the documentation describes.
Failure scenario: any structured dataset — the returned statistics mix
between-population allele-frequency differences into every pair; confirmed
`n = 81` (all three platypus populations) on every row.
Proposed change: decide the contract. Option A (bounded, docs-only): retitle
and describe as whole-dataset LD, pointing users to `gl.report.ld.map` for
within-population LD. Option B (behaviour): loop over `seppop(x)` and add a
`pop` column. The custodian owns this choice; Option A is proposed here.

**F2 [HIGH, confidence: high] — chunk-restart crashes below verbose 2, and looks for chunk files in the wrong directory (DOC5 (proposed rule), VRB1)**
`R/gl.report.ld.r:186-196` — when all pairs are already done, the early
`return(lddone)` sits *inside* `if (verbose >= 2)`; at `verbose 0/1` the
function falls through, `allp <- allp[,-c(1:done)]` empties the pair matrix
and the parallel loop dies. Additionally `:159-160` discovers chunk files
with `list.files(pattern = ...)` in the *working directory* while `:170-173`
loads them from `outpath` — with the default `outpath = tempdir()` a rerun
from a normal working directory never finds its chunks.
Failure scenario: rerunning with the same `chunkname` at `verbose = 0`
errors "task 1 failed - subscript out of bounds" (confirmed); the identical
call at `verbose = 2` returns the cached 300 pairs. With `outpath` left at
its default the restart feature silently recomputes everything.
Proposed change: return unconditionally (message still gated); pass
`path = outpath` to `list.files`.

**F3 [MEDIUM, confidence: high] — genind input is documented but rejected (DOC5 (proposed rule))**
`R/gl.report.ld.r:15-16,92-98` — `@param x` accepts "a genlight or genind
object"; `utils.check.datatype(x)` stops on a genind ("found genind
expecting genlight...", confirmed), so the `is(x, "genlight")` branch that
was meant to convert only genlight objects never sees a genind.
Failure scenario: `gl.report.ld(gl2gi(x))` errors despite the docs.
Proposed change: docs-only — accept genlight only (the internal conversion
sentence also needs rewording).

**F4 [MEDIUM, confidence: high] — dependency guards use cat + return(-1) (DEP1)**
`R/gl.report.ld.r:53-80` — three guards print with `cat(error(...))` and
`return(-1)` instead of `stop(error(...))`; a scripted caller receives `-1`
as if it were a result. `data.table` and `foreach` are in Imports, so two of
the three guards are dead code; only `doParallel` (Suggests) needs a guard.
Failure scenario: without doParallel, `res <- gl.report.ld(x)` puts `-1` in
`res` and the script continues.
Proposed change: single DEP1 guard for `doParallel` using `stop(error(...))`.

**F5 [MEDIUM, confidence: high] — internal gl2gi call prints at verbose = 0 (VRB5)**
`R/gl.report.ld.r:96-98` — `gl2gi(x)` is called without `verbose = 0`, so
its four progress lines print at the caller's `verbose = 0` (confirmed).
Proposed change: `gl2gi(x, verbose = 0)` (or pass the caller's verbosity).

**F6 [MEDIUM, confidence: medium] — chunk files are written regardless of save = FALSE, with collision-prone names (DOC5 (proposed rule))**
`R/gl.report.ld.r:374-379` — every chunk is saved to
`LD_chunks_<chunkname>_<i>.rdata` unconditionally; with the default
`chunkname = NULL` the name is `LD_chunks__i.rdata`, identical across runs
and datasets in the same `outpath`.
Failure scenario: `save = FALSE` still writes files (confirmed); two
analyses sharing an `outpath` overwrite each other's restart state, and a
later restart resumes from the wrong dataset's chunks.
Proposed change: write chunks only when `save = TRUE`; derive a default
chunkname from the data when saving.

**F7 [LOW, confidence: high] — SilicoDArT admitted (DAT7 (proposed rule))**
`R/gl.report.ld.r:92-93` — the default `accept` admits presence/absence
data; the pooled dosage math then emits LD-like numbers from 0/1 calls
without warning (confirmed: runs to completion on testset.gs).
Proposed change: `accept = "SNP"`.

**F8 [LOW, confidence: high] — roxygen block out of date (DOC1, DOC2, DOC7 (proposed rule))**
`R/gl.report.ld.r:4` — the `@family` tag is indented behind the title
continuation, so roxygen absorbs it into the title: `man/gl.report.ld.Rd`'s
`\title{}` literally ends with "@family matched report" and the function has
no family cross-links (confirmed in the generated Rd). `:26-27` documents
`probar` "[default = TRUE]" while the signature default is FALSE; the
`verbose` text deviates from the DOC2 standard; `@author` names an author
with no Custodian line; progress text says "No. of Simulations" for what are
exact calculations, and "Loooking".
Proposed change: fix the tag indentation (and regenerate the Rd), correct
the `probar` default, standard verbose text, `Author(s)/Custodian` per DOC7,
message wording.

**F9 [INFO, confidence: high] — LD.fast estimates differ from snpStats by ~1%**
`R/gl.report.ld.r:228-314` — the inline maximum-likelihood estimator
(`optimize` over the log-likelihood) and `snpStats::ld` (EM) agree to about
two significant figures on spot checks (0.0389 vs 0.0394). No action;
recorded because `gl.report.ld.map` reports the snpStats numbers, so the two
reports give slightly different values for the same pair by design.

## Proposed changes

1. Align the contract with the computation: retitle/redescribe as
   whole-dataset (pooled) LD and cross-reference `gl.report.ld.map` for
   within-population LD (F1, Option A; Option B — a `seppop` loop with a
   `pop` column — is the custodian's alternative and is NOT drafted here).
2. Fix the restart path: unconditional early return; discover chunk files in
   `outpath`; correct the `length(chunkfiles > 0)` guard (F2).
3. Replace the three cat/return(-1) guards with one `stop(error(...))` guard
   for doParallel (F4).
4. Silence the internal `gl2gi` call (F5).
5. Write chunk files only when `save = TRUE` (F6).
6. Restrict the datatype check to SNP data (F7).
7. Roxygen repairs: family-tag indentation + Rd regeneration, probar
   default, DOC2 verbose text, DOC7 author/custodian, message wording (F3,
   F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: pooled-vs-per-population claim — run; restart path — run (crash
  reproduced and cached return verified); genind claim — run; SilicoDArT —
  run; estimator spot-check vs snpStats — run
- ncores > 1 path: SKIPPED — logic identical to ncores = 1 but split across
  workers; not exercised to keep the test suite fast
- FBM path (DAT6): SKIPPED — `gl2gi` densifies by design; no FBM fixture
- Google Group search: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Arthur | Option A (docs-only): describe pooled behaviour, cross-reference gl.report.ld.map |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |

Cross-package caller grep (API3): no calls to `gl.report.ld` in the local
dartR.* clones outside dartR.base. All clear.

## Outcome

Changes 1-7 applied on branch review-gl.report.ld (commit e16f9bb), PR
green-striped-gecko/dartR.base#392. Change 1 applied as the approved
docs-only Option A; Option B (per-population computation) remains recorded
here as the custodian's alternative.

- Characterization suite green (13 assertions); every diff from the
  pre-review baseline maps to an approved change: silence at verbose 0
  (F5), restart returning cached results at any verbosity and with wd !=
  outpath (F2), no chunk files with save = FALSE (F6), SilicoDArT rejected
  (F7).
- Returned statistics unchanged: pooled 25-locus platypus baseline (300
  pairs; first-pair D/R2/n to 1e-6).
- End-to-end run at verbose = 3 clean; chunk/result files named
  LD_chunks_LDallp_*.rdata / LDallp.rdata.
- API notes recorded in NEWS and the PR body: doParallel guard now stops,
  save = FALSE writes nothing, default chunk-file names changed.

```json
{
  "function": "gl.report.ld",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "1.0.0",
  "commit": "f5e7b72",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "approved", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DAT7", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "no-action", "change": null}
  ],
  "coverage_skipped": ["ncores>1 not exercised", "DAT6: gl2gi densifies by design", "Google Group: not queried"],
  "status": "pr-open",
  "pr": 392
}
```
