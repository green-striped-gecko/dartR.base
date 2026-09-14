# Review: gl2paup.parsimony (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gs (SilicoDArT), testset.gl (rejection check)
- Baseline: tests/testthat/test-gl2paup.parsimony.R (snapshot captured pre-review; 21 assertions, all passing against current code)

## Verdict

**Standards: Needs work** — the FS skeleton is present and datatype restriction to SilicoDArT is correct, but warnings and one fatal error bypass verbosity gating, and `verbose = 0` is not silent.
**Spec: Needs work** — the core NEXUS output is valid (`datatype = standard`, `?` for missing, correct taxpartition on the reference data), but the bash-mode subsystem writes scripts to the wrong directory, embeds another user's hardcoded path, and a failed precondition check does not stop the run.

## Findings

**F1 [HIGH, confidence: high] — "Fatal Error" that does not stop (FS5, VRB2)**
`R/gl2paup.parsimony.r:158-162` — the bootstraps/ncpus divisibility check prints `cat(error("Fatal Error: ..."))` with no `stop()`. Execution continues.
Failure scenario: `out.type = "bash"`, `nbootstraps = 10`, `ncpus = 3` prints the fatal message, then writes `bootstrap nreps=3.33333333333333` into the NEXUS PAUP block and generates all three shell scripts. PAUP receives a fractional replicate count. Verified empirically; snapshotted in the baseline test.
Proposed change: replace the `cat(error(...))` with `stop(error(...))` (the neighbouring `base.dir.name` check at lines 164-169 already stops, though via a bare `stop()` after the cat — normalise both to `stop(error(...))`).

**F2 [HIGH, confidence: high] — bash-mode scripts written to the working directory, not outpath (FS7)**
`R/gl2paup.parsimony.r:308,314,369,372,424,427` — `outfilespec2/3/4` are bare filenames (`paste0("generator_", ...)`) opened relative to the current working directory, while the NEXUS file goes to `outpath`. Verified: with `outpath = tempdir()` and a different cwd, the three `generator_*.sh` files land in the cwd.
Failure scenario: any user with a default `outpath` (tempdir per CRAN policy) finds the companion scripts scattered in whatever directory R happens to be running in — or, for a package check, polluting the check directory (a CRAN policy violation).
Proposed change: build all three with `file.path(outpath, ...)`.

**F3 [HIGH, confidence: high] — maketrees script hardcodes a personal directory and ignores base.dir.name (DOC5)**
`R/gl2paup.parsimony.r:399` — the generated PBS job contains `cd /g/data/xl04/ag3760/parsimony/` where the roxygen promises `base.dir.name` is "the base directory on the server to act as the workspace". The consensus script (line 443) does use `base.dir.name`; the maketrees jobs do not. `#PBS -P xl04` (lines 385, 430) and the PAUP binary path are likewise hardcoded.
Failure scenario: any user other than the original author submits bootstrap jobs that cd into another person's gdata directory and fail (or worse, run there).
Proposed change: substitute `base.dir.name` at line 399; consider surfacing the PBS project as a parameter or deriving it from `base.dir.name`.

**F4 [MEDIUM, confidence: high] — single-population input crashes (DAT5)**
`R/gl2paup.parsimony.r:252-255` — the taxpartition loop `for (i in 2:length(b))` runs `2:1` when only one population exists, and `a[i] <- b[i-1] + 1` indexes `b[0]`, giving "replacement has length zero". Verified empirically.
Failure scenario: a genlight with one population (common for a quick single-population export) dies with an uninformative subscript error after the expensive genotype-string construction.
Proposed change: guard the loop with `if (length(b) > 1)` and write the single-partition line directly (or skip the sets block for one population).

**F5 [MEDIUM, confidence: high] — consensus generator script is malformed (DOC5)**
`R/gl2paup.parsimony.r:427-437` — the consensus script begins directly with `#PBS -P xl04` (no `#!/bin/bash` shebang, unlike the other two scripts) and its `#PBS -N` line interpolates the literal `${i}` from a loop that does not exist in this script.
Failure scenario: `qsub generator_<prefix>_consensus.sh` submits a job whose name contains a literal `${i}`; depending on shell, the missing shebang changes interpretation of the here-doc-free body.
Proposed change: write the shebang first and give the consensus job a fixed `-N` name.

**F6 [MEDIUM, confidence: high] — the documented examples rely on a nonexistent `outfile` parameter (DOC5)**
`R/gl2paup.parsimony.r:81-82` — both `@examples` call `gl2paup.parsimony(gg, outfile = "test.nex", ...)`. The function has no `outfile`; R's partial matching binds it to `outfileprefix`, producing `test.nex_bootstrap.nex` (a doubled extension). Verified.
Failure scenario: users copy the example, pass `outfile`, and get files named `<name>.nex_bootstrap.nex`; if a future edit adds a genuine `outfile` parameter the examples silently change meaning.
Proposed change: update the examples to `outfileprefix = "test"`.

**F7 [MEDIUM, confidence: high] — `verbose = 0` is not silent (VRB3, VRB5)**
`R/gl2paup.parsimony.r:146,154,160,166,173,179,182` — the monomorph warning, all parameter-validation warnings, and both fatal-error cats are ungated. Verified: `capture.output(..., verbose = 0)` returns three warning lines plus a printed `NULL`.
Failure scenario: pipelines running quietly still emit console noise; conversely VRB4-class content (parameters coerced) cannot be distinguished from chatter.
Proposed change: gate advisory warnings with `if (verbose >= 2)` (or `>= 1` where results-affecting, VRB4); keep fatal errors as `stop(error(...))` which needs no gate.

**F8 [LOW, confidence: high] — visible NULL return (FS10)**
`R/gl2paup.parsimony.r:489` — `return(NULL)` prints `NULL` at the console on every interactive call. Proposed change: `invisible(NULL)` (as `gl2plink` already does).

**F9 [LOW, confidence: high] — roxygen order and wording drift (DOC1, DOC2, DOC7 (proposed rule))**
`R/gl2paup.parsimony.r:1-85` — `@return` sits after `@export` and `@details` after `@param` (the pre-2026-08-27 order; per DOC1, flag the file); the `verbose` default clause reads "[default 2 or as specified using gl.set.verbosity]" rather than the ratified NULL-cascade wording (DOC2); the author block names only a custodian with no `Author(s):` line (DOC7, proposed rule — default fix: add `Author(s): Arthur Georges.`).

**F10 [LOW, confidence: medium] — test-mode messaging and parameter handling (STY1, DOC5)**
`R/gl2paup.parsimony.r:127-140` — `test = TRUE` reports "... scored for N SNP loci" in a SilicoDArT-only function, silently overrides a user-supplied `ncpus` to 10, and passes `qtl` (a possibly fractional quantile) as `n` to `gl.subsample.ind`. Not exercised empirically (see Coverage).
Proposed change: correct the message; respect user `ncpus`; `floor(qtl)`.

## Proposed changes

1. Convert the divisibility check to `stop(error(...))`; normalise the `base.dir.name` check to the same idiom (F1).
2. Route `outfilespec2/3/4` through `file.path(outpath, ...)` (F2).
3. Substitute `base.dir.name` for the hardcoded `cd` path in the maketrees jobs; parameterise or document the PBS project (F3).
4. Guard the taxpartition construction for a single population (F4).
5. Add the shebang and fix the `-N` line in the consensus script (F5).
6. Update `@examples` to use `outfileprefix` (F6).
7. Gate the ungated warnings behind verbosity per VRB3/VRB4 (F7).
8. Return `invisible(NULL)` (F8).
9. Reorder/reword the roxygen block per DOC1/DOC2/DOC7 and run `devtools::document()` (F9, DOC4).
10. Fix test-mode message, `ncpus` override, and fractional `qtl` (F10).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plots), STY — run
- Spec: standard and bash out.types executed on testset.gs subset; NEXUS inspected (header, dimensions, format line, matrix records, taxpartition, PAUP block) — run
- SNP rejection (DAT7): run — correctly refused with `accept = "SilicoDArT"`
- One-population and verbose = 0 behaviour: run
- `test = TRUE` subsetting path: SKIPPED — stochastic subsampling touches `gl.subsample.ind`/`gl.keep.pop` at runtime; reviewed by reading only (F10)
- Generated bash/PBS scripts on an actual Gadi queue: SKIPPED — no server access; scripts reviewed statically
- FBM path (DAT6): SKIPPED — no FBM fixture; converter densifies via `as.matrix` by design

## Approval (Phase B)

All findings at every severity approved by Arthur Georges, 2026-09-05, via
the formal approval boxes — a blanket class approval covering this batch,
explicitly acknowledging that converted outputs change where they were
wrong, and ratifying the DAT7 fatal-gate policy (already satisfied here:
the function carries `accept = "SilicoDArT"`).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 2 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 3 | Approved | Arthur Georges | 2026-09-05; PBS project surfaced as `pbs.project` parameter |
| 4 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 5 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 6 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 7 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 8 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 9 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 10 | Approved | Arthur Georges | 2026-09-05, blanket class approval |

## Outcome (Phase C)

All ten changes applied on branch `review-gl2paup.parsimony`
(base `upstream/dev` = ddaed27), 2026-09-05. The PBS project is surfaced as
a new `pbs.project` parameter [default 'xl04']; the storage directive is
built from `pbs.project`, the gdata project of `base.dir.name`, and if89
(the PAUP module home). Verification: baseline characterization test
updated — every diff maps to an approved finding (F1 stop, F2 outpath
scripts, F3 base.dir.name/storage, F4 single-population, F5 shebang and
fixed -N, F7/F8 silence at verbose = 0) — 23 assertions pass. End-to-end
runs on testset.gs at verbose = 3 (standard and bash modes) completed
cleanly; generated PBS directives inspected. Sibling grep across the 8
clones: no package-code callers of `gl2paup.parsimony` — all clear.
PR: https://github.com/green-striped-gecko/dartR.base/pull/336

```json
{
  "function": "gl2paup.parsimony",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "medium", "rule": "STY1", "status": "applied", "change": 10}
  ],
  "coverage_skipped": ["test=TRUE path: stochastic, reviewed by reading", "Gadi execution of generated scripts: no server access", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": 336
}
```
