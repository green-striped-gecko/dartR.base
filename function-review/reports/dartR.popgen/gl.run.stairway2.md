# Review: gl.run.stairway2 (dartR.popgen)
- Family mode: analysis (wrapper around Stairway Plot 2, Java)
- Custodian: Bernd Gruber (git author; the roxygen header names no author or
  custodian, see F10). STY5: this report is the discussion record; changes
  approved by Luis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 7ddbca9 (origin/dev)
- Datasets: possums.gl[1:30, 1:200] after gl.filter.monomorphs (30
  individuals, 176 loci); testset.gs[1:20, 1:50]; platypus.gl (SEVERN_ABOVE)
  for the missing-data check
- External software: Stairway Plot 2 (`~/programs/stairway_plot_es`), Java
  26.0.1, macOS
- Baseline: tests/testthat/test-gl.run.stairway2.R (6 tests, 23
  expectations; the 2 real-run tests need STAIRWAY2_DIR and Java, otherwise
  skip)

## Verdict

**Standards: Needs work** — verbosity is hand-rolled and duplicated, Java's
output prints at `verbose = 0`, the SilicoDArT check only runs at
`verbose >= 2`, the global `future` plan is changed and not restored, and
failures stop with empty or misleading messages.
**Spec: Rework** — with the default `cleanup = TRUE` the function deletes
the whole R session temporary folder, including files it did not create
(the binaries downloaded by the documented example, other functions'
saved plots); `run = FALSE`, documented as the cluster mode, always errors;
`plot.file` makes a completed run fail.

What works well: the blueprint follows the Stairway Plot 2 template (folded
SFS from `gl.sfs`, bins 1..n, `nrand` defaults equal to the template's
(nseq-2)/4 ... nseq-2); runs with the same `seed` return identical
results; `parallel > 1` works on macOS and returns the same table layout.

## Findings

**F1 [BLOCKER, confidence: high] — cleanup deletes the whole session tempdir (FS7; principle: a function removes only files it created)**
`R/gl.run.stairway2.r:133,439` — the run folder is `tempd <- tempdir()`,
and `cleanup = TRUE` (default) calls `unlink(tempd, recursive = TRUE)`.
Failure scenario (baseline test 6, and probe): a marker file written to
`tempdir()` before the call is gone afterwards, and `tempdir()` itself no
longer exists. The documented example downloads the binaries with
`gl.download.binary()` into `tempdir()`, so the first run deletes them and
a second call stops with "Cannot find stairway_plot_es". Every other file
in the session tempdir is lost too: plot RDS files saved by other dartR
functions (their default `plot.dir`), and in dartR Shiny (which calls this
function) the session's uploaded files. `future` workers left in the
deleted folder print "shell-init: error retrieving current directory".
Proposed change: run in a folder of its own,
`tempd <- tempfile("stairway2_")`, and `unlink()` only that folder.

**F2 [HIGH, confidence: high] — `run = FALSE` always errors (DOC5, FS10)**
`R/gl.run.stairway2.r:357-462` — `res` and `p1` are only created inside
`if (run == TRUE)`, so the final `return(list(history = res, plot = p1))`
fails.
Failure scenario (baseline test 1): `run = FALSE` writes the blueprint to
`tempdir()/blueprint` and then stops with "object 'res' not found". The
user gets no return value and is not told where the blueprint is, so the
documented use ("might be useful to run on a cluster") does not work. The
documentation also gives the default as `FALSE`; the code default is
`TRUE`.
Proposed change: when `run = FALSE`, return `list(history = NULL, plot =
NULL)` plus the path of the run folder (blueprint and Stairway classes),
report that path at `verbose >= 2`, and never delete it; fix the
documented default to `TRUE`.
**Consequence: the returned list gains a third element (the run folder
path) for every call.**

**F3 [HIGH, confidence: high] — `plot.file` makes a completed run fail (PLT3)**
`R/gl.run.stairway2.r:439-449` — the plot is saved after F1's `unlink()`
has removed `plot.dir`'s default (`tempdir()`).
Failure scenario (probe): `plot.file = "swplot"` with default `plot.dir`
runs Stairway to the end (minutes with `nreps = 200`), then stops with
"cannot open the connection"; the history table is lost.
Proposed change: fixed by change 1 (the session tempdir is no longer
removed).

**F4 [MEDIUM, confidence: high] — SilicoDArT rejected only at verbose >= 2 (FS4, VRB3)**
`R/gl.run.stairway2.r:192-200` — the ploidy check sits inside
`if (verbose >= 2)`.
Failure scenario (baseline test 4): `testset.gs` at `verbose = 2` stops
with "Detected Presence/Absence (SilicoDArT) data"; at `verbose = 0` it
passes, `gl.sfs` builds a spectrum from presence/absence scores and
Stairway runs on it.
Proposed change: stop on `datatype == "SilicoDArT"` (from
`utils.check.datatype`) at every verbosity.
**Consequence: SilicoDArT input now errors at verbose 0 and 1 as well.**

**F5 [MEDIUM, confidence: high] — messaging ignores verbosity (FS2, FS3, VRB1, VRB2)**
`R/gl.run.stairway2.r:123-184,240,254,353-363,412,423,457` —
`utils.flag.start` is commented out and a second, older verbosity block
(lines 161-184) re-implements it; `gl.sfs` is called without `verbose`;
Java and the shell scripts print to the console unconditionally; all
messages use plain `cat()`; "Completed:" prints twice; "Check plots ... in
folder" names a folder that `cleanup` then deletes; line 250 `set.seed =
as.numeric(Sys.time())` assigns a variable and seeds nothing.
Failure scenario (probe): `verbose = 0` prints "Starting gl.sfs",
"Completed: gl.sfs" and about 50 lines of Stairway progress for
`nreps = 2` (it grows with `nreps`).
Proposed change: use `utils.flag.start` and remove the old block; pass
`verbose` to `gl.sfs`; send Java/shell output to the console only at
`verbose >= 3` (`system(..., ignore.stdout = , ignore.stderr = )`); use the
crayon helpers; print "Completed" once; point the "Check plots" message at
the kept folder only when `cleanup = FALSE`; drop the dead `set.seed =`
line.

**F6 [MEDIUM, confidence: high] — global `future` plan changed and not restored (principle: no side effects on session state)**
`R/gl.run.stairway2.r:373-374` — `future::plan(multisession, ...)` is set
and left in place; `detectCores() - 1` is 0 on a one-core machine.
Failure scenario (probe): after `parallel = 2`, `future::plan()` is still
a two-worker multisession plan, so the user's later `future`/`furrr` code
runs with it; its workers sit in the deleted folder (F1). On a one-core
machine `workers = 0` errors.
Proposed change: `oplan <- future::plan(...)` and restore it with
`on.exit(future::plan(oplan), add = TRUE)`; use `max(1, ...)` workers.

**F7 [MEDIUM, confidence: medium] — error-recovery and Windows code paths do not work (principle: platform-specific calls)**
`R/gl.run.stairway2.r:388-422` — `er` holds the value of the last `if` in
its block: on Linux/macOS that is `NULL`, so the low-memory retry never
runs; on Windows it is the exit status (length 1), so the plot step is
always run twice. The retry reads and edits `blueprint.plot.bat` even on
Linux/macOS, where only `.plot.sh` exists. With `parallel > 1` on Windows
the `MOVE /y` lines are rewritten to `cp`, which is not a Windows command.
Failure scenario: code reading; the Windows paths were not run (no Windows
machine). On macOS a failed Stairway step is not detected, and the user
gets `read.csv` "cannot open file ... Ne.final.summary".
Proposed change: check each `system()` exit status; move result files
with `file.rename()` in R instead of shell `mv`/`MOVE`; if the final
summary is missing, stop with a clear `error()` (and keep the run folder
for inspection); apply the `-Xmx` fallback to the script for the current
OS, or drop it.

**F8 [MEDIUM, confidence: high] — failures stop with empty or misleading messages (FS5, VRB2)**
`R/gl.run.stairway2.r:143-156` — a missing binary folder prints with
`cat()` and then calls `stop()` with no message; Java is never checked.
Failure scenario (baseline test 3): `stairway2.path` without
`stairway_plot_es` gives `Error:` with an empty message (the reason is on
stdout, lost in logs and in Shiny). Without Java the run continues until
`read.csv` fails on the missing summary file.
Proposed change: `stop(error(...))` naming the folder and pointing to
`gl.download.binary("stairway2")`; check `Sys.which("java")` before
running and stop with a message if it is empty.

**F9 [LOW, confidence: high] — plot ignores `plot.theme`; one column name has a leading space (PLT1)**
`R/gl.run.stairway2.r:434-436` — `+ plot.theme` is commented out; the
columns are renamed positionally with `" low75"` (leading space).
Failure scenario (baseline test 5): `res$history$low75` is `NULL`; the
user must write `` res$history$` low75` ``. Any `plot.theme` passed
(dartR Shiny passes one) has no effect.
Proposed change: add `+ plot.theme`; rename the column to `low75`; label
the x axis "Years ago".
**Consequence: the history column `" low75"` is renamed `low75`; the plot
takes `plot.theme`.**

**F10 [LOW, confidence: high] — documentation gaps and mismatches (DOC1, DOC2, DOC5, DOC7 (proposed rule))**
`R/gl.run.stairway2.r:1-83` — no `@author`/custodian and no `@family`;
`@references` given twice, with a stray blank line before `@export`;
`run` default documented as `FALSE`; `@return` says "generation, median,
low and high", the table has 11 columns (time in years, Ne median, 95% and
75% limits); `plot_title` default documented as `"Ne"+filename` (it is
`"Ne"`); `filename` "also used for the plot" (it is only the `popid`);
`stairway_plot_dir` is described as the output folder, but it is the name
of the class folder and any other value breaks the run (`java -cp
stairway_plot_es` is hard-coded); `pct_training` is a proportion, not a
percentage; `verbose` text is not the standard one; typos (minumum,
githubh, whethr, Extacting).
Failure scenario: a user sets `run = FALSE` or `stairway_plot_dir =
"out"` from the documentation and the call fails.
Proposed change: add Author(s)/Custodian and `@family`, fix the items
above, and state that `stairway_plot_dir` must stay `"stairway_plot_es"`
(or remove the argument: see change 8).

**F11 [LOW, confidence: high] — world-writable permissions (principle: least privilege)**
`R/gl.run.stairway2.r:350,365,419` — `chmod 777` on the blueprint and
scripts via `system()`.
Failure scenario: on a shared Linux server another user can edit the
script between creation and execution. Low risk because the folder is the
private session tempdir, but no reason for write permission to all.
Proposed change: `Sys.chmod(..., mode = "0755")` (fold into change 6).

**F12 [INFO, confidence: medium] — missing data enters the SFS as if every site had 2n sequences**
`R/gl.run.stairway2.r:222,240` — `nseq = 2 * nInd(x)`, and `gl.sfs` sums
allele counts with `na.rm = TRUE`, so a locus with missing calls is placed
in a lower frequency bin than its true one.
Failure scenario: platypus.gl SEVERN_ABOVE has 3.5% missing calls and 86 of
487 loci with at least one; the singleton share is 0.136 with all loci
and 0.150 with complete-call loci only. The direction and size of the bias
on the inferred history was not measured.
Proposed change: none here; handle in the `gl.sfs` review (pending) and,
meanwhile, mention in `@details` that loci should be filtered to full call
rate or imputed before running.

## Proposed changes

1. Run in a folder created for the call (`tempfile("stairway2_")`) and
   let `cleanup` remove only that folder (F1, F3).
   **Consequence: with `cleanup = FALSE` the Stairway files are in a
   subfolder of `tempdir()`, not in `tempdir()` itself.**
2. Make `run = FALSE` return `list(history = NULL, plot = NULL)` plus the
   run folder path, report the path, and fix the documented default (F2).
   **Consequence: the returned list gains a third element, the run folder
   path, for every call.**
3. Reject SilicoDArT input at every verbosity (F4).
   **Consequence: SilicoDArT input errors at verbose 0 and 1 as well.**
4. Verbosity and messaging: `utils.flag.start`, remove the duplicated
   block, pass `verbose` to `gl.sfs`, Java/shell output only at
   `verbose >= 3`, crayon helpers, single "Completed", drop the dead
   `set.seed =` line (F5).
5. Restore the user's `future` plan on exit; at least one worker (F6).
6. Failure handling: check exit statuses, `file.rename()` instead of
   shell moves, clear error when the summary is missing, `Sys.chmod(0755)`
   (F7, F11). Windows paths reviewed by reading only.
7. Clear errors for a missing binary folder and missing Java (F8).
8. Plot and table: apply `plot.theme`, rename `" low75"` to `low75`, x
   label "Years ago" (F9).
   **Consequence: the history column `" low75"` is renamed `low75`.**
9. Documentation fixes listed in F10, plus the missing-data note from F12
   (docs only).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour vs roxygen with the real binary on possums.gl (default,
  `run = FALSE`, `parallel = 2`, `plot.file`, `plot_title` with a space,
  same seed twice, wrong path, missing `mu`, SilicoDArT at verbose 0 and
  2) — run
- Callers (API3): dartR Shiny (`dartr2shiny/shiny_fun/Fun_gl.run.stairway2.R`)
  calls it with `parallel = TRUE` and `plot.theme`; no sibling dartR.*
  package calls it — run
- Windows code paths (F7): SKIPPED — no Windows machine; findings from code
  reading, confidence medium
- Linux: SKIPPED — same code path as macOS except `detectCores`; not run
- DAT1-DAT4, FS8: not applicable — the function does not modify or return
  a genlight
- DAT6 (FBM): SKIPPED — no FBM fixture; `gl.sfs` densifies with
  `as.matrix()`, to be covered in the `gl.sfs` review
- Numerical check of Stairway's estimates: not applicable — computed by
  the external program; checked only that `year = mutation_per_site / mu *
  gentime` (baseline test 5)
- Google Group / GitHub issues: not searched in this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |
| 9 | approved | Luis | |

## Outcome

Branch `review-stairway2` (from origin/dev 7ddbca9). Evidence: 
`tests/testthat/test-gl.run.stairway2.R`, 7 tests, 35 expectations, all
passing with the real binary (`STAIRWAY2_DIR=~/programs`, Java 26, macOS);
without the binary 4 tests pass and 3 skip. Every assertion changed from the
baseline is tagged with its approved change.

- 1 applied: runs in `tempfile("stairway2_")`; a marker file in `tempdir()`
  survives a run and `plot.file` saves after it (test 7).
- 2 applied: `run = FALSE` returns `history = NULL`, `plot = NULL`,
  `run.dir` (test 1); the list element is named `run.dir`.
- 3 applied: SilicoDArT errors at verbose 0 and 2 (test 4).
- 4 applied: no R output at `verbose = 0` (tests 1, 7); Java output checked
  by probe (a child R process at `verbose = 0`, `parallel = 2`: nothing
  printed between the start and end markers except R package-version load
  warnings from the environment).
- 5 applied: plan restored (test 7). Found while applying and fixed within
  the change: the plan is restored right after the parallel step, not only
  on exit, because workers shut down after `cleanup` printed
  "shell-init: error retrieving current directory".
- 6 applied: Stairbuilder failure stops with "could not create the run
  script" (test 5); moves done with `file.rename()`; `Sys.chmod(0755)`.
  Windows paths not run.
- 7 applied: missing folder error names the folder and
  `gl.download.binary` (test 3); Java checked with `Sys.which()` (not
  testable here: Java installed).
- 8 applied: column `low75` (test 6); `plot.theme` added; x label
  "Years ago".
- 9 applied: roxygen rewritten, `devtools::document()` run; NAMESPACE
  unchanged.
- Estimates unchanged: same data, `seed = 1`, `nreps = 4`: history before
  and after identical (`all.equal`, after renaming `low75`); `parallel = 2`
  identical to sequential.
- `devtools::check()`: 1 test failure in `test-gl.ld.haplotype.R:182`
  (unrelated, present on dev); install WARNING from packages built under
  R 4.4.3; NOTE on `gl.plot.popcluster` globals (fixed in #103). A
  `tests` WARNING from the first test version (callr/pkgload) was removed by
  rewriting that test in-process.
- NEWS entry added (escalation gate: return value and column name change).
  Callers: dartR Shiny only; it does not read `" low75"` or the list length.
- PR: #104

## Machine block

```json
{
  "function": "gl.run.stairway2",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "7ddbca9",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "FS7", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 1},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "principle: no session side effects", "status": "approved", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "medium", "rule": "principle: platform-specific calls", "status": "approved", "change": 6},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 7},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 8},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 9},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "principle: least privilege", "status": "approved", "change": 6},
    {"id": "F12", "severity": "INFO", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 9}
  ],
  "coverage_skipped": ["Windows paths: no Windows machine", "Linux: not run", "DAT6: no FBM fixture", "Google Group / issues not searched"],
  "status": "done",
  "pr": 104
}
```
