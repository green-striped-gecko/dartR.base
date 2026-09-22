# Review: gl.run.structure (dartR.popgen)

- Family mode: analysis (external-software wrapper: STRUCTURE 2.3.4)
- Date: 2026-09-18
- Reviewer: Claude (claude-fable-5-1), dartr-function-review v2.0.0
- Package commit: 31fb40c (`origin/dev`; reviewed in worktree branch
  `review-gl.run.structure`)
- Scope: `R/gl.run.structure.r` plus the private plumbing it calls,
  `R/utils.structure.run.r` and `R/utils.structure.genind2gtypes.r`
  (copied from strataG). `utils.structure.evanno` is read for its contract
  only.
- Datasets: testset.gl (dartR.data 1.2.5) restricted to its first three
  populations, monomorphs removed, call rate >= 0.9: 31 individuals, 18 loci;
  testset.gs for the SilicoDArT path
- Binary: `~/programs/structure` (STRUCTURE 2.3.4, macOS x86_64), run in a
  scratch working directory
- Baseline: `tests/testthat/test-gl.run.structure.R` (new file; 9 blocks,
  all passing; they run only when `STRUCTURE_EXEC` points to the executable
  and skip otherwise). Pins: result class and names pattern, one run per K,
  `summary` and `q.mat` shape, ids sorted alphabetically, `orig.pop` mapped
  by id, rows of `q.mat` summing to 1, reproducibility with
  `randomize = FALSE`, the error for fewer than three K values, run files in
  the working directory with `delete.files = FALSE`, RDS location ignoring
  `gl.set.wd()`, name truncation, failure on names with spaces, failure on an
  executable path with a space, input checks.
- Community check (Google Group, 2026-09-18): "Problem with gl.run.structure
  using large data" (STRUCTURE memory on 9,000 x 27,000; one reply notes the
  full executable path had to be given), "Structure speed and plotting error"
  (the error is in `gl.plot.structure`). No thread on the defects below.

## Verdicts

**Standards: Needs work** -- the backbone is out of order (executable check
and a package guard before the verbosity check), a Suggests package is used
unguarded while an Imports package is guarded with `return(-1)`, the run
files go to the user's working directory, `verbose` controls neither the
filter step nor STRUCTURE's output, and `plot.dir`'s default bypasses
`gl.set.wd()`.

**Spec: Needs work** -- STRUCTURE is driven correctly and the parsed
`q.mat` matches individuals and populations by id, but the function throws
away a finished run when `k.range` has fewer than three values, silently
collapses individual names longer than 11 characters to one label (the
documentation says it fails), and fails on names with spaces through its
own popflag lookup.

What works well: the STRUCTURE input files are written correctly
(MARKERNAMES, POPDATA, POPFLAG, MISSING -9, two rows per individual), the
model parameters are passed through faithfully, `randomize = FALSE` with a
seed reproduces results exactly, and `orig.pop` is mapped back to the
population names correctly.

## Findings

**F1 [HIGH, confidence: high] -- a finished STRUCTURE run is discarded when
`k.range` has fewer than three values (PLT3)**
`R/gl.run.structure.r:255-263`. `utils.structure.evanno(sr)` is called
unconditionally after the runs and stops with "must have at least two values
of k" (the check is `< 3`); the plot is built even when `plot.out = FALSE`
and `plot.file = NULL`.
Failure scenario: `gl.run.structure(x, k.range = 2)` or `k.range = 2:3`
runs STRUCTURE to completion (minutes to hours on real data) and then errors;
`sr` is never returned (observed for both).
Proposed change: return `sr` regardless; build, print and save the Evanno
plot only when it is asked for and at least three K values exist, otherwise
say so at `verbose >= 1`; wrap the plot step so a plotting failure cannot
lose the results.

**F2 [HIGH, confidence: high] -- run files are written to the working
directory, named by a timestamp, and left behind on failure (FS7)**
`R/utils.structure.run.r:443-458, 514-516, 535-537`. The run directory is
`<gtypes description>.structureRun` created relative to `getwd()`; the
description is a timestamp with one-second resolution; any existing directory
of that name is deleted first (`unlink(label, recursive = TRUE)`). STRUCTURE
writes `seed.txt` there too. Cleanup happens only on success.
Failure scenario: with `delete.files = FALSE` the directory
`gtypes.created.on.2026-09-18.19.27.57.structureRun` appears in the user's
working directory (observed); after a STRUCTURE failure the directory stays
even with `delete.files = TRUE` (observed); two calls started in the same
second in the same directory erase each other's files. Writing to the
working directory by default is against CRAN policy.
Proposed change: run in a per-call directory under `tempdir()`, delete it on
exit in every case, and when `delete.files = FALSE` keep the files under
`plot.dir` and print the path at `verbose >= 2`; label runs `k<K>.r<rep>`.

**F3 [HIGH, confidence: high] -- individual names with spaces produce an
NA popflag and STRUCTURE rejects the input (Spec axis: input construction)**
`R/utils.structure.run.r:225-239`. `popflag` is named with the original ids,
then `mutate()` rewrites `id` with underscores and looks `popflag` up by the
new ids, so every popflag is NA and is written as "NA".
Failure scenario: `indNames(x) <- sub("^(..)", "\\1 ", indNames(x))`; the
data file rows read `AA_011723 3 NA 0 1 ...` and STRUCTURE exits with code 1
("Probable error in the input file") (observed). The roxygen note tells the
user to remove spaces themselves.
Proposed change: see change 3 (index ids) which removes the lookup; at
minimum rename `popflag` with the sanitised ids.

**F4 [HIGH, confidence: high] -- names longer than 11 characters collapse
to one id without any error (DOC5 (proposed rule), Spec axis)**
`R/utils.structure.run.r:323-439` reads the ids back from STRUCTURE's
output, which truncates labels at 11 characters. `@details` states "The
function will fail if the names of individuals are not unique after
truncation".
Failure scenario: `indNames(x) <- paste0("individual_number_", 1:31)`: the
call succeeds and all 31 rows of `q.mat` carry `id = "individual_"`
(observed); `gl.plot.structure` uses `sr[[1]]$q.mat$id` as its labels, so
the plot cannot be matched to individuals.
Proposed change: write individuals to STRUCTURE as an index (1..n), restore
the full `indNames(x)` in `q.mat$id` and `prior.anc` names, and return
`q.mat` in genlight order. This also removes F3.

**F5 [MEDIUM, confidence: high] -- unquoted command line: an executable
path with a space is reported as "You do not have STRUCTURE installed"
(DEP2, platform handling)**
`R/utils.structure.run.r:496-509`. `paste0(exec, " -m ", ...)` is passed to
`system()` unquoted; a space in `exec` or in the working-directory path
splits the command.
Failure scenario: the binary copied to `<tmp>/my tools/structure` gives exit
code 127 and the misleading message (observed); a Windows install under
`Program Files` behaves the same. The Google Group thread on large data also
records a user who could only run after giving the full path.
Proposed change: `system2(exec, shQuote(args))`, check `file.exists` and
executability up front (already done for existence), and report the exit
status with the command on failure.

**F6 [MEDIUM, confidence: high] -- `plot.dir = tempdir()` as the default
bypasses `gl.set.wd()` (DOC5 (proposed rule), FS7)**
`R/gl.run.structure.r:169, 202`. The roxygen says "[default as specified by
the global working directory or tempdir()]", but because the default value
is `tempdir()` rather than `NULL`, `gl.check.wd()` never consults
`options()$dartR_wd`.
Failure scenario: `gl.set.wd(<dir>)` then `plot.file = "evanno_test"`: the
RDS lands in `tempdir()`, not in the set directory (observed).
Proposed change: `plot.dir = NULL`.

**F7 [MEDIUM, confidence: high] -- dependency guards: `purrr` (Imports) is
guarded with `return(-1)`, `tidyr` and `gridExtra` (Suggests) are not
guarded (DEP1)**
`R/gl.run.structure.r:176-184`; `tidyr::spread` at
`R/utils.structure.run.r:165`, `tidyr::gather` at
`R/utils.structure.genind2gtypes.r:196`, `gridExtra::grid.arrange` via
`utils.structure.evanno`.
Failure scenario: without `tidyr` the call fails inside the gtypes
conversion with a namespace error; the `purrr` branch can never trigger and
would return `-1` instead of stopping.
Proposed change: `stop(error())` guards for `tidyr` and `gridExtra`; drop
the `purrr` guard.

**F8 [MEDIUM, confidence: high] -- `verbose` controls neither the filter
step nor STRUCTURE's output (VRB1, VRB3)**
`R/gl.run.structure.r:224` calls `gl.filter.allna(x)` without `verbose`, so
its "Starting/Completed" lines print at `verbose = 0`; `system(cmd)` at
`R/utils.structure.run.r:501` streams STRUCTURE's full log to the console at
every verbosity (the roxygen says this "cannot be switched off currently";
`ignore.stdout` does it). No per-run progress line exists.
Failure scenario: `verbose = 0` still prints the filter messages and several
hundred STRUCTURE lines per K (observed).
Proposed change: pass `verbose` through; show STRUCTURE's own output only at
`verbose >= 3`; print "Running STRUCTURE K = 2, replicate 1 of 3" at
`verbose >= 2`.

**F9 [MEDIUM, confidence: medium] -- the Evanno layout draws the mean LnP(K)
panel twice and never the LnP'(K) panel (PLT1, Spec axis)**
`R/gl.run.structure.r:258, 261`: `ev$plots$mean.ln.k + ev$plots$mean.ln.k`
in both branches; `ev$plots$ln.pk` is never used.
Failure scenario: the printed and saved plot shows two identical top panels.
Proposed change: `mean.ln.k + ln.pk` on the top row.

**F10 [LOW, confidence: high] -- structure order, idioms and documentation
(FS2, FS3, FS5, DOC2, DOC5 (proposed rule), DOC7 (proposed rule))**
The executable check and the package guard run before SET VERBOSITY
(`:176-197`); `utils.flag.start(build = "Jody")`; typo "exex path"; the
datatype check could use `accept = "SNP"` instead of a second `stop()`;
`k.range` is documented `[required]` but `NULL` is handled (one run per
number of populations); `seed` is ignored unless `randomize = FALSE`
(STRUCTURE's rule, not stated); `verbose` text non-standard; `@author` has no
Custodian label; the `@details` claim about failing on truncated names is
false (F4); `@return` could state that `q.mat` rows are ordered by id.
Proposed change: docs and ordering only.

**F11 [LOW, confidence: high] -- deprecated idioms in the plumbing emit
warnings on every run (STY3, DEP2)**
`R/utils.structure.run.r:163-167, 234-240`: `.data$col` inside
`dplyr::select()` (deprecated in tidyselect 1.2.0), `tidyr::spread()` and
`tidyr::gather()` (superseded); `ggplot2::aes_string()` in
`utils.structure.evanno` (deprecated in ggplot2 3.0.0).
Failure scenario: a warning block on each call; a future tidyselect that
errors breaks the conversion.
Proposed change: use string column names in `select()`, `pivot_wider()` /
`pivot_longer()`; leave `aes_string()` to the `utils.structure.evanno`
review.

**F12 [INFO, confidence: medium] -- `structureRead` parses STRUCTURE's
text output by token position (Spec axis, robustness)**
`R/utils.structure.run.r:331-354`: `result[loc[1] + 6]` and
`result[loc + 2]` after `grep("Estimated")` / `grep("likelihood")` depend on
the 2.3.4 output layout. No change proposed; noting that only 2.3.4 is
supported and that a format change would produce NA summaries rather than an
error.

## Proposed changes

1. Return `sr` in every case; build, print and save the Evanno plot only
   when `plot.out` or `plot.file` asks for it and at least three K values
   exist, otherwise report at `verbose >= 1`; guard the plot step so it
   cannot lose the results (F1).
   **Consequence: calls with one or two K values now succeed and return the
   runs instead of erroring after STRUCTURE has finished.**
2. Run STRUCTURE in a per-call directory under `tempdir()`, removed on exit
   whether or not the run succeeds; with `delete.files = FALSE` keep the
   files under `plot.dir` and print the path; label runs `k<K>.r<rep>` (F2).
   **Consequence: with `delete.files = FALSE` the files move from the
   working directory to `plot.dir`; `names(sr)` change from
   `<timestamp>.structureRun.k2.r1` to `k2.r1`.**
3. Send individuals to STRUCTURE as an index, restore the full names in
   `q.mat$id` and `prior.anc`, and return `q.mat` rows in genlight order
   (F3, F4).
   **Consequence: `q.mat` rows are in `indNames(x)` order instead of
   alphabetical by id; names longer than 11 characters and names with
   spaces now work.**
4. Quote the command line (`system2` + `shQuote`) and report the exit status
   with the command; keep the existence check (F5).
5. `plot.dir = NULL` so `gl.set.wd()` is honoured (F6).
   **Consequence: with `gl.set.wd()` set, `plot.file` is saved there.**
6. Dependency guards: `stop(error())` for `tidyr` and `gridExtra`; drop the
   `purrr` guard (F7).
7. Verbosity: pass `verbose` to `gl.filter.allna`; STRUCTURE's own output
   only at `verbose >= 3`; one progress line per run at `verbose >= 2` (F8).
   **Consequence: at the default `verbose = 2` STRUCTURE's log no longer
   streams to the console.**
8. Evanno layout: top row `mean.ln.k + ln.pk` (F9).
9. Structure order (verbosity first), drop `build =`, `accept = "SNP"`,
   documentation fixes listed in F10. Docs and ordering only.
10. Replace deprecated tidyselect / tidyr idioms in the plumbing (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY -- run
- Spec: behaviour vs roxygen with the real binary on testset.gl -- run
- Input file inspection (`delete.files = FALSE`): data, mainparams,
  extraparams read and checked against the documented STRUCTURE format --
  run
- Population prior paths: `pop.prior = "usepopinfo"` and `"locprior"` with
  `k.range = 1:3` -- run; both return a 31 x 6 `q.mat` with rows summing to
  1, and `usepopinfo` returns 31 `prior.anc` matrices (Pop x Gen.0..2)
- Windows path handling: reasoned from `system()` quoting; not executed
- FBM path (DAT6): SKIPPED -- no fixture; the conversion goes through
  `gl2gi` and densifies in any case
- `num.k.rep > 1` (delta K branch): run once (`k.range = 1:3`,
  `num.k.rep = 2`, 6 runs returned), plot not inspected visually
- Downstream callers: `gl.plot.structure`, `gl.map.structure`, `gl.evanno`
  in this package consume `sr` positionally (`sr[[1]]$q.mat`, `x$summary`);
  dartr2shiny wraps the function. Change 2 (names) and change 3 (row order)
  were checked against `gl.plot.structure`, which orders by `Label` itself
  (`R/gl.plot.structure.r:282`) and uses `x[[2]]` positionally.

## Approval

All ten changes approved by Luis on 2026-09-18 (AskUserQuestion selections;
changes 1, 2, 3, 5 and 7 approved with their stated consequences).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | behaviour change accepted |
| 2 | approved | Luis | file location and run-name change accepted |
| 3 | approved | Luis | q.mat row order change accepted |
| 4 | approved | Luis |  |
| 5 | approved | Luis | plot.file location change accepted |
| 6 | approved | Luis |  |
| 7 | approved | Luis | STRUCTURE log hidden below verbose 3 accepted |
| 8 | approved | Luis |  |
| 9 | approved | Luis |  |
| 10 | approved | Luis |  |

## Outcome

Applied on branch `review-gl.run.structure` (cut from `origin/dev`), all ten
changes, in `R/gl.run.structure.r`, `R/utils.structure.run.r` and
`R/utils.structure.genind2gtypes.r`. Evidence:

- Characterization test: the pre-change baseline (9 blocks) run against the
  revised code moved exactly the pins tagged with approved findings: run
  names (change 2), `q.mat` order (change 3), the error for fewer than three
  K values (change 1), run files in the working directory (change 2), RDS
  location (change 5), truncated names and names with spaces (change 3),
  executable path with a space (change 4). The updated file (12 blocks, 51
  assertions) passes with `STRUCTURE_EXEC` set; it skips without it.
- Change 1: `k.range = 2` returns one run with the note "The Evanno plot
  needs at least three values of K; 1 found"; `k.range = 2:3` returns two.
- Change 2: nothing is written to the working directory in any call,
  including a forced failure (fake executable exiting 2: "STRUCTURE exited
  with status 2 for run k2.r1 ... fake structure: cannot allocate", no
  `structureRun_*` directory left in `tempdir()`); `delete.files = FALSE`
  keeps `k1.r1_data`, `_mainparams`, `_extraparams`, `_out_f`, `_log` and
  `seed.txt` under `<plot.dir>/structureRun_<timestamp>`, with the paths in
  each run's `files`.
- Change 3: `q.mat$id` equals `indNames(x)` and `orig.pop` equals `pop(x)`
  without reordering; names like "individual number 12" come back intact;
  with `pop.prior = "usepopinfo"` and an unnamed `popflag` with zeros at
  individuals 1 and 5, the data file carries popflag 0 for indices 1 and 5
  only, and `prior.anc` is named by the 29 flagged individuals in genlight
  order (the original code also returned 29 entries).
- Change 4: the binary copied to `<tmp>/my tools/structure` runs.
- Change 5: with `options(dartR_wd)` set, `plot.file` lands there and reads
  back as a patchwork object.
- Change 7: `verbose = 0` prints 0 R-level lines (before: the two
  `gl.filter.allna` lines plus STRUCTURE's log); `verbose = 2` prints one
  "Running STRUCTURE: K = k, replicate r (run i of n)" line per run;
  `verbose = 3` streams STRUCTURE's output.
- Change 8: with `num.k.rep = 2` the layout builds with the delta K panel;
  with one replicate the bottom row is LnP''(K) alone.
- `devtools::document()` regenerated `man/gl.run.structure.Rd` and
  `man/utils.structure.run.Rd`; NAMESPACE unchanged.
- `devtools::check()`: 0 errors with the STRUCTURE tests enabled; the 1
  warning (dependencies built under R 4.4.3 at install) and 4 notes
  (worktree `.git` file, timestamps, `gl.plot.popcluster` globals, NEWS
  parsing) are pre-existing and unrelated.
- Addendum F13 (found while applying change 3): an unnamed `popflag` was
  matched to the ids in sorted order, not in `indNames(x)` order, because
  `names(popflag) <- unique(g$data$id)` ran on a data.table keyed by id.
  The index mapping in change 3 names an unnamed `popflag` by
  `indNames(x)`; approval of this reading is requested with the PR.
- PR: green-striped-gecko/dartR.popgen#92, branch `review-gl.run.structure`
  (cut from `origin/dev`), commit 81edbea.

```json
{
  "function": "gl.run.structure",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "31fb40c",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "spec-input-construction", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DOC5/spec", "status": "approved", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DEP2/platform", "status": "approved", "change": 4},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5/FS7", "status": "approved", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 6},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1/VRB3", "status": "approved", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "medium", "rule": "PLT1/spec", "status": "approved", "change": 8},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "FS2/FS3/FS5/DOC2/DOC5/DOC7", "status": "approved", "change": 9},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "STY3/DEP2", "status": "approved", "change": 10},
    {"id": "F12", "severity": "INFO", "confidence": "medium", "rule": "spec-robustness", "status": "no-change", "change": null}
  ],
  "coverage_skipped": ["Windows quoting not executed", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": 92
}
```
