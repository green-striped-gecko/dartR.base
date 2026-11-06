# Review: gl.LDNe (dartR.popgen)

- Family mode: analysis (external-tool wrapper returning a list of tables;
  the genlight is not modified, so no history)
- Date: 2026-09-18
- Reviewer: Claude (claude-fable-5-1), dartr-function-review v2.0.0
- Package commit: a4fd6ed (dev_luis, synced with origin/dev)
- Datasets: possums.gl[1:60, 1:100] (the roxygen example: populations A and
  B of 30, 100 SNPs); subsets with a two-individual population placed
  first, second and in the middle of three; testset.gs for the SilicoDArT
  path. NeEstimator V2.1 Mac binary (`Ne2-1M`, x86_64 under Rosetta on an
  arm64 Mac).
- Baseline: tests/testthat/test-gl.LDNe.R (dartR.popgen), 10 tests, all
  passing on the reviewed state; skipped unless `NEEST_DIR` names the
  folder holding the binary. Six tests pin defects (F1, F3, F4, F5, F6, F7,
  F8, F9).

**Standards: Needs work** - the house structure is present, but four
argument checks print a message and carry on, `verbose` is ignored by the
table print, the gl2genepop call and the binary, and two error paths raise
an empty error.
**Spec: Needs work** - the LD estimate, parametric CIs, Waples correction
(checked against eq 1a) and `pairing = "separate"` are right, but the
documented `mating = "monogamy"` value cannot run, and a first population
with fewer than three individuals takes the next population's jackknife
CIs.

What works well: the tables agree line by line with NeEstimator's own
output for the example, the Waples correction reproduces eq 1a to the
printed precision, `pairing = "separate"` restricts the comparisons as the
ChrMap file intends, and the working directory is restored on error.

## Findings

**F1 [BLOCKER, confidence: high] - `mating = "monogamy"`, the documented
value, errors (DOC5)**
`R/gl.LDNe.r:207-213` - `pmatch(mating, c("random", "mono"))` looks for
`mating` as a prefix of the table entries, so "monogamy" (longer than
"mono") is NA and the next `if` fails with "missing value where TRUE/FALSE
needed". Only the undocumented "mono" runs.
Failure scenario (reproduced): `gl.LDNe(pops, mating = "monogamy")` stops.
The dartr2shiny module offers exactly `random`/`monogamy`, so the Shiny
monogamy option has never worked.
Proposed change: `match.arg(mating, c("random", "monogamy"))`, keeping
"mono" as an accepted abbreviation, mapped to NeEstimator's 0/1.

**F2 [BLOCKER, confidence: high] - jackknife CIs shift to the wrong
population when the first population has fewer than three individuals
(DAT2 principle: values must track the population they describe)**
`R/gl.LDNe.r:322-329` - NeEstimator prints no jackknife line for a
population of fewer than three individuals (verified in the raw output).
The fix-up inserts NA at position `r` with `c(CI[1:(r - 1)], NA, CI[r:n])`;
for `r = 1` the index `1:0` is `c(1, 0)`, which selects element 1, so the
result is (pop2, NA, pop2).
Failure scenario (reproduced): populations A (2 individuals) and B (30):
A's table shows jackknife CI 11.0 to 24.7, which are B's values, and B
shows NA. With the small population second or in the middle of three the
placement is correct. The April 2025 forum thread "Error during gl.LDNe"
reports CIs that do not contain the estimate; this defect produces that
symptom for the affected layout.
Proposed change: `append(CI, NA, after = r - 1)`.
**Consequence: jackknife CI values change for datasets whose first
population has fewer than three individuals (they were another
population's values).**

**F3 [HIGH, confidence: high] - Waples correction arguments are not
validated; an invalid type applies eq 1b (DOC5, FS5)**
`R/gl.LDNe.r:141-150`, `:408-415` - both checks use `message()` and
continue; the second test reads `length(Waples.correction.value == 1)`,
which is always at least 1. Any type other than "nChromosomes" falls into
the `else` branch and applies the genome-length formula.
Failure scenario (reproduced): `Waples.correction = "foo",
Waples.correction.value = 10` prints the message and appends five
"Waples' corrected" rows with negative Ne (-37 for 15). A missing value
stops with "non-numeric argument to mathematical function"; a length-2
value recycles. The dartr2shiny module passes the string "NULL" as its
default type, so with the code as it stands every default Shiny run
appends eq 1b rows (Shiny path inferred from `Fun_gl.LDNe.R`, not run).
Proposed change: `match.arg` on the type and a real length-1 numeric check
on the value, both with `stop(error())`.
**Consequence: calls with an unrecognised correction type, including the
Shiny default "NULL" string, stop instead of appending rows; dartr2shiny
must pass NULL.**

**F4 [HIGH, confidence: high] - `plot.file` without `plot.out` errors
(PLT3)**
`R/gl.LDNe.r:436-529` - the plot object `p3` exists only inside
`if (plot.out)`, but `utils.plot.save(p3, ...)` runs whenever `plot.file`
is set.
Failure scenario (reproduced): `plot.out = FALSE, plot.file = "ne"` stops
with "object 'p3' not found" after the binary has run.
Proposed change: build the plot when `plot.out` or `plot.file` is set;
print it only when `plot.out`.

**F5 [MEDIUM, confidence: high] - `plot_colors_pop` must have exactly one
colour per population, and the documented default is not the real one
(DOC5, PLT1)**
`R/gl.LDNe.r:452` - `rep(plot_colors_pop, each = nFreq)` is assigned to a
column of `nPop * nFreq` rows, so any other length errors. The roxygen says
"[default discrete_palette]" (a palette function); the signature default
is `gl.select.colors(x)`. Passing the documented palette function errors
as well.
Failure scenario (reproduced): four colours for two populations stop with
"replacement has 8 rows, data has 4". The 2025 forum thread reports
"replacement has 15 rows, data has 18" on this line.
Proposed change: accept a palette function or a vector, take the first
`nPop` colours, stop with a clear message when fewer are supplied; fix the
doc.

**F6 [MEDIUM, confidence: high] - `naive = TRUE` drops the population
names from the result (DOC5)**
`R/gl.LDNe.r:422-433` - `lapply(seq_along(pop_list), ...)` returns an
unnamed list.
Failure scenario (reproduced): `res$A` is NULL; code indexing by
population name (as the roxygen example implies) breaks only when `naive`
is on.
Proposed change: keep `names(pop_list)`.

**F7 [MEDIUM, confidence: high] - `pairing` and data-type checks message
and continue (FS5, VRB2)**
`R/gl.LDNe.r:126-137` - both use `message(error())`.
Failure scenario (reproduced): `pairing = "foo"` reaches `option[10] <-
setPairs` and stops with "object 'setPairs' not found"; SilicoDArT input
prints the message, then stops inside gl2genepop with a second message.
Proposed change: `match.arg(pairing)` and `stop(error())` for the data
type.

**F8 [MEDIUM, confidence: high] - output ignores `verbose` (VRB1, VRB3)**
`R/gl.LDNe.r:190`, `:286`, `:519` - gl2genepop is called without
`verbose` (it prints at its default 2), the binary's stdout is never
redirected, and `print(pop_list)` is unconditional.
Failure scenario (reproduced): `verbose = 0` prints 30 lines from R plus
the binary's full log (locus list, per-population estimates).
Proposed change: `gl2genepop(..., verbose = 0)`; `system(cmd,
ignore.stdout = verbose < 3)`; print the tables at `verbose >= 2`.
**Consequence: at `verbose` 0 or 1 the tables are no longer printed (the
list is still returned invisibly).**

**F9 [MEDIUM, confidence: high] - the saved output file is never
refreshed (FS7)**
`R/gl.LDNe.r:386` - `file.copy(outfile, file.path(outpath, outfile))`
without `overwrite = TRUE` returns FALSE when the target exists.
Failure scenario (reproduced): two runs with the same `outfile` and
`outpath`, the second with three critical values; the saved file still
holds the first run (same size, same mtime) while the message says "The
results are saved in".
Proposed change: `overwrite = TRUE`.

**F10 [MEDIUM, confidence: medium] - fixed file names in `tempdir()`, no
exit-status check (FS5; platform lens)**
`R/gl.LDNe.r:190-286` - `dummy.gen`, `infodummy`, `option` and the outfile
are written with fixed names directly in `tempdir()`, the binary's exit
status is discarded, and `read.delim(outfile)` reads whatever is there.
Failure scenario: forked parallel calls (`mclapply` shares the parent's
`tempdir()`) overwrite each other's `dummy.gen`; the forum thread "Issues
parallelizing gl.LDNe()" reports "Data of sample XXX end too soon" on
random samples, which is the binary reading a file another fork is
rewriting. A binary failure after a successful run in the same session
would return the previous outfile. Not reproduced here (the binary did not
fail in any probe).
Proposed change: run in a fresh `tempfile("LDNe_")` directory under
`tempdir()`, check the exit status and that the outfile exists, stop
otherwise; copy the outfile to `outpath` as now.

**F11 [LOW, confidence: high] - error paths raise an empty error, and an
unknown OS is unhandled (VRB2, FS5)**
`R/gl.LDNe.r:211-213`, `:270-280` - `cat(error(...)); stop()` prints the
text but the condition message is empty, so `tryCatch` and Shiny toasts
show nothing. `:248-261` - on an OS other than Windows/Linux/Darwin `prog`
is undefined.
Proposed change: `stop(error(...))`; stop with a message on an unsupported
OS.

**F12 [LOW, confidence: high] - documentation drift (DOC1, DOC2, DOC5,
DOC7 proposed, FS3)**
`R/gl.LDNe.r:63` - `@return` says "Dataframe"; the function returns,
invisibly, a named list with one data frame per population, 10 rows plus
5 with `Waples.correction` and 1 with `naive`. `:57-59` - `plot.file` has
a stray leftover line; `plot.dir` "[default = working directory]" but the
default is NULL (tempdir via `gl.check.wd`). `:55-56` - default colours
(F5). `:40` - "monogamy" (F1). `:60-62` - non-standard verbose text.
`:81-85` - the second example lacks a closing parenthesis and writes to
"./TestNe" in the working directory. `:64` - no Author(s) line (DOC7).
`:116` - `build = "Jody"` (FS3). `man/gl.LDNe.Rd` is in sync (DOC4).
Proposed change: docs-only edits, then `devtools::document()`.

**F13 [INFO, confidence: high] - parser notes, no change proposed**
NeEstimator appends a trailing "0+" column to every table (raw: "0.050
0+ No S* 0+"); the whole-column `duplicated()` at `:376` removes it, which
is why the "Frequency N" count uses `!duplicated(freq)`. The 2023-2025
forum error "factor level [2] is duplicated" came from a `levels =` argument
now commented out at `:450` and no longer occurs. Temporary-disk use by the
binary (forum thread "Temporary files overload") is documented behaviour.

## Proposed changes

1. Accept `mating = "monogamy"` (and "mono") via `match.arg` (F1).
2. Insert the jackknife NA with `append(after = r - 1)` (F2).
   **Consequence: jackknife CI values change when the first population has
   fewer than three individuals.**
3. Validate `Waples.correction`, `Waples.correction.value`, `pairing` and
   the data type with `stop(error())` (F3, F7). **Consequence: an
   unrecognised correction type, including the Shiny default "NULL"
   string, stops; dartr2shiny must pass NULL.**
4. Build the plot whenever `plot.out` or `plot.file` is set; print only
   when `plot.out` (F4).
5. Accept a palette function or vector for `plot_colors_pop`, use the first
   `nPop` colours, clear error when fewer (F5).
6. Keep population names when `naive = TRUE` (F6).
7. Honour `verbose`: silent gl2genepop, binary stdout at `verbose >= 3`,
   tables at `verbose >= 2` (F8). **Consequence: no table print at verbose
   0 or 1.**
8. Run in a fresh per-call directory, check the exit status and the
   outfile, overwrite the saved copy (F9, F10).
9. `stop(error())` on the mating and executable error paths; explicit stop
   on an unsupported OS (F11).
10. Documentation fixes and `devtools::document()`; drop `build =` (F12).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API - run
- Spec: 20 probes on possums.gl[1:60, 1:100] and subsets (default, verbose
  0, mating values, critical sets, plot.out/plot.file combinations,
  colours, naive, Waples type/value cases, pairing, small population first,
  second and middle, non-alphabetical population levels, SilicoDArT,
  repeated outfile, missing executable, underscore population names) - run
- Numerical: Waples eq 1a recomputed from the Ne row (22.415 vs printed
  22.4); tables compared with the binary's own log for the example - run
- DOC4 Rd sync: roxygenise in a scratch copy, no diff - run
- Callers: `gl.check.panel` (dartR.popgen) calls with `verbose = 0`,
  `mating = "random"` and reads `x$\`Frequency 1\`[6]`; dartr2shiny
  `Fun_gl.LDNe.R` passes `mating` in {random, monogamy} and
  `Waples.correction` in {"NULL", nChromosomes, genomeLength} - run
- Google Group: threads "Error during gl.LNDe" (2TecbYjEu1E, fixed),
  "Error during gl.LDNe" (qIgrNXlBgkQ, F5 symptom and CI complaint),
  "Issues parallelizing gl.LDNe()" (QVT2w6115y0, F10), "Temporary files
  overload" (WYFkUlWHzrM, design); GitHub dartR#673 (pairing option,
  present) - read
- Windows and Linux binaries: SKIPPED - only the Mac binary is available;
  the Windows CI complaint in qIgrNXlBgkQ is a binary matter, out of scope
- Binary failure path (F10 stale outfile): SKIPPED - could not make the
  binary fail on the reference data; confidence set to medium
- FBM path (DAT6): not applicable - the conversion is delegated to
  gl2genepop

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |
| 8 | approved | Luis |  |
| 9 | approved | Luis |  |
| 10 | approved | Luis |  |

## Outcome

Applied on branch `review-gl.LDNe`, cut from `origin/dev` in a separate
worktree so that the shared checkout's index stayed out of the commit.

- Change 1 (F1): `mating` accepts "random" and "monogamy" through
  `match.arg`, with "mono" mapped to "monogamy" before the check; the
  NeEstimator code is now `ifelse(mating == "random", 0, 1)`. Evidence:
  "monogamy", "mono" and "random" all run, "monogamy" and "mono" give the
  same Ne (31.6, 30.1), "clonal" stops.
- Change 2 (F2): `append(CI, NA, after = r - 1)`. Evidence: with a
  2-individual population A first, A's jackknife limits are NA and B keeps
  11.0 to 24.7; with the small population second or in the middle of three
  the values are unchanged from the baseline.
- Change 3 (F3, F7): `pairing`, `mating`, `Waples.correction`, its value
  and the data type are all rejected with `stop(error())` before the binary
  runs. Evidence: five error expectations, including the string "NULL" as
  the correction type and a length-2 value.
- Change 4 (F4): the plot is built when `plot.out` or `plot.file` is set.
  Evidence: `plot.out = FALSE, plot.file = "ldne"` writes `ldne.RDS` and
  `ldne_tab.RDS`.
- Change 5 (F5): a palette function is called with the number of
  populations, a longer vector is truncated, fewer colours than populations
  stops with the counts. Evidence: three expectations.
- Change 6 (F6): `names(pop_list) <- pops` after the naive block. Evidence:
  `res$A` is a data frame with `naive = TRUE`.
- Change 7 (F8): `gl2genepop` silent below `verbose 3`, `system(cmd,
  ignore.stdout = verbose < 3)`, tables at `verbose >= 2`. Evidence:
  `capture.output()` at `verbose = 0` has length 0 (was 30 lines plus the
  binary's log); at `verbose = 2` the tables print and the gl2genepop
  messages do not.
- Change 8 (F9, F10): each call runs in `tempfile("LDNe_")` under
  `tempdir()`, removed on exit; the exit status and the existence of the
  output file are checked; the copy to `outpath` uses `overwrite = TRUE`.
  Evidence: a second run with three critical values enlarges the saved
  file (was byte-identical to the first run); no `LDNe_` directory is left
  behind.
- Change 9 (F11): `stop(error())` on the missing-executable path; an
  unsupported operating system stops with a message. Evidence: the missing
  binary now raises "Cannot find Ne2-1M ..." as a readable condition
  message (it was empty).
- Change 10 (F12): documentation corrected as listed and
  `devtools::document()` run (31 insertions, 17 deletions in
  `man/gl.LDNe.Rd`); `build = "Jody"` dropped from `utils.flag.start`.

Addenda found while applying, not part of the approved list:

- **F14 [MEDIUM]** - the per-call run directory is removed when the call
  ends, so NeEstimator's byproduct files (`dummyLoc.txt`, `dummyBur.txt`
  and the tabular `<outfile>xLD.txt`) no longer remain in `tempdir()`
  afterwards. Without the removal, every call would leave a full set of
  files behind rather than overwriting one set, which is worse for the
  disk-space problem reported in the forum thread "Temporary files overload
  with Ne2-1L/gl.LDNe()". No dartRverse function reads these files (grep
  over the six packages and dartr2shiny). The documented output,
  `file.path(outpath, outfile)`, is copied out before the removal.
- **F15 [LOW]** - `outpath` is resolved with `normalizePath()` before the
  working directory changes. The roxygen says `outpath = '.'` directs
  output to the working directory; because the copy ran from inside
  `tempdir()`, it landed in `tempdir()` instead. Evidence: with
  `outpath = "."` the file now appears in the caller's working directory.

Verification:

- `tests/testthat/test-gl.LDNe.R`: 11 tests, 52 expectations, all passing.
  Skipped unless `NEEST_DIR` names the folder holding the binary for the
  running OS.
- Unchanged from the baseline: the estimates, parametric CIs, harmonic mean
  sample sizes, independent-comparison counts and column layout for
  `possums.gl[1:60, 1:100]`; the Waples eq 1a and eq 1b values; the
  `pairing = "separate"` comparison counts.
- Changed, each mapped to an approved change or an addendum: "monogamy"
  runs (1); jackknife limits for a small first population (2); five
  argument errors (3); `plot.file` alone saves (4); colour handling (5);
  names with `naive` (6); silence at `verbose 0` (7); refreshed output file
  and no leftover directory (8, F14); readable error message (9); relative
  `outpath` (F15).
- `gl.check.panel`, the only in-package caller, reads
  `x$\`Frequency 1\`[6]` from the result; that access still returns 13.5
  and 16.1 for the two populations.
- `devtools::check()`: 1 error, in `test-gl.ld.haplotype.R:182` ("verbose = 0
  is silent"), which expects no output and gets three "Processing genlight
  object with SNP data" lines. It reproduces on a pristine `origin/dev`
  worktree, so it is not from this change and was left alone; a nested call in
  that function is not passing `verbose = 0` to the datatype check. The
  remaining 1 warning and 4 notes come from dependencies built under a newer R
  and from the checkout, not from this change. The gl.LDNe tests skip under
  `R CMD check` because `NEEST_DIR` is unset there.
- Addenda F14 and F15 approved by Luis with the pre-push OK.
- PR: green-striped-gecko/dartR.popgen#93, commit c7c612f on branch
  `review-gl.LDNe` (cut from origin/dev be0a8a1).

```json
{
  "function": "gl.LDNe",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "a4fd6ed",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "approved", "change": 8},
    {"id": "F10", "severity": "MEDIUM", "confidence": "medium", "rule": "FS5", "status": "approved", "change": 8},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 9},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 10},
    {"id": "F13", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "note", "change": null}
  ],
  "coverage_skipped": ["Windows/Linux binaries: not available", "binary failure path: not reproducible", "DAT6: not applicable"],
  "status": "pr-open",
  "pr": 93
}
```
