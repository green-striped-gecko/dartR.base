# Review: gl.infer.sex (dartR.sexlinked)

- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 179080b (`origin/dev`, version 1.2.2, loaded with `devtools::load_all()` under R 4.4.2), branch `review-gl.infer.sex`
- Datasets: LBP, platypus.gl (dartR.data 1.2.5), both passed through `gl.keep.sexlinked()` from the same commit
- Baseline: `tests/testthat/test-gl.infer.sex.R` (20 expectations, all pass on 179080b)
- Related: reviews of `gl.report.sexlinked` (PR #33), `gl.keep.sexlinked` (PR #34), `gl.filter.sexlinked` (PR #35)

## Verdict

**Standards: Needs work**: the function has no `verbose` argument, prints unconditionally, never checks its input, and resets the caller's random number stream.

**Spec: Needs work**: an individual with no usable genotypes is reported as `*F`, a one-against-one disagreement is always reported as `*F`, and a `system` that differs from the one used in `gl.keep.sexlinked()` gives every individual `*F` with no error. The k-means step assumes both sexes are present, which the documentation does not say.

What works well: on LBP with `system = "xy"`, 345 of the 373 individuals with a recorded sex get a definite assignment and 344 of those match the recorded sex. Individuals with conflicting evidence are flagged with `*` as the documentation promises.

## Findings

**F1 [HIGH, confidence: high]: no check that the input matches `system` (FS5, DAT5)**
`R/gl.infer.sex.r:107-131`: the subsets are read by name (`$y.linked`, `$w.linked`, ...). A name that is not in the list returns `NULL`, and `NULL` is treated as "no loci of this type".
Failure scenario: `k <- gl.keep.sexlinked(LBP, system = "xy")` then `gl.infer.sex(k, system = "zw")` returns 376 rows, all `agreed.sex = "*F"`, and `id` replaced by `"1"` to `"376"`. Nothing is printed except the per-type "not enough loci" notes and `***SUCCESS***`.
Proposed change: at the start, stop with an error if `gl_sexlinked` is not a list containing `results.table` and the two subset names for the chosen `system`, and say which system the list was built for.

**F2 [HIGH, confidence: high]: individuals with no information are assigned `*F` (spec)**
`R/gl.infer.sex.r:193-203`: when all three preliminary assignments are `NA`, both counts in `yy` are 0, `order()` keeps `F` first, and the function returns `"*F"`.
Failure scenario: an individual with every X-linked genotype missing, in a dataset with no Y-linked loci, gets `agreed.sex = "*F"` (tested on LBP, row `Y2` with its X-linked genotypes set to `NA`). The same happens to all individuals in the F1 scenario.
Proposed change: return `NA` when no preliminary assignment is available.

**F3 [MEDIUM, confidence: high]: ties are always resolved as `*F` (spec)**
`R/gl.infer.sex.r:195-201`: with one `F` and one `M`, both counts are 1 and the stable `order()` puts `F` first.
Failure scenario: on LBP (`xy`), 28 individuals have a Y-linked call of `F` and an X-linked call of `M`; all 28 get `*F`. Had the calls been reversed, they would still get `*F`. The label reflects alphabetical order, not evidence. The documentation says `*` entries are for human checking, so the letter after `*` suggests a lean that does not exist.
Proposed change: return `"*?"` when the counts tie. Entries still start with `*`, so code that selects indefinite assignments with `grepl("^\\*", ...)` keeps working.

**F4 [MEDIUM, confidence: high]: k-means always returns two sexes (spec, DOC5)**
`R/gl.infer.sex.r:294, 352`: `kmeans(Z, 2)` splits any sample into two groups.
Failure scenario: 162 LBP females only, X-linked loci only: 68 are assigned `M` as definite assignments. With gametologs only (`zw` subsets), 64 of 162 females are assigned `M`. Samples from one sex are plausible (a subset by population, a captive colony).
Proposed change: document in `@details` that the X/Z-linked and gametolog assignments assume both sexes are present, and that a single-sex sample is split in two regardless. A warning based on cluster separation was tested and rejected: the centre heterozygosity of the gametolog clusters is 0.07 and 0.67 for females only, against 0.02 and 0.68 for the mixed sample, so no threshold separates the two cases.

**F5 [MEDIUM, confidence: high]: opaque error when no sex-linked loci exist (FS5)**
`R/gl.infer.sex.r:124-131, 145`: `n_ind` stays `NULL` when all three subsets are `NULL`.
Failure scenario: `gl.keep.sexlinked()` output with all subsets `NULL` fails with `Error in rep(NA, n_ind) : invalid 'times' argument`.
Proposed change: covered by the input check in change 1: stop with "No sex-linked loci found in gl_sexlinked".

**F6 [MEDIUM, confidence: high]: the caller's random number stream is reset (principle: no side effects on global state)**
`R/gl.infer.sex.r:293, 351`: `set.seed(seed)` is called on the global stream and not restored.
Failure scenario: `set.seed(1); gl.infer.sex(k, "xy", seed = 100); runif(1)` gives a different number from `set.seed(1); runif(1)`. A script that sets a seed, calls `gl.infer.sex()`, then simulates or bootstraps does not get the stream it asked for.
Proposed change: save `.Random.seed` on entry (if it exists) and restore it with `on.exit()`.

**F7 [MEDIUM, confidence: high]: no `verbose` argument; messages are ungated (FS2, FS3, FS9, VRB1, VRB2, VRB3)**
`R/gl.infer.sex.r:86-88, 137-186, 242`: every run prints the per-type notes and `***SUCCESS***` with `message()`. `verbose = 0` is rejected as an unused argument. The notes contain a typo ("gamtologous") and the multi-line strings print their source indentation.
Failure scenario: a pipeline running with `gl.set.verbosity(0)` still prints four lines per call and cannot silence them without `suppressMessages()`.
Proposed change: add `verbose = NULL` as the last argument, with `gl.check.verbosity()`, `utils.flag.start()`, notes at `verbose >= 2` via `cat(report(...))` / `cat(warn(...))`, and the standard "Completed:" line at `verbose >= 1`. Fix the typo and the indentation.
Downstream: `dartr2shiny` wraps this function (`config/functions.csv`, `input_generator/dartR.sexlinked/gl.infer.sex.r`). A new argument at the end does not break existing calls, but the generator config needs a row for it (API3).

**F8 [LOW, confidence: high]: rows are matched by individual name (DAT2)**
`R/gl.infer.sex.r:311-325, 369-383`: X/Z-linked and gametolog assignments are written back by looking up row names.
Failure scenario: with two individuals sharing a name (adegenet allows it), every X-linked assignment becomes `NA` and `id` becomes row numbers.
Proposed change: keep assignments in row order (assign `NA` at the rows removed by `na.omit()` by position). This also removes the `"STOP"` placeholder loop.

**F9 [MEDIUM, confidence: high]: documentation does not match behaviour (DOC5)**
`R/gl.infer.sex.r:1-80`:
- Only the 5 gametologs with the smallest adjusted p-value are used, and `#Het.g`/`#Hom.g` count only those 5. Not documented.
- The minimum numbers of loci (2 X/Z-linked, 5 gametologs) below which a type is skipped are not documented.
- `@details` says indefinite assignments are "denoted as "M" or "F""; the output uses `*M`/`*F`.
- The column list gives `#Hom` before `#Het`; the output has `#Het` first.
- `seed` says "chosen randomly by default" but does not say the chosen value is not reported.
- `@return` says only "A dataframe".
Failure scenario: a user reading `#Het.g = 4` as 4 of 66 gametologs misreads the evidence; a user searching for `"M"` to find indefinite calls finds none.
Proposed change: correct each point above, and add the F4 assumption.

**F10 [LOW, confidence: high]: roxygen house style (DOC1, DOC2, DOC7 (proposed rule))**
`R/gl.infer.sex.r:13-18, 62-63, 66-67`: no `@family`; `seed` lacks `[default NULL]`; `@author` has a custodian but no `Author(s):` line; the example builds a plot (`plot.display = TRUE`) it does not need.
Failure scenario: `?gl.infer.sex` has no link to the other sexlinked functions; the example spends time on a plot during `R CMD check`.
Proposed change: add `@family`, the default tag, the `Author(s):` line, the `verbose` parameter text (DOC2), and `plot.display = FALSE` in the example.

## Proposed changes

1. Check the input at entry: stop if `gl_sexlinked` lacks `results.table` or the subset names for `system`, naming the system it was built for; stop if all three subsets are `NULL` (F1, F5). **Consequence: calls that now return all-`*F` or fail with "invalid 'times' argument" stop with an explanatory error instead.**
2. Return `NA` in `agreed.sex` when no preliminary assignment exists (F2). **Consequence: `agreed.sex` changes from `*F` to `NA` for individuals with no usable genotypes.**
3. Return `*?` in `agreed.sex` for a one-against-one tie (F3). **Consequence: `agreed.sex` changes from `*F` to `*?` for tied individuals; 28 individuals in LBP (`xy`).**
4. Document that the X/Z-linked and gametolog assignments assume both sexes are present (F4). Documentation only.
5. Restore the caller's random number stream on exit (F6). **Consequence: code run after `gl.infer.sex()` draws different random numbers than it does today.**
6. Add `verbose = NULL` with the standard verbosity flags and gated, coloured notes; fix the typo and indentation (F7). **Consequence: new argument (API1, API3); at the default verbosity the printed text changes; `dartr2shiny` config needs the new argument.**
7. Match rows by position, not by name (F8). No change to output for datasets with unique names.
8. Correct the roxygen documentation and house style, then run `devtools::document()` (F9, F10).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API: run. DEP: only `stats` is used, imported. PLT: no plots. FS8: returns a data frame, so no history entry is expected.
- Analysis-family checks: numerical correctness against an independent computation: the count columns were not recomputed independently; the final assignment was checked against recorded sex on LBP (344 of 345 definite assignments match). SNP vs SilicoDArT dispatch: not applicable (input is `gl.keep.sexlinked()` output, SNP only). NA handling: run (F2, F8).
- Spec: behaviour vs roxygen on LBP (`xy`, `zw`, females-only subsets) and platypus.gl (`xy`): run.
- Stochasticity: with `seed = NULL`, two runs on LBP gave identical output, because k-means converges to the same partition; the seed matters only for less separated data. Not tested further.
- Label choice in `Z.sex()`/`g.sex()` uses the individual with the largest heterozygote count rather than proportion (`R/gl.infer.sex.r:296, 354`). No failure found on LBP or platypus.gl; noted, not raised as a finding.
- Downstream callers: grep of dartR.base, dartR.popgen, dartR.captive, dartRverse and dartr2shiny: only dartr2shiny uses it (F7).
- dartR Google Group: not searched. GitHub issues for dartR.sexlinked: none filed.
- FBM path (DAT6): SKIPPED, because no FBM fixture is available and the input subsets are small.
- Environment note: `Rscript` on `PATH` is Homebrew R 4.6.1 with no packages installed; all runs used `/usr/local/bin/Rscript` (R 4.4.2).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis | consequence approved |
| 4 | approved | Luis | |
| 5 | approved | Luis | consequence approved |
| 6 | approved | Luis | consequence approved |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

Branch `review-gl.infer.sex` from `origin/dev` (179080b). Characterization test: 29 expectations pass; full suite (report, keep, filter, infer) 107 expectations, 0 failures.

- 1 applied: mismatched `system` and all-`NULL` input now stop with an error (tests "input is checked").
- 2 applied: an individual with no genotypes gets `NA` (was `*F`).
- 3 applied: ties give `*?`. On LBP `xy`, the 28 tied individuals change `*F` to `*?`; this is the only difference from 179080b on LBP `xy`, LBP `zw` and platypus.gl `xy` (all other columns identical).
- 4 applied: `@details` states the both-sexes assumption.
- 5 applied: `runif()` after the call equals `runif()` without it, with and without `seed`.
- 6 applied: `verbose` added; `verbose = 0` is silent; `verbose = 3` run end to end on LBP.
- 7 applied: with a duplicated name, X-linked assignments equal those with unique names (were all `NA`).
- 8 applied: roxygen corrected, `devtools::document()` run. An unrelated diff in `man/gl.drop.sexlinked.Rd` (existing drift between its source and man page) was reverted to keep the PR to one function.
- NEWS entry added. Follow-up outside this PR: the roxygen copy in `dartr2shiny/input_generator/dartR.sexlinked/gl.infer.sex.r` is out of date.
- PR: dartR.sexlinked#36.

```json
{
  "function": "gl.infer.sex",
  "package": "dartR.sexlinked",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "179080b",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "principle: no global side effects", "status": "approved", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS2", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group not searched"],
  "status": "pr-open",
  "pr": 36
}
```
