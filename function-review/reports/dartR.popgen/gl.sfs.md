# Review: gl.sfs (dartR.popgen)
- Family mode: analysis
- Custodian: Bernd Gruber & Carlo Pacioni (STY5: this report is the
  discussion record; changes approved by Luis)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 461af97 (origin/dev)
- Datasets: possums.gl, two populations (4 + 3 individuals, 193 loci after
  full call rate and monomorph filters); platypus.gl SEVERN_ABOVE (23
  individuals, 1,000 loci, 237 with missing calls); testset.gs[1:20, 1:50]
- Baseline: tests/testthat/test-gl.sfs.R (9 tests, 17 expectations; runs
  anywhere). Reference values are independent counts from `as.matrix(x)`

## Verdict

**Standards: Needs work** — SilicoDArT data is accepted, `plot.file` fails
whenever no bar plot was drawn, and the warnings ignore `verbose` or stop
with an empty message.
**Spec: Rework** — the single-population spectrum is right for complete
data, but loci with missing calls are placed in the wrong class (a
monomorphic locus can become polymorphic), the folded multi-population
spectrum folds each population on its own minor allele, and `minbinsize`
on a multi-population spectrum removes polymorphisms private to one
population.

What works well: single-population folded and unfolded spectra, and the
unfolded two-population spectrum, match independent counts cell by cell on
complete data; the function feeds `gl.run.stairway2` and `gl.run.epos`
correctly under those conditions.

## Findings

**F1 [HIGH, confidence: high] — missing calls put loci in the wrong class (DOC5)**
`R/gl.sfs.r:103-111,127` — allele counts are summed with `na.rm = TRUE` and
classed against the full sample size 2n. A locus with k individuals
missing is scored out of 2(n - k) sequences but folded around n.
Failure scenario (baseline test 5): platypus.gl SEVERN_ABOVE, 237 of 1,000
loci with missing calls. Among the called individuals 487 loci are
polymorphic; `gl.sfs` reports 525. A locus fixed for the alternative
allele in all called individuals, with one individual missing, is folded
into class 2 as if it were polymorphic (38 such loci here). The warning is
shown only at `verbose >= 2`, so `gl.run.stairway2` and `gl.run.epos`
(which call it at lower verbosity) never show it.
Proposed change: build the spectrum only from loci scored in every
individual (every individual of every population for a multi-population
spectrum) and report how many loci were dropped at `verbose >= 1`.
Alternatives, for the member to choose: (b) keep all loci but raise the
warning to `verbose >= 1` with the number of affected loci; (c) project
each locus down to a common sample size (hypergeometric projection, as in
dadi/easySFS), which keeps most loci but is a larger change.
**Consequence (recommended option): the spectrum changes for every
dataset with missing calls, and so do the Stairway Plot 2 and EPOS
results built on it.**

**F2 [HIGH, confidence: high] — folded multi-population spectrum folds each population separately (DOC5)**
`R/gl.sfs.r:122-133` — each population's count is folded around its own
sample size.
Failure scenario (baseline test 3): a locus with 6 of 8 alternative
alleles in population A and 1 of 6 in population B is recorded as (2, 1)
instead of (6, 1): the minor allele of the combined sample is the
reference allele in A and the alternative allele in B. 27 of 63 cells
differ from folding on the minor allele of all populations combined (the
folded joint spectrum used by fastsimcoal2, which the documentation
names).
Proposed change: decide the minor allele from the combined count across
all populations and fold every population on that allele.
**Consequence: folded multi-population spectra change.**

**F3 [HIGH, confidence: high] — `minbinsize` on a multi-population spectrum removes private polymorphisms (DOC5)**
`R/gl.sfs.r:158-169` — the first `minbinsize` classes are removed in every
dimension, so every site with a count below `minbinsize` in any one
population is dropped.
Failure scenario (baseline test 4): unfolded, `minbinsize = 1`: 10 of 193
polymorphic sites (all polymorphic in one population and absent from the
other) are removed along with the monomorphic cell.
Proposed change: keep the full array and set to 0 the cells whose total
count over all populations is below `minbinsize` (for `minbinsize = 1`,
only the all-zero cell); replaces the `eval(parse())` code.
**Consequence: multi-population spectra with `minbinsize > 0` keep their
full dimensions and include private polymorphisms.**

**F4 [MEDIUM, confidence: high] — `plot.file` fails when no bar plot was drawn (PLT3, PLT1)**
`R/gl.sfs.r:174-199` — `gp` is only created inside `if (plot.out)` and only
for a vector.
Failure scenario (baseline test 6): `plot.out = FALSE, plot.file = "sfs"`
stops with "object 'gp' not found" and the spectrum is not returned; the
same for any multi-population spectrum. The bars are placed at 1, 2, ...
whatever the class, so after `minbinsize = 2` the bar for class 2 is
labelled 1; the plot uses no theme.
Proposed change: build the plot for a vector spectrum whenever it is
needed (display or save), with class numbers on the x axis and
`plot.theme`; for an array, skip the save with a message.
**Consequence: new argument `plot.theme` (default `theme_dartR()`).**

**F5 [MEDIUM, confidence: high] — SilicoDArT data is accepted (FS4)**
`R/gl.sfs.r:67` — `datatype` is never used.
Failure scenario (baseline test 7): `testset.gs` returns a 21-class
"spectrum" of presence/absence scores.
Proposed change: stop on SilicoDArT input.
**Consequence: SilicoDArT input errors.**

**F6 [LOW, confidence: high] — messages ignore `verbose` or are empty (VRB2, VRB3, FS5)**
`R/gl.sfs.r:80-96,185` — "No population definition provided" prints at
every verbosity; the dimension check prints with `cat()` and calls
`stop()` with no message; "more than 2 dimensions" is printed for a
two-population spectrum.
Failure scenario (baseline tests 8, 9): at `verbose = 0` the population
message prints; `gl.sfs(possums.gl)` stops with an empty error.
Proposed change: gate the population message at `verbose >= 2`; put the
dimension message in `stop(error(...))`; correct the plot message.

**F7 [LOW, confidence: high] — the genotype matrix is built twice (STY2, DAT6 (proposed rule))**
`R/gl.sfs.r:71,103,127` — `as.matrix(x)` once to count missing values and
again for the counts.
Failure scenario: on a large dataset the full matrix is materialised
twice.
Proposed change: build it once and reuse it.

**F8 [LOW, confidence: high] — documentation gaps (DOC1, DOC5, DOC6 (proposed rule), DOC7 (proposed rule))**
`R/gl.sfs.r:1-42` — no `@name`/`@title`/`@family`; `@author` has no
Author(s) part; non-ASCII "Sanchez" accent in `@references`; `plot.dir`
documented as the working directory (it is the dartR default, usually
`tempdir()`); missing-data handling and the folding rule for several
populations are not described; "through an error".
Failure scenario: a user reads that missing data is tolerated and that
the folded joint spectrum matches fastsimcoal2.
Proposed change: fix the items above and describe missing-data handling
and joint folding.

## Proposed changes

1. Build the spectrum from loci scored in every individual and report the
   number dropped (F1). Alternatives: warn only, or projection.
   **Consequence: the spectrum changes for every dataset with missing
   calls, and so do Stairway Plot 2 and EPOS results built on it.**
2. Fold a multi-population spectrum on the minor allele of all
   populations combined (F2).
   **Consequence: folded multi-population spectra change.**
3. `minbinsize` on a multi-population spectrum zeroes the cells whose
   total count is below it, keeping the full array (F3).
   **Consequence: multi-population spectra with `minbinsize > 0` keep
   their dimensions and include private polymorphisms.**
4. Plot built whenever needed, class numbers on the x axis, new
   `plot.theme`, arrays skip the save with a message (F4).
   **Consequence: new argument `plot.theme`.**
5. Reject SilicoDArT input (F5).
   **Consequence: SilicoDArT input errors.**
6. Messages: gate, fix the empty stop and the wrong plot message (F6).
7. Build the genotype matrix once (F7).
8. Documentation (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: single-population (folded, unfolded, `minbinsize`) and
  two-population (folded, unfolded, `minbinsize`) spectra against
  independent counts; missing data on platypus.gl; `plot.file`;
  SilicoDArT; no population; too many dimensions — run
- Three or more populations: not run (same code path as two)
- Callers (API3): `gl.run.stairway2` and `gl.run.epos` (single population,
  folded and unfolded); dartR Shiny (`minbinsize`, `folded`,
  `singlepop`); dartR.sim `gl.sim.Neconst` (example only) — run
- DAT6 (FBM): SKIPPED — no FBM fixture
- fastsimcoal2 file format: not checked; the function returns an R array,
  not an `.obs` file
- Google Group / GitHub issues: not searched in this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | option: drop loci with missing calls (recommended) |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

Branch `review-sfs` (from origin/dev 461af97). Evidence:
`tests/testthat/test-gl.sfs.R`, 9 tests, 24 expectations, all passing (no
external software). Every changed assertion is tagged with its approved
change; reference values are independent counts from the genotype matrix.

- Snapshot of pre-change outputs (single population folded, unfolded,
  `minbinsize = 2`; two populations unfolded with `minbinsize` 0; complete
  data): identical after the change. Differences only where approved:
  folded two-population spectrum (change 2) and bandicoot.gl, 873 loci
  with missing calls now excluded (change 1).
- 1 applied: loci with any missing call are excluded before counting;
  "237 of 1000 loci have missing calls and were excluded" at
  `verbose >= 1`; the spectrum equals an independent count on the 763
  complete loci (test 5).
- 2 applied: counts flipped to the reference allele when the alternative
  allele is the more frequent over all individuals; equals the
  independent joint folding cell by cell (test 3).
- 3 applied: cells with total class below `minbinsize` set to 0; array
  keeps its dimensions; the 10 private polymorphisms are kept (test 4).
  Replaces the `eval(parse())` code.
- 4 applied: plot built for display or save, x axis = class number,
  `plot.theme` argument added after `plot.dir` (positional callers of
  existing arguments unaffected); `plot.file` works with `plot.out =
  FALSE`; arrays skip the save with a warning (test 6).
- 5 applied: SilicoDArT errors (test 7).
- 6 applied: population message at `verbose >= 2`; dimension error in
  `stop(error())` (tests 8, 9).
- 7 applied: `as.matrix(x)` once; the per-population fill loop replaced by
  one `table()` of factors.
- 8 applied: roxygen rewritten; `devtools::document()` also rewrote seven
  unrelated Rd files (stale family links on dev), reverted; committed Rd:
  `gl.sfs.Rd` plus the new "demographic history" family links in
  `gl.run.epos.Rd` and `gl.run.stairway2.Rd`.
- Downstream: `test-gl.run.stairway2.R` and `test-gl.run.epos.R` pass with
  their binaries (38 and 37 expectations); their fixtures have no missing
  calls.
- `devtools::check()`: 0 errors; install WARNING (packages built under
  R 4.4.3) and 3 NOTEs present before this change.
- NEWS entry added. Callers: gl.run.stairway2, gl.run.epos, dartR Shiny,
  dartR.sim example.
- PR: #106

## Machine block

```json
{
  "function": "gl.sfs",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "461af97",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["3+ populations not run", "DAT6: no FBM fixture", "fastsimcoal2 file format not checked", "Google Group / issues not searched"],
  "status": "pr-open",
  "pr": 106
}
```
