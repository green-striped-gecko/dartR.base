# Review: gl.select.panel and gl.check.panel (dartR.popgen)

Reviewed as a pair: `gl.check.panel` compares a panel made by
`gl.select.panel` with the full data.

- Family mode: modify (`gl.select.panel` returns a locus subset) and
  analysis (`gl.check.panel`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 9ea1286 (origin/dev, reviewed state)
- Datasets: possums.gl (300 x 200, 10 populations, locus names X1..X200,
  no SNP metadata; the roxygen examples), bandicoot.gl (96 x 1,000, 5
  populations, DArT names and metadata), testset.gs; NeEstimator V2
  (`Ne2-1M`, macOS) for the Ne check
- Baseline: tests/testthat/test-gl.select.panel.R (11 tests, both
  functions, snapshot captured pre-review; defects marked
  `BASELINE (F<n>)`)

**Standards: Needs work** — neither function validates its arguments;
`gl.check.panel` has no verbosity, start/end flags or datatype check;
`gl.select.panel` prints at `verbose = 0` and documents plot arguments it
never uses.
**Spec: Rework** — three of the nine selection methods do not work on
the package's own example data. `dapc` silently returns every locus, and
`pahigh` and `monopop` stop with errors. `gl.check.panel` fails on the
parameter name its help page gives, and accepts a "full" data set that
holds different individuals.

What works well: `random`, `stratified`, `hafall`, `hafpop`, `pic` and
`picdart` return `nl` loci with metadata kept in step. On bandicoot.gl,
`pahigh` picks only private-allele loci and `dapc` returns 30 loci.
`gl.check.panel` computes Fst, He, Ho, Fis, allelic richness and Ne
correctly for panels whose individuals match the full data; for a random
100-locus panel of possums.gl the correlations are 0.97-0.99.

## Findings

**F1 [BLOCKER, confidence: high] — `dapc` matches loci by name, finds
none on possums.gl and silently returns every locus (DOC5)**
`R/gl.select.panel.R:108-119, 274` — the loci are taken from the row
names of `dapc()$var.contr` and split at the first dot. For possums.gl
those row names are `1, 2, 3, …` while the locus names are `X1, X2, X3,
…`, so `gl.keep.loc()` finds none of the 50 names, prints 50 warnings
and "no loci listed to keep! Genlight object returned unchanged", and the
function returns all 200 loci. Names that contain a dot would be cut at
the dot on any data set.
Failure scenario: the roxygen example
`gl.select.panel(possums.gl, method = "dapc", nl = 5)` returns the full
200-locus object as the "panel". A user who does not count the loci
designs a panel from everything.
Proposed change: take the loci by position (the rows of `var.contr`
follow the locus order of the object passed to `dapc`) (change 1).

**F2 [HIGH, confidence: high] — `pahigh` stops on SNP data without SNP
metadata (DAT5)**
`R/gl.select.panel.R:139-148` — allele frequencies of each population
subset come from `gl.alf()`. On a one-population subset of the
private-allele loci, all genotypes are often 0 or 1, and with no SNP
metadata `utils.check.datatype()` in dartR.base then stops: "object has
ploidy 2 (SNP) but all non-missing genotypes are scored 0 or 1 and no SNP
metadata … is present". `[1:nl2]` also yields NA names when a population
pair has fewer private alleles than `nl2`.
Failure scenario: possums.gl, `method = "pahigh"`, any `nl`: error.
Simulated or VCF-derived genlights have no SNP metadata either.
Proposed change: compute the per-population allele frequencies directly
from the genotype matrix (`colMeans(as.matrix(p), na.rm = TRUE) / 2`),
and keep at most the loci available (`head()`) (change 2).

**F3 [HIGH, confidence: high] — `monopop` fails whenever a population has
fewer monomorphic loci than its share (DOC5)**
`R/gl.select.panel.R:169-175` — `sample(dl, nl2)` with
`nl2 = ceiling(nl / nPop)`. Loci with no calls in a population give NA in
`index` and an NA name.
Failure scenario: possums.gl, `nl = 50` (5 per population): "cannot take
a sample larger than the population when 'replace = FALSE'".
Proposed change: take `min(nl2, available)` per population, drop NA, and
let `exact` top up as for the other methods (change 3).

**F4 [MEDIUM, confidence: high] — the returned object has its individuals
reordered by population (DAT2)**
`R/gl.select.panel.R:97` — `x <- x[order(pop(x)), ]` before selection,
and that object is returned. A panel is a locus subset; individual order
should not change.
Failure scenario: a genlight in field-collection order comes back sorted
by population; any vector aligned to the input individuals (phenotypes,
coordinates kept outside `ind.metrics`) is now misaligned. The order is
not needed by `seppop()`.
Proposed change: select on the input order (change 4).
**Consequence: the returned object keeps the input order of individuals
(today sorted by population).**

**F5 [MEDIUM, confidence: high] — no argument checks (FS4, FS5)**
`R/gl.select.panel.R:55-99`.
- An unknown `method` (e.g. "Random") fails with "object 'selloc' not
  found".
- `nl` above `nLoc(x)` gives "cannot take a sample larger…" (`random`)
  or "Subsetting resulted in zero loci" (`hafall`).
- An object without populations fails with "argument 1 is not a
  vector".
- SilicoDArT input passes the datatype check but fails inside `gl.alf()`
  for most methods.
Failure scenario: the errors name no argument.
Proposed change: `match.arg()` on the nine methods; `nl` a whole number
from 1 to `nLoc(x)`; populations required for the methods that use them
(`dapc`, `pahigh`, `monopop`, `stratified`, `hafpop`); SNP data only
(change 5).

**F6 [LOW, confidence: high] — output at `verbose = 0`; outdated start
flag (VRB1, VRB2, FS3)**
`R/gl.select.panel.R:75-79, 189, 224, 262` — the inner `gl.alf()` calls
run at the default verbosity, so `method = "hafpop"` prints ten
"Processing genlight object with SNP data" lines at `verbose = 0`. The
"not enough loci" message uses `report()` instead of `warn()`.
`utils.flag.start()` still gets `build = "Jody"`.
Proposed change: pass `verbose = 0` to inner calls; `warn()`; drop
`build` (change 6).

**F7 [LOW, confidence: high] — `plot.out`, `plot.file`, `plot.dir` are
documented but no plot exists (DOC5; API2, proposed rule)**
`R/gl.select.panel.R:21-23, 60-62, 70` — no plot is made or saved.
Failure scenario: `plot.file = "panel"` saves nothing, with no message.
Proposed change: remove the three arguments (change 7).
**Consequence: calls that pass `plot.out`, `plot.file` or `plot.dir` stop
with "unused argument"; an existing dartr2shiny app built from the old
header needs regenerating.**

**F8 [LOW, confidence: high] — `gl.select.panel` roxygen gaps (DOC1,
DOC5, DOC7 (proposed rule))**
`R/gl.select.panel.R:1-53`.
- No `@name`, `@family` or `@author`/custodian.
- `pic` and `picdart` are missing from `@param method`.
- `hafall`/`hafpop` are described as "highest allele frequencies", but
  they rank by minor allele frequency: the loci closest to 0.5 (on
  possums.gl the top 20 all have MAF >= 0.487).
- `exact` does not say that surplus loci are dropped at random, which
  discards the ranking for the per-population methods.
- The example comment says "Select 20 loci" with `nl = 50`.
Proposed change: rewrite the header (change 8).

**F9 [MEDIUM, confidence: high] — `gl.check.panel`: the documented
parameter "Nall" and any unknown value fail with "object 'res' not
found" (DOC5, FS5)**
`R/gl.check.panel.r:8, 81` — the help page lists "Nall"; the code tests
`"Na"`. Any other value falls through all branches.
Failure scenario: `gl.check.panel(p, possums.gl, parameter = "Nall")`
errors.
Proposed change: `match.arg()` on Fst, He, Ho, Fis, Na, Ne, accepting
"Nall" as an alias of "Na"; fix the help page (change 9).

**F10 [MEDIUM, confidence: high] — `gl.check.panel` checks that the two
objects hold the same individuals only by population label (DAT2)**
`R/gl.check.panel.r:42-47` — both objects are sorted by population and
the population labels compared position by position. A different data
set with the same population sizes passes; a panel with fewer
individuals fails with "level sets of factors are different".
Failure scenario: comparing a panel against the wrong full data set
(same populations and sizes, different animals) runs without complaint.
Proposed change: require the same individual names in both, and align
the panel to `xorig` by name (change 10).

**F11 [LOW, confidence: high] — `gl.check.panel` structure, plotting and
documentation (FS2, FS3, FS4, FS9, PLT2 (proposed rule), DOC1, DOC5)**
`R/gl.check.panel.r:1-139`.
- There is no `gl.check.verbosity()`, start/end flag or datatype check.
- `print(gg)` runs even with `plot.out = FALSE`, and `plot.file` and
  `plot.dir` are ignored.
- `@return` promises "the result of the linear regression" but the
  function returns a data frame of original vs panel values.
- The Ne value is taken by row position (`[6]`) rather than by its label
  "Estimated Ne^".
- There are no `@name`, `@family` or `@author`.
Failure scenario: `plot.out = FALSE` still draws; `plot.file = "fst"`
saves nothing.
Proposed change: standard structure; honour `plot.out`, `plot.file` and
`plot.dir` via `utils.plot.save()`; take Ne by label; rewrite the header
(change 11).

**F12 [INFO] — `gl.check.panel` runs `gl.LDNe()` with
`singleton.rm = FALSE` for the panel only**
`R/gl.check.panel.r:119-121` — the compared column ("Frequency 1",
critical 0.05) is unaffected (checked on four possums populations: panel
Ne 12.5-15.7 against full data 10.7-16.2), so no change is proposed.

## Proposed changes

1. `dapc`: select loci by position (F1).
   **Consequence: on data sets where `dapc` silently returned every locus
   (possums.gl), it returns an `nl`-locus panel.**
2. `pahigh`: allele frequencies from the genotype matrix; at most the
   loci available (F2).
3. `monopop`: take what is available per population; `exact` tops up
   (F3).
4. Keep the input order of individuals (F4).
   **Consequence: the returned object's individuals are in input order
   instead of sorted by population.**
5. Argument checks for `method`, `nl`, populations and datatype (F5).
6. Silence inner calls at `verbose = 0`; `warn()`; drop `build` (F6).
7. Remove the unused plot arguments from `gl.select.panel` (F7).
   **Consequence: calls passing them stop with "unused argument";
   dartr2shiny apps built from the old header need regenerating.**
8. `gl.select.panel` roxygen rewrite (F8). Docs only.
9. `gl.check.panel`: validate `parameter`, accept "Nall" (F9).
10. `gl.check.panel`: match individuals by name (F10).
    **Consequence: a full data set with different individuals now stops
    with an error; a panel with a subset of individuals is compared on the
    shared individuals of `xorig` only if they match by name — otherwise
    it errors as today, with a clear message.**
11. `gl.check.panel`: standard structure, working `plot.out` /
    `plot.file` / `plot.dir`, Ne by label, header rewrite (F11).

Callers: dartr2shiny generator copies of both functions (take the
roxygen headers); no `dartR.*` sibling calls them.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run for both.
- Spec: all nine methods on possums.gl; `dapc`/`pahigh` on bandicoot.gl;
  every `parameter` of `gl.check.panel` on a random possums panel,
  including Ne with NeEstimator — run.
- `pahigh` correctness on bandicoot.gl: every selected locus is a private
  allele in some population pair — run. Pair order of
  `gl.report.pa()` matches `combn()` — checked.
- `gl.check.panel` He equals `gl.report.heterozygosity()` of the full
  data — run.
- dartR Google Group: not searched for these functions. GitHub issues not
  searched.
- FBM path (DAT6): SKIPPED — no FBM fixture.
- Randomness: `random`, `monopop`, `stratified` and the `exact` top-up use
  R's RNG; reproducible with `set.seed()`, not documented (folded into
  change 8).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved: dapc returns nl loci |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis | consequence approved: input individual order kept |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis | consequence approved: plot arguments removed |
| 8 | approved | Luis |  |
| 9 | approved | Luis |  |
| 10 | approved | Luis | consequence approved: mismatched individuals now error |
| 11 | approved | Luis |  |

## Outcome

- Changes 1-11 applied (branch `review-panel`): `R/gl.select.panel.R`,
  `R/gl.check.panel.r`; two Rd files regenerated; NEWS entry added. New
  `@family panel selection` (the two functions); custodian Bernd Gruber,
  who wrote both (git history), since neither file named one.
- Snapshot diffs against the pre-review baseline: 10, all mapped. dapc
  returns 50 loci on possums.gl (change 1); pahigh and monopop no longer
  error (changes 2, 3); input individual order kept (change 4, two
  expectations); bad `method` and `nl` give the new messages (change 5, two
  expectations); verbose = 0 silent (change 6); "Nall" works (change 9);
  different individuals rejected (change 10). The untagged tests (ranking
  methods, bandicoot snapshots, check.panel He and Fst) passed unchanged.
- Old vs new code on bandicoot.gl with the same seed: identical loci for
  dapc, pahigh, hafall, hafpop, pic, stratified and random.
- Tests rewritten for the approved behaviour: 15 tests, 49 expectations,
  all pass with NeEstimator present (`NEEST_DIR`); the Ne test is skipped
  without it.
- possums.gl at `verbose = 3`: dapc panel of 50 loci; its pairwise Fst
  correlates 0.977 with the full data (45 population pairs).
- `devtools::check()`: 0 errors, 1 warning and 2 notes, all present before
  this change.
- PR: dartR.popgen#113 (commit f6b9f2a, branch `review-panel`).

## Machine block

```json
{
  "function": "gl.select.panel",
  "package": "dartR.popgen",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "9ea1286",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "API2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 9},
    {"id": "F10", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 10},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 11},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "DOC5", "status": "no_change", "change": null}
  ],
  "coverage_skipped": ["forum and GitHub issues not searched", "DAT6: no FBM fixture"],
  "status": "done",
  "pr": 113
}
```
