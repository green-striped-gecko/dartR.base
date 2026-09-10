# Review: gl.report.ld.map (dartR.base)
- Family mode: report
- Date: 2026-09-10
- Reviewer: Claude (Claude Fable 5), dartr-function-review v1.0.0
- Package commit: f5e7b72 (upstream/dev; working copy verified identical by `git diff upstream/dev -- R/gl.report.ld.map.r R/utils.read.ped.r R/gl2plink.r`)
- Datasets: platypus.gl (mapped, 383 loci after callrate/monomorph filtering), testset.gl (unmapped subset), testset.gs (SilicoDArT) — dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.report.ld.map.R (new file, snapshot captured pre-review; 27 assertions, all passing)

## Verdicts

**Standards: Needs work** — the FS backbone is present and the input object
comes back untouched, but the function is not silent at `verbose = 0` (notes,
warnings and the plot all still appear), and saving a plot crashes when the
plot was not displayed.

**Spec: Needs work** — the core report is trustworthy where it reports:
independently recomputed `snpStats` R-squared values match every sampled
reported pair exactly, and locus names and positions are correctly paired.
But for signed LD statistics the function silently discards every pair with
a non-positive value — for `ld.stat = "R"` that is roughly half of all
computed pairs — and several documented claims are wrong.

## Independent verification (spec axis)

For the largest platypus population (TENTERFIELD, 41 individuals, 256 loci
after the function's own per-population MAF filter), a `SnpMatrix` was built
directly from the genlight dosages and `snpStats::ld(stats = "R.squared")`
computed independently. Eight randomly sampled reported pairs, keyed by the
reported locus names, match to within 5e-16; reported `pos_loc_a`/`pos_loc_b`
match the object's `@position` slot 8/8. The report's pairing of names,
positions and statistics is correct.

## Findings

**F1 [HIGH, confidence: high] — pairs with non-positive LD are silently discarded (DOC5 (proposed rule))**
`R/gl.report.ld.map.r:268` — `ld_columns <- ld_columns[-ld_columns$Freq < 0, ]`
keeps only rows with `Freq > 0`. The zero entries of the sparse matrix
returned by `snpStats::ld` (pairs outside the computed depth band) are
conflated with genuinely computed values that are zero or negative, so the
filter throws both away.
Failure scenario: `gl.report.ld.map(x, ld.stat = "R")` on the platypus set
returns 234 pairs, minimum +0.0029, although 16,771 of the 32,640 computed
all-pairs R values in TENTERFIELD alone are negative. Every signed statistic
("R", "Covar", "Q", "LLR" excepted as it is non-negative) loses roughly the
half of its distribution below zero; even for R-squared, exact zeros are
dropped. Boxplots, histograms and the pops-in-LD plot are all computed from
the truncated distribution.
Proposed change: select pairs by membership of the computed band (index
offset within the `snpStats::ld` depth, or explicit tracking of computed
entries) instead of by the sign of the value.

**F2 [HIGH, confidence: high] — plot.file with plot.display = FALSE crashes after the full computation (PLT3)**
`R/gl.report.ld.map.r:322-408` — `p1`, `p2`, `p3` and `p4` are only built
inside `if (plot.display)`; the save branch `if (!is.null(plot.file))`
references `p4` unconditionally.
Failure scenario: `gl.report.ld.map(x, plot.display = FALSE,
plot.file = "ld")` errors with "object 'p4' not found" after the whole LD
computation has run; the returned data frame is lost.
Proposed change: build the plots whenever `plot.display` or `plot.file`
requires them; display and save independently.

**F3 [MEDIUM, confidence: high] — not silent at verbose = 0: ungated messages and plot display (VRB5, VRB3, VRB4 (proposed rule))**
`R/gl.report.ld.map.r:144-150` (no-chromosome note), `:180-182`
(skipped-population warning), `:322`/`:404` (plot gating).
Failure scenario: on an unmapped testset.gl subset at `verbose = 0` the
function prints 57 lines (the chromosome note plus 27 skipped-population
warnings), and with default `plot.display = TRUE` the plot renders
(Rplots.pdf created in a script session). There is no
`if (verbose == 0) plot.display <- FALSE` line.
Proposed change: gate the skipped-population warning at `verbose >= 1`
(it means the results omit populations — VRB4), the chromosome note at
`verbose >= 2`, and force `plot.display <- FALSE` at `verbose == 0`.

**F4 [MEDIUM, confidence: high] — documentation contradicts behaviour (DOC5 (proposed rule), DOC1)**
`R/gl.report.ld.map.r:9-14` says `ld.max.pairwise` "should be set as NULL
(the default)" while the signature default is 1000000 (the unmapped case is
actually detected from the empty `@chromosome` slot, so the parameter text is
doubly misleading); `:40-41` describes `plot.display` as "histograms of base
composition" (copied from another function); `@family graphics` although this
is the matched report of `gl.filter.ld` (its filter twin is family "matched
filter"); `:34-35` documents `ind.limit` as a "minimum number of individuals"
but `:179` skips populations with exactly `ind.limit` individuals (`<=`).
Failure scenario: a user with mapped SNPs who follows the description and
passes `ld.max.pairwise = NULL` silently gets fabricated positions 1..n and
all-pairs LD instead of map-distance LD.
Proposed change: correct the four documentation points (behaviour unchanged);
`@family matched report`.

**F5 [MEDIUM, confidence: high] — loci sharing a map position are silently excluded (VRB4 (proposed rule))**
`R/gl.report.ld.map.r:233-238` — within each chromosome, loci with a
duplicated `loc_bp` are removed from the analysis with no message at any
verbosity.
Failure scenario: DArT data mapped to a reference commonly carries multiple
tags at the same position; all but the first are absent from the report and
from any downstream `gl.filter.ld` decision, and nothing tells the user.
Proposed change: count the exclusions and report them at `verbose >= 1`.

**F6 [LOW, confidence: high] — SilicoDArT admitted at entry, rejected only mid-loop (DAT7 (proposed rule))**
`R/gl.report.ld.map.r:113` — the datatype check uses the default `accept`,
so presence/absence data passes the entry check and dies later inside an
internal call ("found SilicoDArT expecting SNP" raised from within the
per-population loop).
Failure scenario: the user gets the right outcome (an error) by accident and
after work has begun; a change to the internal functions would silently
admit silico data to an R-squared calculation on 0/1 calls.
Proposed change: `utils.check.datatype(x, accept = "SNP", verbose = verbose)`.

**F7 [LOW, confidence: high] — dead and misleading code (STY1, DEP1)**
`R/gl.report.ld.map.r:128-136` guards package `fields`, which the function
never uses; `:336-342` contains a self-assignment branch
(`boxplot.colors <- boxplot.colors`) and a duplicated `is()` test; `:259`
overwrites the `ld.max.pairwise` parameter inside the loop; `:388-401` and
`:419-436` are large commented-out blocks.
Failure scenario: none at run time; the `fields` guard makes users install a
package the function does not need.
Proposed change: remove the `fields` guard, the dead branches and the
commented-out blocks.

## Proposed changes

1. Keep computed zero and negative LD values: select pairs by the computed
   band rather than by value sign (F1).
   **Consequence: numerical output changes — row counts increase for every
   statistic whose computed values include zeros or negatives; for signed
   statistics the report gains the entire negative half of the
   distribution.**
2. Build plots independently of display so `plot.file` works with
   `plot.display = FALSE` (F2).
3. Verbosity gating: skipped-population warning at `verbose >= 1`,
   chromosome note at `verbose >= 2`, `plot.display` forced off at
   `verbose == 0` (F3).
4. Documentation corrections: `ld.max.pairwise` semantics, `plot.display`
   text, `@family matched report`, `ind.limit` boundary wording (F4).
5. Report the count of duplicate-position loci excluded per chromosome at
   `verbose >= 1` (F5).
6. Restrict the datatype check to SNP data (F6).
7. Remove the unused `fields` guard, dead branches, the in-loop parameter
   overwrite, and commented-out blocks (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen on platypus.gl (mapped) and testset.gl
  (unmapped) — run; independent recomputation of reported statistics — run
- Signed-statistic path (`ld.stat = "R"`) — run
- GitHub issues: #210 ("gl.filter.ld breaks on subsets") checked — NOT
  reproducible on dev f5e7b72 (50-locus subset of the example runs clean end
  to end: report 1056x11, filter to 46 loci); candidate for closing
- FBM path (DAT6): SKIPPED — the function routes through gl2plink/file
  round-trip which densifies regardless; no FBM fixture exercised
- Google Group search: SKIPPED — not queried this session

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | rejected | Arthur | keep discarding pairs with LD <= 0; the truncation is now stated in the documentation (folded into change 4) |
| 2 | approved | Arthur | |
| 3 | approved | Arthur | |
| 4 | approved | Arthur | includes documenting the change-1 truncation |
| 5 | approved | Arthur | |
| 6 | approved | Arthur | |
| 7 | approved | Arthur | |

Cross-package caller grep (API3): the only consumer of this function's
output across the local dartR.* clones is `dartR.popgen::gl.ld.distance`
(takes the report data frame). The approved set leaves the returned data
frame unchanged (change 1 rejected), so no caller is affected. All clear.

## Outcome

(pending Phase C)

```json
{
  "function": "gl.report.ld.map",
  "package": "dartR.base",
  "family": "report",
  "skill_version": "1.0.0",
  "commit": "f5e7b72",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "proposed", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "proposed", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "proposed", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "proposed", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DAT7", "status": "proposed", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "proposed", "change": 7}
  ],
  "coverage_skipped": ["DAT6: no FBM path exercised", "Google Group: not queried"],
  "status": "awaiting-approval",
  "pr": null
}
```
