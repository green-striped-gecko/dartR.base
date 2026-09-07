# Review: gl.report.heterozygosity (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.report.heterozygosity
- Datasets: testset.gl (30 pops, many below n.limit = 10), platypus.gl
  (all pops above n.limit)
- Family mode: report
- Checks skipped: bootstrap CI path exercised only for mechanics, not
  distributional correctness (delegated to package boot); parallel path
  not exercised (single-CPU run); FBM path not exercised (not available:
  no FBM fixture); Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — method-coercion, negative-n.invariant and
  secondaries-history warnings all print at verbose = 0 (VRB2); the
  nboots/CI validation uses the cat(error)+stop() split; header breaks
  DOC1/DOC2/DOC7 and the `@family unmatched report` label ignores its
  matched filter.
- **Spec: FAIL** — three reproducible crashes in documented feature
  combinations: (1) `subsample.pop = TRUE` crashes whenever ANY
  population is below n.limit ("Input is logical but should be a plain
  list…" — utils.subsample.pop stores NA placeholders that rbindlist
  rejects), although n.limit is documented as a skip threshold; (2)
  `method = 'ind'` with `plot.display = FALSE` and `verbose >= 3`
  crashes ("object 'outliers' not found" — the outlier table is built
  only inside the plotting block); (3) `subsample.pop = TRUE` with
  `method = 'ind'` crashes ("object 'res_sub' not found" — referenced in
  the ind branch and the return path but computed only under 'pop').
  The core statistics are exact: Ho and He match hand computation from
  genotype counts; FIS = (uHe − Ho)/uHe per locus matches the
  documented formula.

## Findings

### F1 — subsample.pop crashes when any population is below n.limit (HIGH, confidence: certain)

Rule: spec axis (documented feature unusable). Location:
`R/utils.het.report.r` utils.subsample.pop (dependency) +
`R/gl.report.heterozygosity.r` subsample plot block.

utils.subsample.pop assigns `NA` for below-limit populations and then
runs every entry through data.table::rbindlist, which rejects the
logical NA. Any dataset with one population under n.limit (testset.gl
has many) crashes; platypus.gl (all pops >= 10) works. The n.limit
parameter is documented as "Minimum number of individuals that should
have a population to perform subsampling" — i.e. skip, not crash.

Proposed change (dependency fix, PR #248 precedent): drop below-limit
populations before stacking, keeping names aligned; in the report's
subsample plot, key colours by population name
(`colors_pops[as.character(pop)]`) instead of `rep(colors_pops, each=5)`
so skipped populations cannot desync the palette. Consequence:
subsample.pop = TRUE works on datasets with small populations,
reporting only the qualifying ones.

### F2 — method='ind' + display off + verbose >= 3 crashes (HIGH, confidence: certain)

Rule: spec axis / PLT3 (results independent of plotting). Location:
`R/gl.report.heterozygosity.r` (outliers block).

`outliers` is computed inside `if (plot.display)` (from the built
boxplot) but consumed in the verbose >= 3 report. Verified crash.
Proposed change: compute outliers from the data
(`boxplot.stats(df$Ho)$out`) regardless of plotting.

### F3 — subsample.pop + method='ind' crashes (MEDIUM, confidence: certain)

Rule: spec axis. Location: `R/gl.report.heterozygosity.r` (ind branch
print; return path).

`res_sub` exists only under method='pop', but the ind branch prints it
at verbose >= 2 and the shared return path returns
`list(res_sub, df)` whenever subsample.pop = TRUE. Verified crash.
Proposed change: subsample.pop is a pop-method feature — under
method='ind' warn (gated at verbose >= 1) that it is ignored and return
the plain dataframe.

### F4 — Warnings ungated at verbose = 0; cat(error)+stop() split (LOW, confidence: certain)

Rule: VRB2. Location: parameter-checking block.

The method-coercion warning, the negative-n.invariant warning and the
gl.filter.secondaries-history warning all print at verbose = 0
(verified); the nboots/CI check prints an error then calls stop() with
no message. Proposed change: gate the warnings at verbose >= 1; single
`stop(error(...))`.

### F5 — subsample.pop return is an unnamed list; @return silent about it (LOW, confidence: certain)

Rule: DOC5. `list(res_sub, df)` (verified unnamed) with `@return`
promising a dataframe. Proposed change: named list
(`list(subsample = res_sub, results = df)`) and document both shapes.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. `@family unmatched report` — but
gl.filter.heterozygosity is its matched filter (via method='ind');
changed to "matched report". DOC2 wording; `@author` Custodian only
(add Author(s), retaining the CP/Carlo Pacioni code credits); FLAG
SCRIPT START after SET WORKING DIRECTORY (FS4).

## Report notes (other functions / not fixed here)

- The results table prints at verbose >= 2 (not the ratified >= 1 used
  in the smaller reports); it is silent at verbose 0, which is the core
  rule — left as-is.
- gl.colors banners can leak when defaults are evaluated with plotting
  on — a gl.colors defect (already noted in the hwe review).
- The commented-out `if (ncpus > 1)` guard around the parallel set-up is
  dead scaffolding; the computed value is harmless — left as-is.
- verbose == 5 comparisons are equivalent to >= 5 (5 is the maximum) —
  not a defect.

## Coverage

Characterization baseline:
`tests/testthat/test-gl.report.heterozygosity.R` — 15 assertions: Ho/He
vs hand computation, output columns, verbose-0 silence both methods,
invisible return + untouched input, the three crashes (baseline),
platypus subsample structure (unnamed list, baseline), ungated warnings
(baseline). All 15 pass on the pre-fix code. Bootstrap mechanics not in
the baseline (runtime); exercised ad hoc in Phase C verification.

## Approval

Held for author consultation (2026-08-30). After the member heard back
from the collaborators, all six findings were approved via the
approval boxes (2026-08-31).

## Addendum (discovered during apply)

The shared plot-save block references p3 unconditionally, but both
method branches build their plots only under plot.display — so
plot.file with plot.display = FALSE crashed on the unbuilt plot (the
p3-not-found class). Below-HIGH; fixed with an existence guard and a
gated (verbose >= 1) warning that nothing was saved.

## Outcome

All six findings applied, plus the addendum. Characterization suite:
21/21 pass; flips map to F1 (subsample on testset.gl now returns
results for exactly the populations at or above n.limit), F2 (ind
results independent of plotting, outliers reported), F3 (gated
warning, plain dataframe), F4 (silent at verbose 0), F5 (named
$subsample/$results list). Reported statistics unchanged (Ho/He still
match the hand computation). End-to-end verified: pop+subsample with
display at verbose 3 (skipped pops absent from the subsample plot,
colours keyed by name), ind at verbose 3 without display, plot.file
without display (gated warning). NEWS entry added.

```json
{
  "function": "gl.report.heterozygosity",
  "package": "dartR.base",
  "family_mode": "report",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/utils.het.report.r + subsample plot block", "status": "applied"},
    {"id": "F2", "severity": "HIGH", "rules": ["spec", "PLT3"], "loc": "R/gl.report.heterozygosity.r outliers block", "status": "applied"},
    {"id": "F3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.report.heterozygosity.r ind branch/return", "status": "applied"},
    {"id": "F4", "severity": "LOW", "rules": ["VRB2"], "loc": "R/gl.report.heterozygosity.r checks", "status": "applied"},
    {"id": "F5", "severity": "LOW", "rules": ["DOC5"], "loc": "R/gl.report.heterozygosity.r return", "status": "applied"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.report.heterozygosity.r header", "status": "applied"}
  ],
  "datasets": ["testset.gl", "platypus.gl"],
  "baseline_test": "tests/testthat/test-gl.report.heterozygosity.R",
  "pr": 291
}
```
