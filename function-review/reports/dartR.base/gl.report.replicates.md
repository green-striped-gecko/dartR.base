# Review: gl.report.replicates (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: 2739215 (dev, synced with upstream/dev);
  Branch: review-gl.report.replicates; Datasets: platypus.gl subset
  (30 x 500) with two exact duplicate individuals injected
- Family mode: report (pairwise identity screening; input untouched
  verified, no history append, correct)
- Custodian: Luis Mijangos
- Checks skipped: RcppParallel kernel not reviewed line-by-line beyond
  verifying its outputs against a direct R computation; Google Group
  not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — no datatype/genlight check at all
  (as.matrix(x) is called directly); Rcpp and RcppParallel are
  compiled/used with no requireNamespace check, and RcppParallel is
  absent even from DESCRIPTION Suggests; no verbose gate on the plot
  (renders at verbose 0, verified); the result is returned visibly
  where house reports return invisibly; @family "report functions"
  instead of "matched report"; markdown-style "##" headings in
  @details render literally (markdown mode is off in this package).
- **Spec: FAIL** — three verified defects. (1) The pair table carries
  BOTH orderings of every pair (3 real pairs → 6 rows), and the drop
  rule `ifelse(ind1_miss > ind2_miss, ind1, ind2)` selects the
  opposite member in each ordering whenever missingness ties — which
  is precisely the exact-duplicate case — so BOTH replicates land in
  ind.list.drop (verified: T27 and T27.1, T35 and T35.1). The
  histograms also double-count every pair. (2) When no pairs pass the
  thresholds the function returns a bare character message instead of
  the documented 3-element list — feeding that to
  gl.filter.replicates crashes (verified). (3) ind.list.rep uses
  `>= perc_geno` while table.rep/ind.list.drop use `> perc_geno` —
  inconsistent outputs at the boundary. The pairwise statistics
  themselves are verified correct (nloc and perc match a direct R
  computation for spot-checked pairs).

## Findings

### F1 — Tied missingness drops both replicates; doubled pair table (HIGH, confidence: certain)

Rule: spec axis. Location: pair table and drop rule.

Proposed change: reduce the table to one row per unordered pair, and
break missingness ties deterministically (drop ind2, the second-named
individual). Consequences: ind.list.drop retains one member of each
tied pair (the bug fix); table.rep halves to one row per pair; the
histograms count each pair once.

### F2 — No-pairs path returns a string (MEDIUM, confidence: certain)

Rule: spec/API. Location: early return.

Proposed change: return the documented list with an empty table.rep
(correct columns), empty ind.list.drop/ind.list.rep, message gated at
verbose >= 1; the Completed flag runs. Downstream,
gl.filter.replicates can then no-op gracefully.

### F3 — Threshold inconsistency; missing checks and gates (MEDIUM, confidence: certain)

Rule: spec/DEP/VRB. ind.list.rep aligned to the strict `>` comparisons
used everywhere else; utils.check.datatype added; Rcpp and
RcppParallel checked with stop(error(...)) (and RcppParallel added to
DESCRIPTION Suggests); `verbose == 0 → plot.out <- FALSE` guard added.

### F4 — Return visibility, results display, and header (LOW, confidence: certain)

Rule: VRB/DOC. Return made invisible (house standard for reports);
table.rep printed at verbose >= 3; @family corrected to "matched
report"; "##" headings in @details converted to \strong{}; curly
quotes; typo "asses".

## Report notes (other functions / not fixed here)

- Parameter names (loc_threshold, perc_geno, plot.out, plot_theme,
  plot_colors) use legacy underscore/dot naming rather than the house
  plot.display/plot.theme/plot.colors set; renaming is an API break
  across the pair and is left to the custodian.
- The C++ kernel compiles on first call per session (~7 s observed),
  as documented.

## Coverage

`tests/testthat/test-gl.report.replicates.R` — 16 assertions: pair
statistics verified against a direct computation, doubled-table and
tied-drop baselines, no-pairs string baseline, visible-return
baseline, verbose-0 plot render baseline, input untouched. All 16
pass on the pre-fix code (Rcpp/RcppParallel skipped if unavailable).

## Approval

All four findings approved via the approval boxes (2026-08-31); F1
approved with its output-shape consequences (one row per pair, tied
pairs keep one member, histogram counts halve).

## Outcome

All four findings applied. In the applied F1 orientation ind1 is the
individual EARLIER in the object and a missingness tie drops ind2 (the
later individual — for suffixed duplicates, the ".1" copy). One
additional latent defect fixed within F3's scope: the NaN diagonal of
the keep_pairs matrix injected NA entries into ind.list.rep for every
individual; it is now treated as no-match, and ind.list.rep names
exactly the individuals with replicate partners. Characterization
suite: 18/18 pass; flips map to F1 (3 rows, drop list T27.1/T35.1/T3),
F2 (documented empty structure), F3 (no page at verbose 0), F4
(invisible return). Pair statistics unchanged (still match the direct
computation). NEWS entry added; RcppParallel added to DESCRIPTION
Suggests.

```json
{"function": "gl.report.replicates", "package": "dartR.base", "family_mode": "report",
 "commit": "2739215", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.report.replicates.r pair table", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "API"], "loc": "R/gl.report.replicates.r early return", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["spec", "DEP", "VRB"], "loc": "R/gl.report.replicates.r checks", "status": "applied"},
  {"id": "F4", "severity": "LOW", "rules": ["VRB", "DOC"], "loc": "R/gl.report.replicates.r", "status": "applied"}],
 "datasets": ["platypus.gl (constructed duplicates)"],
 "baseline_test": "tests/testthat/test-gl.report.replicates.R", "pr": null}
```
