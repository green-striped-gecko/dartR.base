# Review: gl.report.pa (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-gl.report.pa.
- Datasets: testset.gl (30 pops), testset.gs (SilicoDArT), pair subset
  via gl.keep.pop.
- Family mode: report (with the analysis lens on the Chao estimator and
  the bootstrap test). The helper utils.pa.Chao is reviewed here as the
  report's Chao engine (the utils.het.report precedent).
- Author/Custodian: Bernd Gruber.
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — 32 warning lines leak at verbose 0 (gl.keep.loc
  called inside utils.pa.Chao); three `cat(error()) + return(-1)`
  package guards; duplicate `@examples` and two `@family` tags; stale
  `save2temp` paragraph; an orphaned "[default FALSE]" line; the
  palette.discrete doc default does not match the signature; datatype2
  is computed and never used.
- **Spec: FAIL** — (1) The Chao1/Chao2 estimates are invalid three ways
  (F1, all verified): the helper recomputes maf across ALL populations
  in the dataset rather than the pair, so a pair's Chao values change
  when unrelated populations are present (full-x 0/0 vs pair-only 0/4
  for the same pair); pairs with zero private alleles fall through
  gl.keep.loc unchanged so their Chao is computed over the full locus
  set; and f1/f2 are taken positionally (table rows 1 and 2) rather
  than by category, so they are misassigned whenever singletons or
  doubletons are absent. (2) plot.file with plot.display = FALSE
  crashes ("object 'p3' not found") after the whole computation.
  (3) A mistyped method crashes with "object 'pall' not found".
  (4) one2rest + plot.display with the default palette crashes
  ("argument is of length zero") — the NULL default is handled in the
  pairwise branch only. (5) The test.asym documentation describes
  shuffling individuals with replacement using the smaller sample size;
  the code performs per-locus permutation without replacement at the
  actual sizes, and its arithmetic halves SilicoDArT frequencies so
  presence-fixed private alleles are missed. The core pairwise
  priv/fixed/AFD arithmetic is verified correct against an independent
  recomputation.

## Findings

### F1 — Chao estimates invalid (HIGH, confidence: certain) [escalation: numbers change]

Rule: spec axis (numerical validity); VRB1 (the gl.keep.loc leak).
Location: utils.pa.Chao.r whole file; gl.report.pa.r call sites.

Proposed change: rewrite the helper with direct matrix arithmetic on
the pair — count each private allele's copies within the pooled pair
sample, take f1/f2 as the counts of private alleles seen exactly
once/twice (category-matched, not positional), return 0 when there are
no private alleles, and use the bias-corrected form f1(f1-1)/2 when
f2 = 0. No gl.keep.loc/gl.recalc.metrics calls (removes the verbose-0
leak and the empty-set fallback). Chao1/Chao2 values change for
essentially every dataset with more than two populations.

### F2 — plot.file crashes without plot.display (MEDIUM, confidence: certain)

Proposed change: guard utils.plot.save on the existence of p3, with a
gated note that there is no plot to save (the established p3 idiom).

### F3 — Unvalidated method and x2 datatype (MEDIUM, confidence: certain)

Proposed change: fatal, informative stop for method outside
c("pairwise", "one2rest"); when x2 is provided, stop if its datatype
differs from x (datatype2 is currently computed and ignored).

### F4 — one2rest Sankey crashes on the default palette (MEDIUM, confidence: certain)

Proposed change: apply the same gl.select.colors(x, verbose = 0)
default used by the pairwise branch.

### F5 — cat(error) + return(-1) package guards (MEDIUM, confidence: certain)

Proposed change: the three package checks (tibble, networkD3, tidyr)
become stop(error(...)) so failures are fatal rather than returning -1
into downstream code.

### F6 — one2rest orientation rests on "zzest" sorting last (LOW, confidence: certain)

seppop sorts alphabetically; a population whose name sorts after
"zzest" (or named "zzest") silently swaps priv1/priv2 for its row.
Proposed change: index the seppop result explicitly by name so the
focal population is always pop1.

### F7 — test.asym docs do not match the test performed (MEDIUM, confidence: certain)

Proposed change: rewrite the details paragraph to describe the
implemented per-locus permutation test (without replacement, actual
sample sizes); skip test.asym for SilicoDArT with a gated warning
(its /2 arithmetic misses presence-fixed private alleles, so the
counts feeding the test are wrong for that datatype).

### F8 — Documentation (LOW, confidence: certain)

Duplicate @examples block removed; single @family; stale save2temp
paragraph removed; orphaned "[default FALSE]" removed and loc.names
given its default tag; palette.discrete default aligned to the
signature (NULL, resolved to gl.select.colors(x)); one2rest documented
as reporting no Chao, no test.asym, and no interactive map.

## Report notes (other functions / not fixed here)

- gl.keep.loc prints "Warning: no loci listed to keep!" even when
  called with verbose = 0 — an ungated warning for its own review.
- The gl.filter.pa example passes a factor element as pop1, which
  works only via factor-level integer coding — noted in that review.

## Coverage

`tests/testthat/test-gl.report.pa.R` — 12 assertions: pairwise table
vs independent recomputation (anchor), verbose-0 leak (baseline), Chao
pair-dependence (baseline), plot.file crash (baseline), method crash
(baseline), one2rest runs + Sankey default-palette crash (baseline),
SilicoDArT report runs with Chao NA. All 12 pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied. Suite: 13/13 pass; flips map to F1 (silent at
verbose 0; Chao identical between full-x and pair-only runs), F2
(plot.file without plot.display completes), F3 (informative method
stop), F4 (one2rest Sankey with default palette). End-to-end verbose 3
run clean; Chao values now scale with the pair's observed private
counts. Caller grep all-clear (gl.sample uses $fixed; gl.select.panel
uses $names_loci). PR #299.

```json
{"function": "gl.report.pa", "package": "dartR.base", "family_mode": "report",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec", "VRB1"], "loc": "R/utils.pa.Chao.r", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "PLT"], "loc": "R/gl.report.pa.r plot.file", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.report.pa.r validation", "status": "applied"},
  {"id": "F4", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.report.pa.r one2rest palette", "status": "applied"},
  {"id": "F5", "severity": "MEDIUM", "rules": ["DEP", "STY"], "loc": "R/gl.report.pa.r pkg guards", "status": "applied"},
  {"id": "F6", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.report.pa.r one2rest sentinel", "status": "applied"},
  {"id": "F7", "severity": "MEDIUM", "rules": ["DOC", "spec"], "loc": "R/gl.report.pa.r test.asym", "status": "applied"},
  {"id": "F8", "severity": "LOW", "rules": ["DOC"], "loc": "R/gl.report.pa.r header", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-gl.report.pa.R",
 "pr": 299}
```
