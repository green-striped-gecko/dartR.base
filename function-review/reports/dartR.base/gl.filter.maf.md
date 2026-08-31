# Review: gl.filter.maf (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Package commit: ddaed27 (dev, synced with upstream/dev);
  Branch: review-gl.filter.maf; Datasets: testset.gl (4-population
  subset)
- Family mode: modify (locus filter; matched filter to gl.report.maf)
- Custodian: Luis Mijangos
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the threshold-coercion warning prints at
  verbose 0 (verified); the return is visible against the ratified
  invisible-filter convention; dead code (`hold <- x` and
  `plot.colors.pop` captured and never used; a nonsense comment
  "Remove loci with NA count <= 1-threshold"); a full
  gl.filter.monomorphs preamble runs solely to warn (the #279 class);
  header drift ("histograms of base composition" copy-pasted from the
  bases function).
- **Spec: FAIL** — (1) the by.pop path collapses when NO loci qualify
  for removal: `x[, -integer(0)]` subsets to zero loci — verified
  crash ("Subsetting resulted in zero loci") when no population meets
  ind.limit; the same empty-negative-index reaches the subset whenever
  the qualifying set is empty. (2) plot.file with plot.display = FALSE
  crashes ("object 'p3' not found" — verified). (3) by.pop = TRUE with
  0 or 1 qualifying populations would crash at display on `p1`, which
  is built only in the global block (code-certain; currently masked by
  crash 1). (4) The per-population plot list is named positionally
  with the qualifying-pop names (`names(plots_pops_merge) <-
  popn.hold` over the full-length list), so plots are mislabeled and
  wrongly selected whenever a non-qualifying population precedes a
  qualifying one. (5) Boundary and NA inconsistencies: by.pop counts
  `maf <= threshold` while the global path (and the docs) removes
  `maf < threshold`; NA-maf loci are removed by the global path but
  retained by by.pop. The global arithmetic itself is verified correct
  (independent recomputation matches; MAC conversion works; history
  single-entry; input untouched).

## Findings

### F1 — by.pop empty-set collapse (HIGH, confidence: certain)

Rule: spec axis (the empty-negative-index class). Location: by.pop
subset.

Proposed change: subset only when `length(loci.list) > 0`, otherwise
return the object unchanged with a gated (verbose >= 1) message; when
no population meets ind.limit, a gated message states that MAF could
not be assessed by population and the object is returned unchanged.

### F2 — Plot crashes (HIGH, confidence: certain)

Rule: PLT (p3-not-found class). Location: printing and save blocks.

plot.file without display crashes on the unbuilt p3 (verified); the
by.pop <= 1-qualifying-pop display path references p1 from the global
block. Proposed change: restructure printing so each path prints the
plots it built; guard the save on the plot existing with a gated
notice.

### F3 — Per-population plots mislabeled (MEDIUM, confidence: certain)

Rule: spec. Location: by.pop plot naming.

Proposed change: name the merged plot list by the population names it
was built from and select the qualifying subset by name.

### F4 — Warning gate and MAC boundary (MEDIUM, confidence: certain)

Rule: VRB/DOC. The threshold-range warning gated at verbose >= 1
(printed at verbose 0, verified); the MAC interpretation applies at
`threshold >= 1` while the docs say "> 1" — code aligned to the
documented `> 1` (threshold == 1 then falls to the range coercion).

### F5 — Boundary alignment, invisible return, and tidy (LOW-MEDIUM, confidence: certain)

Rule: spec/STY/DOC. by.pop aligned to the documented strict
`maf < threshold` (was `<=`; no observable difference on the test
data); NA policy documented (global removes NA-maf loci — all-NA loci
— itemised at verbose >= 3; by.pop treats NA as non-qualifying);
invisible return (ratified filter convention); dead code and the
monomorphs preamble removed; docs corrected ("base composition",
@return invisible note, plot.dir wording).

### F6 — Datatype restriction (MEDIUM, decision, confidence: certain)

Rule: DAT dispatch. MAF is undefined for presence/absence data, yet
SilicoDArT objects are accepted and run through utils.recalc.maf
producing meaningless values. Proposed change: `accept = "SNP"`.
Behaviour change: SilicoDArT input becomes a fatal error. (The
matched report gl.report.maf shares this looseness — noted there for
a follow-up if approved here.)

## Report notes (other functions / not fixed here)

- gl.report.maf also accepts SilicoDArT (same looseness as F6) — if
  F6 is approved, the report side should follow in a small amendment.

## Coverage

`tests/testthat/test-gl.filter.maf.R` — 13 assertions: global-path
counts vs independent recomputation (17 kept at 0.05), NA-maf removal
count, MAC conversion (threshold 5 → 15 loci), by.pop removal count
(221) + history + input untouched, empty-set crash (baseline), p3
crash (baseline), verbose-0 warning leak (baseline), visible return
(baseline). All 13 pass on the pre-fix code.

## Approval

All six findings approved via the approval boxes (2026-08-31); F6
approved as a behaviour change (SilicoDArT rejected), with the matched
report to follow by amendment to PR #294.

## Outcome

All six findings applied. Characterization suite: 16/16 pass; flips
map to F1 (unchanged object instead of the zero-loci crash), F2 (no
p3 crash; silent gated notice), F4 (warning gated), F5 (invisible
return), F6 (SilicoDArT fatal). Global and by.pop removal counts
unchanged (17 kept at 0.05 globally; 221 removed under by.pop
pop.limit=1 — the <= to < boundary alignment produced no observable
difference on the test data). End-to-end at verbose 3/2 clean on the
global path, the multi-pop by.pop path, and the previously-crashing
few-qualifying-pop display path. NEWS entry added.

```json
{"function": "gl.filter.maf", "package": "dartR.base", "family_mode": "modify",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec"], "loc": "R/gl.filter.maf.r by.pop subset", "status": "applied"},
  {"id": "F2", "severity": "HIGH", "rules": ["PLT"], "loc": "R/gl.filter.maf.r plots", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.maf.r plot naming", "status": "applied"},
  {"id": "F4", "severity": "MEDIUM", "rules": ["VRB", "DOC"], "loc": "R/gl.filter.maf.r threshold", "status": "applied"},
  {"id": "F5", "severity": "LOW", "rules": ["spec", "STY", "DOC"], "loc": "R/gl.filter.maf.r", "status": "applied"},
  {"id": "F6", "severity": "MEDIUM", "rules": ["DAT"], "loc": "R/gl.filter.maf.r datatype", "status": "applied"}],
 "datasets": ["testset.gl"],
 "baseline_test": "tests/testthat/test-gl.filter.maf.R", "pr": null}
```
