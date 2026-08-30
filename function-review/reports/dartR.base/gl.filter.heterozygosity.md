# Review: gl.filter.heterozygosity (dartR.base)

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.filter.heterozygosity
- Datasets: testset.gl (SNP), crafted 3-individual fixture with an
  all-NA individual, flag-less variant
- Family mode: modify (filters INDIVIDUALS by observed heterozygosity)
- Checks skipped: FBM path not exercised (not available: no FBM
  fixture); Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the monomorphs warning prints at verbose = 0 on
  every call whose monomorphs flag is FALSE (testset.gl included);
  visible return; the indented `@family matched filter` leaks into the
  rendered help title (0 concepts in the Rd, verified); t.lower's range
  error names t.upper; DOC1/DOC2/DOC7 drift.
- **Spec: FAIL** — two robustness gaps and one data-integrity gap: (1)
  an all-NA individual (NaN heterozygosity) crashes the subsetting
  ("Length of the provided vector does not match nInd(x)"); (2)
  removing individuals leaves the locus-metric flags stale (CallRate
  remains TRUE although the per-locus metrics no longer reflect the
  retained individuals — the gl.drop.ind class fixed in PRs #237/#238),
  and no recalc option is offered; (3) t.lower > t.upper falls through
  to a cryptic zero-individuals error; a flag-less object crashes on
  the unguarded `!flags$monomorphs`. The selection logic itself is
  exact: retained individuals match independent recomputation
  (inclusive boundaries), nLoc unchanged, ind.metrics stays in sync
  (dartR's subset operator handles it), one history entry.

## Findings

### F1 — Stale locus-metric flags after individual removal (MEDIUM, confidence: certain)

Rule: DAT4 (flags must reflect metric validity); precedent
gl.drop.ind/gl.keep.ind (PRs #237/#238). Location:
`R/gl.filter.heterozygosity.r` (post-filter section).

Verified: CallRate flag TRUE after an individual is removed. Every
individual-composition-dependent locus metric (CallRate, FreqHets,
FreqHomRef, FreqHomSnp, OneRatios, PICs, AvgPIC) is invalidated by the
removal, but the flags stay TRUE and downstream functions will trust
stale metrics. Proposed change: when individuals are removed, reset the
same flag set gl.drop.ind resets, at every verbosity. (Adding
recalc/mono.rm parameters as in gl.drop.ind would be an API extension —
offered as an option.)

### F2 — All-NA individuals crash the filter (MEDIUM, confidence: certain)

Rule: spec axis (NaN/NA-handling class). Location: heterozygosity
computation + subsetting.

An individual with all genotypes missing has c.hets = NaN; the NA in
the logical index crashes the genlight subset (verified). Proposed
change: individuals whose heterozygosity cannot be computed are removed
and itemised at verbose >= 3 (the campaign's NA policy applied to
individuals).

### F3 — Monomorphs warning ungated; unguarded flag access (LOW, confidence: certain)

Rule: VRB2; DEP guard class. Location: monomorphs check.

The warning prints at verbose = 0 (verified — on testset.gl every call
emits it), and `!x@other$loc.metrics.flags$monomorphs` crashes with
"invalid argument type" when the flags slot is absent (verified).
Proposed change: `isFALSE(...)` guard, gated at verbose >= 2.

### F4 — Threshold validation gaps (LOW, confidence: certain)

Rule: spec axis. Location: parameter checks.

t.lower's range error says "Parameter 't.upper' must lie between 0 and
1" (copy-paste, verified); t.lower > t.upper is not caught and falls
through to "Subsetting resulted in zero individuals" (verified).
Proposed change: correct the message; swap reversed thresholds with a
gated warning (the gl.filter.taglength precedent).

### F5 — Visible return (LOW, confidence: certain)

Rule: FS/VRB5 (campaign precedent). Proposed change:
`return(invisible(x.kept))`.

### F6 — Header conformance (LOW, docs-only, confidence: certain)

Rules: DOC1, DOC2, DOC7. The indented `@family matched filter` leaks
into the rendered help title (verified: 0 \concept entries); DOC2
wording; `@author` Custodian only; the deleted-individuals listing uses
paste0 with glued commas (trailing-comma wart, tidied in passing);
@seealso to gl.report.heterozygosity present via docs but not tagged.
Proposed change: rewrite the header to the ratified template.

## Report notes (other functions / not fixed here)

- ind.metrics synchronisation through `x[index, ]` works correctly
  (verified) — no DAT2 issue here.
- The matched report (gl.report.heterozygosity) review is on hold
  pending author consultation; nothing in this filter's findings
  depends on that outcome.

## Coverage

Characterization baseline:
`tests/testthat/test-gl.filter.heterozygosity.R` — 14 assertions:
retained individuals vs independent recomputation, ind.metrics sync,
nLoc unchanged, history +1, input untouched, stale flags (baseline),
all-NA crash (baseline), ungated warning (baseline), flag-less crash
(baseline), reversed-threshold error (baseline), t.upper message typo
(baseline), visible return (baseline). All 14 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — all six findings
**approved**. F1 taken as the minimal option (flag reset only, matching
the gl.drop.ind precedent; the recalc/mono.rm API extension was offered
and deferred).

## Outcome

Applied F1–F6. Verification: all 19 characterization assertions pass
post-fix; every diff from baseline maps to an approved finding (flags
reset when individuals removed, untouched on a no-op run [F1]; all-NA
individuals removed cleanly and itemised [F2]; monomorphs warning gated
with isFALSE guard — flag-less objects no longer crash [F3]; reversed
thresholds swapped silently at verbose 0 with counts matching
independent recomputation, t.lower message corrected [F4]; invisible
return [F5]). The @family concept now renders; sibling matched-filter
Rd cross-references update mechanically (F6). A legitimately empty
retention range still raises dartR's "zero individuals" subset error —
accurate, left as-is. Caller grep across dartR.base + 7 siblings: no
code callers. dartr2shiny: not present in the workspace. NEWS entry
added. Note: the matched report's review (gl.report.heterozygosity)
remains on hold pending author consultation; nothing here depends on
it.

```json
{
  "function": "gl.filter.heterozygosity",
  "package": "dartR.base",
  "family_mode": "modify",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["DAT4"], "loc": "R/gl.filter.heterozygosity.r post-filter", "status": "proposed"},
    {"id": "F2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.heterozygosity.r het computation", "status": "proposed"},
    {"id": "F3", "severity": "LOW", "rules": ["VRB2", "DEP"], "loc": "R/gl.filter.heterozygosity.r monomorphs check", "status": "proposed"},
    {"id": "F4", "severity": "LOW", "rules": ["spec"], "loc": "R/gl.filter.heterozygosity.r checks", "status": "proposed"},
    {"id": "F5", "severity": "LOW", "rules": ["VRB5"], "loc": "R/gl.filter.heterozygosity.r return", "status": "proposed"},
    {"id": "F6", "severity": "LOW", "rules": ["DOC1", "DOC2", "DOC7"], "loc": "R/gl.filter.heterozygosity.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "crafted all-NA fixture"],
  "baseline_test": "tests/testthat/test-gl.filter.heterozygosity.R",
  "pr": null
}
```
