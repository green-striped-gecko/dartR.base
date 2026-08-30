# Review: gl.Ho (dartR.base) — reviewed jointly with gl.He

## Provenance

- Model: claude-sonnet-5 (Claude Code)
- Skill: dartr-function-review v2.0.0
- Package commit: 2739215 (dev, synced with upstream/dev)
- Branch: review-gl.Ho-gl.He (joint branch — the two 15-line siblings
  share every finding; joint report for gl.He in gl.He.md)
- Datasets: testset.gl, testset.gs, crafted all-NA-locus fixture
- Family mode: analysis (pure per-locus accessor)
- Checks skipped: FBM exercised implicitly via as.matrix only; Google
  Group not searched (not available: no browser session).

## Verdicts

- **Standards: FAIL** — header far below template: no @name/@title/
  @family/@description tags, "A simple vector whit Ho" typo, parameter
  named `gl` rather than the house `x` (retained — renaming would break
  named-argument callers). The pure-function style itself (no verbosity
  scaffold, silent, visible vector return) is judged appropriate for an
  accessor whose vector IS the product, and is documented as deliberate
  rather than "fixed".
- **Spec: FAIL** — no datatype guard: SilicoDArT presence/absence data
  silently returns meaningless "heterozygosities" (verified: values
  spanning 0-1 from 0/1 presence codes). The computation itself is
  exact: matches hand computation, carries locus names, and the mean of
  per-locus gl.Ho over a population equals gl.report.heterozygosity's
  population Ho (verified). All-NA loci yield NaN (documented, not
  changed).

## Findings

### F1 — No datatype guard: SilicoDArT returns nonsense (MEDIUM, confidence: certain)

Location: `R/gl.Ho.r`. Heterozygosity is undefined for presence/absence
data; the function happily averages `m == 1` over presence codes.
Proposed change: `utils.check.datatype(gl, accept = "SNP", verbose = 0)`
— SilicoDArT input becomes a fatal error.

### F2 — Header far below template (LOW, docs-only, confidence: certain)

No @name/@title/@family; description-only line with the "whit" typo.
Proposed change: minimal compliant header — @name/@title,
`@family utilities` (the gl.alf precedent), corrected @return, @details
noting: whole-dataset pooling (per-population values via seppop or
gl.report.heterozygosity), NaN for all-NA loci, SNP-only, and that the
function is a silent pure accessor by design; @seealso gl.He, gl.alf,
gl.report.heterozygosity.

## Report notes

- Zero callers across the eight packages and tutorials — these are
  standalone conveniences; nothing depends on the guard change.
- The `gl` parameter name is retained deliberately (API stability).

## Coverage

Characterization baseline: `tests/testthat/test-gl.Ho.He.R` (shared
with gl.He) — 13 assertions: hand-computation identity, names/length,
aggregation consistency with the population report, pure-function
behaviour, NaN for all-NA loci, silico numeric passthrough (baseline).
All 13 pass on the pre-fix code.

## Approval

Decisions recorded 2026-08-30 via approval boxes — both findings
**approved** for both functions (F1 datatype guard; F2 documentation).

## Outcome

Applied F1-F2 to both functions. Verification: all 13 shared
characterization assertions pass post-fix; the only behavioural diff is
the approved SilicoDArT fatal error. Values re-verified exact against
hand computation, gl.alf, and the population report. Joining
@family utilities updates the family's Rd cross-references mechanically.
Caller grep: zero callers anywhere in the eight packages and tutorials.
dartr2shiny: not present in the workspace. NEWS entry added.

```json
{
  "function": "gl.Ho",
  "package": "dartR.base",
  "family_mode": "analysis",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.Ho.r", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["DOC1", "DOC7"], "loc": "R/gl.Ho.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs", "crafted all-NA fixture"],
  "baseline_test": "tests/testthat/test-gl.Ho.He.R",
  "pr": 273
}
```
