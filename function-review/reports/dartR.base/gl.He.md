# Review: gl.He (dartR.base) — reviewed jointly with gl.Ho

## Provenance

As for gl.Ho.md (joint branch review-gl.Ho-gl.He, shared baseline
tests/testthat/test-gl.Ho.He.R).

## Verdicts

- **Standards: FAIL** — as gl.Ho, plus the `@return` documents the WRONG
  statistic ("A simple vector whit Ho" — it returns He).
- **Spec: FAIL** — no datatype guard (SilicoDArT silently returns
  meaningless values, verified); and the docs do not state that He here
  is the plain 2p(1-p) with p estimated from pooled mean dosage — i.e.
  neither sample-size-corrected (not uHe) nor per-population. The
  computation itself is exact: matches hand computation and equals
  2*alf1*alf2 from gl.alf (verified).

## Findings

### F1 — No datatype guard: SilicoDArT returns nonsense (MEDIUM, confidence: certain)

Location: `R/gl.He.r`. Mean presence/2 is not an allele frequency.
Proposed change: `utils.check.datatype(gl, accept = "SNP", verbose = 0)`.

### F2 — @return documents the wrong statistic; header below template (LOW, docs-only, confidence: certain)

`@return "A simple vector whit Ho for each loci"` — it returns expected
heterozygosity. Proposed change: minimal compliant header as for gl.Ho,
with @details clarifying: He = 2p(1-p), p from pooled mean allele
dosage; no sample-size correction (for unbiased uHe use
gl.report.heterozygosity); NaN for all-NA loci; SNP-only.

## Coverage / Approval / Outcome

Shared with gl.Ho (see gl.Ho.md): both findings approved and applied;
all 13 assertions pass post-fix.

```json
{
  "function": "gl.He",
  "package": "dartR.base",
  "family_mode": "analysis",
  "commit": "2739215",
  "skill_version": "2.0.0",
  "verdict_standards": "FAIL",
  "verdict_spec": "FAIL",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl.He.r", "status": "proposed"},
    {"id": "F2", "severity": "LOW", "rules": ["DOC1", "DOC5", "DOC7"], "loc": "R/gl.He.r header", "status": "proposed"}
  ],
  "datasets": ["testset.gl", "testset.gs"],
  "baseline_test": "tests/testthat/test-gl.Ho.He.R",
  "pr": 273
}
```
