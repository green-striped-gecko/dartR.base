# Review: gl.filter.pa (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-gl.filter.pa.
- Datasets: testset.gl, testset.gs.
- Family mode: modify (matched filter for gl.report.pa).
- Author: Bernd Gruber & Ella Kelly; Custodian: Luis Mijangos.
- Checks skipped: Google Group not searched (not available: no browser
  session).

## Verdicts

- **Standards: FAIL** — the filtered object returns visibly (filters in
  this family return invisibly); no results summary at verbose >= 2
  (loci before/after); `T` for TRUE; the loc.metrics subset lacks
  drop = FALSE; the description is informal ("not sure why yet") with
  an unmatched parenthesis; "[FALSE]" instead of "[default FALSE]".
- **Spec: FAIL** — (1) SilicoDArT frequencies are halved as if diploid,
  so presence-fixed private alleles are missed: on testset.gs the
  filter keeps 45 loci where the true private/fixed union is 66
  (verified). The datatype is checked and then ignored. (2) A bogus
  population name crashes with "'data' must be of a vector type, was
  'NULL'" from as.matrix(NULL). The SNP path itself is verified
  correct: kept loci match an independent recomputation of the
  private/fixed union exactly, loc.metrics stay aligned, invers
  returns the exact complement, and history is appended as a call.

## Findings

### G1 — SilicoDArT arithmetic wrong (HIGH, confidence: certain) [escalation: output changes]

Rule: spec axis; DAT (datatype dispatch). Location: the /2 divisions.

Proposed change: divide by 2 only for SNP data (the arithmetic
gl.report.pa already uses). For SilicoDArT users the kept locus set
changes (45 -> 66 on testset.gs) — previously missed presence-fixed
private loci are now retained.

### G2 — Bogus population name crashes cryptically (MEDIUM, confidence: certain)

Proposed change: validate pop1/pop2 against popNames(x) with an
informative fatal stop naming the offender; coerce factor input to
character (the documented example passes a factor element, which works
only via level coding).

### G3 — Return and reporting conventions (LOW, confidence: certain)

Proposed change: return invisibly; report loci before/after at
verbose >= 2; T -> TRUE; drop = FALSE on the loc.metrics subset.

### G4 — Documentation (LOW, confidence: certain)

Description tidied (informal aside and unmatched parenthesis removed;
the gl.nhybrids motivation kept); "[default FALSE]"; @return mentions
the invers complement; the example passes a character population name.

## Report notes (other functions / not fixed here)

- None.

## Coverage

`tests/testthat/test-gl.filter.pa.R` — 8 assertions: SNP union vs
independent recomputation + invers complement (anchor), silico
undercount (baseline), bogus-pop crash (baseline), visible return
(baseline), history-is-call. All 8 pass pre-fix.

## Approval

All four findings approved via the approval boxes (2026-09-01).

## Outcome

All four findings applied. Suite: 8/8 pass; flips map to G1 (silico
kept set now equals the true union, 66 on testset.gs), G2 (informative
stop naming the offending population), G3 (invisible return).
End-to-end verbose 3 run clean with the new summary line. Caller grep
all-clear (no callers of gl.filter.pa in the family). PR recorded on
merge branch.

```json
{"function": "gl.filter.pa", "package": "dartR.base", "family_mode": "modify",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "G1", "severity": "HIGH", "rules": ["spec", "DAT"], "loc": "R/gl.filter.pa.r frequencies", "status": "applied"},
  {"id": "G2", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl.filter.pa.r validation", "status": "applied"},
  {"id": "G3", "severity": "LOW", "rules": ["VRB", "STY"], "loc": "R/gl.filter.pa.r return", "status": "applied"},
  {"id": "G4", "severity": "LOW", "rules": ["DOC"], "loc": "R/gl.filter.pa.r header", "status": "applied"}],
 "datasets": ["testset.gl", "testset.gs"],
 "baseline_test": "tests/testthat/test-gl.filter.pa.R",
 "pr": null}
```
