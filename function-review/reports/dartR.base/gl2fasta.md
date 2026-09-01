# Review: gl2fasta (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-gl2fasta.
- Nominated with a user bug report (Rachel, dartR Google Group):
  gl2fasta on a merged/co-analysed two-lab DArT dataset outputs data
  for a single locus with every method; other datasets fine; the
  TrimmedSequence/loc.all/position slots verified present by the
  user.
- Datasets: testset.gl subsets; constructed pathologies (gutting
  SnpPosition, factor SnpPosition, all-overshoot, flag-less object).
- Family mode: io (export).
- Custodian: Luis Mijangos.

## Verdicts

- **Standards: FAIL** — cat(error()) + return(-1) package guard;
  NULL-unsafe monomorphs flag check (verified crash on flag-less
  objects); `matrix` shadows base::matrix; sequences written with
  cat(result, " \\n") (trailing spaces); @return documents "a new gl
  object with all loci rendered homozygous" but the function returns
  NULL (verified); two callee lines leak at verbose 0 (from
  gl.filter.overshoot's ungated else-branch — already fixed on the
  open PR #262).
- **Spec: FAIL — the user's symptom REPRODUCED.** (1) gl2fasta
  silently runs gl.filter.overshoot(verbose = 0) as a pre-filter.
  On a dataset where SnpPosition is inconsistent with TrimmedSequence
  lengths — precisely what a merged/co-analysed dataset can produce —
  the pre-filter guts the data and gl2fasta writes single-locus FASTA
  with NO message at any verbosity (the removal count prints only at
  the callee's verbose >= 3, which is hardcoded to 0): verified,
  method 3 wrote 1-character sequences and method 1 a single tag,
  console silent. This is the reported bug's mechanism. (2) A factor
  SnpPosition silently corrupts method-1 output (substr consumes
  factor LEVEL CODES: 118/59/155-char sequences where 523 were
  expected). (3) The documented method = 0 help path crashes on a
  cat():cat() colon typo ("argument of length 0"). (4) sink() is
  opened without on.exit: an all-overshoot dataset crashes mid-write
  and leaves the sink open, hijacking the session's console output
  (verified: "Error in stdout(): invalid connection"). All four
  methods verified CORRECT on clean data (sequence lengths exact for
  concatenated tags and SNP-only outputs).

## Findings

### F1 — Silent gutting pre-filter (HIGH, confidence: certain) [the user-report fix]

Proposed change: compare nLoc before/after the overshoot pre-filter;
report the removed count at verbose >= 2; WARN at verbose >= 1
whenever loci are removed; and STOP with an informative message when
fewer than 2 loci survive — naming the SnpPosition-vs-TrimmedSequence
inconsistency and noting merged/co-analysed datasets as a known
trigger. The user's mystery becomes a one-line diagnosis.

### F2 — SnpPosition type not validated (MEDIUM, confidence: certain)

Proposed change: coerce via as.numeric(as.character(...)) up front;
fatal, informative stop if that produces NAs. Silent corruption
becomes impossible.

### F3 — sink without on.exit (MEDIUM, confidence: certain)

Proposed change: on.exit(sink(), add = TRUE) guarding both sink
blocks, so any mid-write error restores the console.

### F4 — method = 0 help path crashes (MEDIUM, confidence: certain)

Proposed change: the cat():cat() colon replaced with separate cat
calls; the listing prints, then the existing informative stop fires.

### F5 — Guards, docs, tidy (LOW, confidence: certain)

Package guard becomes a fatal stop; monomorphs check NULL-safe
(!isTRUE); @return corrected (NULL, invisibly); FASTA lines written
without trailing spaces (sep = ""); `matrix` renamed.

## Report notes (other functions / not fixed here)

- The two verbose-0 leak lines originate in gl.filter.overshoot's
  ungated else-branch cat — already gated on the open PR #262; the
  leak disappears when it merges.
- The zero-loci escape through `[` with negative indices (which let
  the all-overshoot case reach the sink) is closed by PR #328.

## Coverage

`tests/testthat/test-gl2fasta.R` — 10 assertions: all four methods'
sequence lengths exact on clean data (anchor), silent gutting
(baseline), factor-corruption (baseline), method-0 crash (baseline),
NULL return (anchor), flag-less crash (baseline). All pass pre-fix.
The sink-hijack case is asserted post-fix only (exercising it
pre-fix corrupts the test session's output streams).

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied, with one refinement discovered during apply:
the F2 coercion (as.numeric(as.character())) RECOVERS factor-coded
numeric positions correctly - verified byte-identical FASTA against
the clean reference - so factor input is silently repaired rather
than fatal; only genuinely non-numeric values stop. Suite: 13/13
pass; flips map to F1 (the reproduced user scenario now warns "7 of
8 loci removed" and stops with the merged-dataset diagnosis), F2, F4,
F5; the all-overshoot case stops cleanly with sinks intact (F1+F3).
All four methods' output verified byte-correct on clean data. PR
recorded below.

```json
{"function": "gl2fasta", "package": "dartR.base", "family_mode": "io",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [
  {"id": "F1", "severity": "HIGH", "rules": ["spec", "VRB"], "loc": "R/gl2fasta.r overshoot pre-filter", "status": "applied"},
  {"id": "F2", "severity": "MEDIUM", "rules": ["spec", "DAT"], "loc": "R/gl2fasta.r SnpPosition type", "status": "applied"},
  {"id": "F3", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl2fasta.r sink", "status": "applied"},
  {"id": "F4", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/gl2fasta.r method 0", "status": "applied"},
  {"id": "F5", "severity": "LOW", "rules": ["STY", "DOC"], "loc": "R/gl2fasta.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-gl2fasta.R",
 "pr": null}
```
