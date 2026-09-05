# Review: gi2gl (dartR.base)

- Family mode: io (reviewed as a pair with `gl2gi`; round-trip fidelity is
  part of the spec axis)
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: testset.gl (via gl2gi), synthetic df2genind fixtures
  (single-locus, biallelic, triallelic)
- Baseline: tests/testthat/test-gi2gl.R (snapshot captured pre-review,
  including round-trip characterization)

## Verdict

**Standards: Needs work** — the body runs `gl.compliance.check` and appends
history correctly, but the file's roxygen header carries `@name gl2gi`, which
merges this function's documentation into the wrong manual page.
**Spec: Needs work** — for a clean biallelic genind the conversion is sound
(verified: structure, names, populations, NA pattern all survive), but a
genind containing any locus with more than two alleles is silently corrupted,
and a single-locus genind crashes.

### Round-trip fidelity (gl -> gi -> gl), verified empirically

- On raw `testset.gl` the round trip **crashes** (root cause: `gl2gi` F1,
  the loc.metrics desync on dropped all-NA loci).
- On `gl.filter.allna(testset.gl)` (752 loci) the round trip returns a
  compliant ploidy-2 dartR object. Survives exactly: dimensions, indNames,
  locNames, pop, NA pattern, `ind.metrics`, `loc.metrics` (752 rows),
  `latlon`; history is appended.
- Lost: dosage polarity — `gi2gl` counts the genind tab's first allele
  column, so 596 of 752 loci come back exactly complemented (0 <-> 2, hets
  unchanged); and allele spellings — the returned `loc.all` is a uniform
  placeholder, not the original `C/T`, `G/A`, ... pairs. The polarity loss
  is declared in `@details`; the allele-spelling loss is not.

## Findings

**F1 [HIGH, confidence: high] — loci with more than two alleles silently corrupt the output (DAT5)**
`R/gi2gl.R:42-55` — the column-offset walk over `gi@tab` advances by 1 or 2
columns per locus (`locna[i-1] == 1` else `+2`), so a locus with three or
more alleles shifts every subsequent locus onto the wrong tab column.
Failure scenario: for a two-locus genind whose first locus has 3 alleles,
the returned "locus 2" genotypes are locus 1's third-allele counts
(verified: expected first-allele counts 2,1,0,2; returned 0,0,1,2). No
error, no warning — plausible wrong data under correct locus names, for any
genind of microsatellite or multiallelic origin, which is exactly the
external input this converter exists to accept.
Proposed change: fail fast when `max(gi@loc.n.all) > 2` with an error naming
the offending loci (or drop them with a warning gated at `verbose >= 1`),
and step the offset by `loc.n.all[i-1]` rather than a hard-coded 2.

**F2 [MEDIUM, confidence: high] — single-locus genind crashes (FS5)**
`R/gi2gl.R:44` — `for (i in 2:length(locna))` iterates `c(2, 1)` when
`length(locna) == 1`, so `locna[0]` (length zero) reaches the `if`.
Failure scenario: `gi2gl` on a one-locus genind fails with "argument is of
length zero". Verified empirically.
Proposed change: `for (i in seq_len(length(locna))[-1])` or guard
`if (length(locna) > 1)`.

**F3 [MEDIUM, confidence: high] — @name collides with gl2gi (DOC1, DOC4)**
`R/gi2gl.R:1` — the header declares `@name gl2gi`. roxygen merges both
functions into one `man/gl2gi.Rd` ("Please edit documentation in
R/gi2gl.R, R/gl2gi.r"): the shared page carries two descriptions and two
`@return` values, `?gl2gi` is titled "Converts a genind object into a
genlight object" (the wrong function), and no standalone `gi2gl.Rd` exists.
The header also has no `@examples` and no `@description` tag.
Failure scenario: a user reading `?gl2gi` is told it converts genind to
genlight; `help(gi2gl)` resolves only through the alias to the merged page.
Proposed change: `@name gi2gl`, add `@description` and a runnable
`@examples` (e.g. `gi2gl(gl2gi(gl.filter.allna(testset.gl)))`), reorder tags,
re-document (DOC4).

**F4 [LOW, confidence: high] — round trip discards recoverable allele information (DOC5, proposed rule)**
`R/gi2gl.R:52-62` — the new genlight is built without `loc.all`, although
`colnames(gi@tab)` carry the true allele labels per locus; the placeholder
later injected by `gl.compliance.check` presents uniform fabricated alleles.
The 0 <-> 2 polarity flip (596/752 loci on the reference data) is likewise
derivable from the tab column names when the genind came from `gl2gi`.
Failure scenario: `gi2gl(gl2gi(x))` loses the real allele spellings even
though the information is present in the input; any downstream export that
prints alleles reports the placeholder.
Proposed change: reconstruct `loc.all` from the tab column names for
biallelic loci; document the residual polarity ambiguity as now.

**F5 [LOW, confidence: medium] — ploidy hard-coded to 2 (DAT1)**
`R/gi2gl.R:52-62` — `ploidy = 2` regardless of `ploidy(gi)`.
Failure scenario: a haploid genind converts to a nominally diploid genlight
whose 0/1 values then read as homozygote/heterozygote dosages.
Proposed change: reject non-diploid genind objects with a clear message.

**F6 [LOW, confidence: high] — author block lacks Custodian line (DOC7) (proposed rule)**
`R/gi2gl.R:11` — `@author` names Bernd Gruber with no `Custodian:` label.
Proposed change: add `Custodian: Bernd Gruber` per the DOC7 default.

What works well: the function runs `gl.compliance.check` on its product
(DAT5 done right for the construction side) and appends history to the
returned object (FS8).

## Proposed changes

1. Guard against loci with more than two alleles and step the tab-column
   offset by the actual allele count (F1).
   **Consequence: genind objects with multiallelic loci now error (or are
   dropped with a warning) instead of returning corrupted genotypes.**
2. Fix the single-locus loop bound (F2).
3. Rename the roxygen topic to `gi2gl`, add description and examples,
   re-document (F3).
4. Reconstruct `loc.all` from the genind tab column names for biallelic
   loci (F4).
5. Reject non-diploid genind input (F5).
6. Add the Custodian line (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (no plot bundle;
  no file output; adegenet is Depends so DEP1 not applicable; FS8 verified
  present on the returned object; input guard `is(gi, "genind")` verified).
- Spec: round trip on raw and all-NA-filtered testset.gl — run (see
  fidelity summary above); biallelic 2-locus sanity fixture — run;
  triallelic corruption fixture — run; single-locus fixture — run;
  non-genind rejection — run; `verbose = 0` silence — verified (zero
  captured lines, `gl.compliance.check` included).
- FBM path (DAT6): SKIPPED — not applicable; input is a genind.
- Haploid genind probe: SKIPPED — behaviour asserted from code reading only
  (F5 confidence medium accordingly).

## Approval (Phase B)

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 2 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 3 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 4 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 5 | Approved | Arthur Georges | 2026-09-05, via approval boxes |
| 6 | Approved | Arthur Georges | 2026-09-05, via approval boxes |

All findings at every severity approved 2026-09-05 (blanket class
approval, explicitly acknowledging that converted outputs change where
they were wrong; the DAT7 `accept = "SNP"` gate is not applicable —
input is a genind, guarded by the existing `is(gi, "genind")` check).

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gi2gl` (base `upstream/dev`,
ddaed27; independent of `review-gl2gi`, no stacking). PR: (recorded
below after opening).

- F1: fatal error naming the offending loci when any locus has more
  than two alleles; the tab-column offset is now `cumsum` over the
  actual per-locus allele counts rather than hard-coded steps of 1 or 2.
- F2: the cumsum offset also removes the `2:length(locna)` loop, fixing
  the single-locus crash.
- F3: roxygen topic renamed to `@name gi2gl` with its own
  `@description`, `@examples` and the house tag order; `man/gi2gl.Rd`
  now exists as a standalone page and `man/gl2gi.Rd` no longer carries
  the merged content (both Rd files regenerated on this branch).
- F4: `loc.all` reconstructed from the genind allele names
  (second/first of the tab pair, so dosage and pair are mutually
  consistent; loci monomorphic in the data come back as `a/a`),
  replacing the uniform placeholder.
- F5: non-diploid genind objects rejected with a clear fatal error.
- F6: `Author(s)`/`Custodian` block added.

Verification: all baseline assertions pass (617, most from the
per-locus round-trip loop). Updated expectations, each mapped to an
approved finding: multiallelic input now errors naming the locus (F1);
single-locus genind now converts (F2); `loc.all` is now set-equal per
biallelic locus to the original pair rather than a placeholder (F4).
The raw-testset round-trip test was rewritten to hold in both merge
states of the independent `review-gl2gi` branch (its F1 fix removes the
crash). End-to-end `gi2gl(gl2gi(gl.filter.allna(testset.gl)), verbose
= 3)` completes with in-register metadata. Sibling caller grep across
the 8 clones: no callers of `gi2gl` anywhere — all-clear. NEWS entry
added.

### Round-trip re-verification with both fixes (gl2gi PR #338 + this branch)

Run 2026-09-05 with `review-gl2gi`'s `gl2gi.r` and this branch's
`gi2gl.R` loaded together:

- `gl.filter.allna(testset.gl)` (752 loci): dimensions, indNames,
  locNames, pop and NA pattern identical; `loc.metrics` 752 rows,
  `ind.metrics` 274 rows (in register); 596 of 752 loci come back as
  exact 0 <-> 2 complements and every locus is identical-or-complement
  (the documented reference-allele ambiguity, unchanged); `loc.all` is
  set-equal to the original pair on all 611 loci biallelic in the data
  (previously a uniform placeholder).
- Raw `testset.gl` (755 loci incl. 3 all-NA): the round trip now
  SUCCEEDS — 752 loci with 752-row `loc.metrics` (previously crashed in
  `gl.recalc.metrics`).

```json
{
  "function": "gi2gl",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "DAT1", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC7", "status": "applied", "change": 6}
  ],
  "coverage_skipped": ["DAT6: not applicable (genind input)", "haploid genind probe: asserted from code reading"],
  "status": "phase-c-complete",
  "pr": null
}
```
