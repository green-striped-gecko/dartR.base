# Review: gl2plink (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: platypus.gl (with and without assigned position/chromosome), testset.gl (NULL position/chromosome state), testset.gs (SilicoDArT)
- Baseline: tests/testthat/test-gl2plink.R (snapshot captured pre-review; 16 assertions, all passing against current code)

## Verdict

**Standards: Needs work** — the FS skeleton and `gl.check.wd` idiom are in place, but the datatype gate admits SilicoDArT into SNP-only math and two warnings bypass verbosity gating.
**Spec: Needs work** — genotype coding in the .ped is correct (hom-ref/het/hom-alt from `loc.all`, missing as `0 0`, verified against the matrix), but the .map fabricates coordinates when `@position`/`@chromosome` are unset, and the sex column reaches the file unconverted for scalar input.

## Findings

**F1 [HIGH, confidence: high] — unset position/chromosome fabricate .map coordinates (DOC5)**
`R/gl2plink.r:123-133` — when `@chromosome` is NULL the code substitutes the literal chromosome `"1"`, and when `@position` is NULL it substitutes `1:nLoc(x)` as base-pair positions. The `@details` claim "The chromosome information for unmapped SNPS is coded as 0" and "SNP position is taken from the accessor x$position", with no mention of either fallback. Under the genome-only convention `@position`/`@chromosome` stay NULL until explicitly assigned — this is the normal state (testset.gl ships this way), not an edge case.
Failure scenario: `gl2plink(testset.gl)` writes a .map asserting every locus sits on chromosome 1 at positions 1..755 (verified: first line `1 100049687-12-C/T 0 1`). Any downstream position-aware PLINK operation (LD pruning, --make-bed ordering, clumping) operates on fictitious coordinates; the only signal is a warning that also prints — and scrolls away — at `verbose = 0`. A related hazard, worth a `@details` sentence rather than a code change: objects whose `@position` holds within-tag offsets (platypus.gl ships with 5-68) get those offsets written as genomic base pairs.
Proposed change: code unmapped loci as PLINK's own convention — chromosome 0, position 0 — matching the existing NA-handling branch at lines 150-152, and document the behaviour; keep (and gate per VRB4) the warning.

**F2 [MEDIUM, confidence: high] — scalar sex.code bypasses conversion (DOC5)**
`R/gl2plink.r:188` — the F/M/U-to-2/1/0 recoding runs only `if (length(sex.code) > 1)`. The default scalar `"unknown"` is written verbatim into the .ped sex column, and a user-supplied scalar `"F"` (an all-female dataset, per the documented "starts with F" contract) is written as `F`. Verified both.
Failure scenario: PLINK treats any value outside 1/2 as missing sex, so an all-female cohort silently loses its sex information; the file also fails strict .ped validators expecting numeric codes.
Proposed change: drop the length guard — recode unconditionally (recycling handles scalars), mapping unknown/empty to "0".

**F3 [MEDIUM, confidence: high] — SilicoDArT admitted, then crashes mid-write (DAT7)**
`R/gl2plink.r:119` — `utils.check.datatype(x, verbose = verbose)` keeps the permissive default although the algorithm needs `loc.all` and 0/1/2 dosage. Verified: testset.gs fails with "number of items to replace is not a multiple of replacement length" after the .map has already been written.
Failure scenario: a SilicoDArT user gets an uninterpretable error plus a half-written file set.
Proposed change: `accept = "SNP"` (the recurring DAT7 defect class: gl2structure, gl2vcf, gl.propShared).

**F4 [MEDIUM, confidence: medium] — bed-mode allele list written but never used (DOC5)**
`R/gl2plink.r:238-247,249-268` — `bed.files = TRUE` writes `mylist.txt` (locus, first allele) to outpath, but the assembled PLINK command passes only `--file/--allow-no-sex/--allow-extra-chr/--out`. Without `--a2-allele mylist.txt` (its evident purpose), PLINK assigns A1/A2 by minor-allele frequency, so the .bed/.bim ref/alt orientation no longer matches the genlight's `loc.all`; `mylist.txt` remains as an unexplained artefact. Static finding — the PLINK branch was not executed (no binary; see Coverage).
Proposed change: pass `--a2-allele` pointing at the written list (and name the file from `outfile`), or remove the dead write.

**F5 [LOW, confidence: high] — ungated warnings leak at verbose = 0 (VRB3, VRB5)**
`R/gl2plink.r:125,130` — both dummy-fallback warnings are raw `cat(warn(...))`. Verified: two warning lines (one split by an embedded newline in the source string, lines 131-132) print at `verbose = 0`.
Proposed change: gate per VRB4 (results-affecting, so `verbose >= 1`) and repair the broken string literal.

**F6 [LOW, confidence: high] — description overstates the default products (DOC5)**
`R/gl2plink.r:6-8` — "This function produces the following PLINK files: bed, bim, fam, ped and map." The default (`bed.files = FALSE`) produces .ped and .map only; no separate .fam file is ever written by R (the six fam columns are embedded in the .ped); bed/bim/fam require the external binary.
Proposed change: state the two-file default and the bed-mode requirement explicitly.

**F7 [LOW, confidence: high] — roxygen order and wording drift (DOC1, DOC2, DOC7 (proposed rule))**
`R/gl2plink.r:1-91` — `@return` after `@export` and `@details` after `@param` (pre-ratification order; per DOC1 flag the file); `verbose` default clause not the ratified DOC2 wording; custodian-only author block (DOC7 — default fix: add `Author(s): Luis Mijangos.`). `@family linker` vs `linkers` elsewhere.

## Proposed changes

1. Write chromosome 0 / position 0 for unset slots, document the fallback and the tag-offset hazard (F1). **Consequence: .map contents change for objects with unset position/chromosome.**
2. Recode sex.code unconditionally (F2). **Consequence: the .ped sex column changes for scalar callers (literal text becomes numeric codes).**
3. Restrict datatype with `accept = "SNP"` (F3).
4. Wire `mylist.txt` into the PLINK call via `--a2-allele`, or delete the write (F4).
5. Gate the dummy-fallback warnings and fix the split string literal (F5).
6. Correct the `@description` file list (F6).
7. Reorder/reword roxygen per DOC1/DOC2/DOC7; run `devtools::document()` (F7, DOC4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plots), STY — run
- Spec: .ped/.map generated and inspected for platypus.gl with assigned coordinates (genotype coding verified cell-by-cell against `as.matrix`; missing as `0 0`; FID/IID from pop/indNames) and for testset.gl in the unset-slot state — run
- sex.code scalar vs vector behaviour: run
- SilicoDArT admission: run — crashes (F3)
- `bed.files = TRUE` / PLINK binary invocation (F4, cross-platform path handling): SKIPPED — no PLINK binary installed; reviewed statically
- `chr.format = "numeric"` with non-numeric chromosome names: SKIPPED — not exercised; the coercion path would produce NA chromosomes silently, worth covering in Phase C tests if change 1 is approved
- FBM path (DAT6): SKIPPED — no FBM fixture; converter densifies via `as.matrix` by design

## Approval (Phase B)

All findings at every severity approved by Arthur Georges, 2026-09-05, via
the formal approval boxes -- a blanket class approval covering this batch,
explicitly acknowledging that converted outputs change where they were
wrong, and ratifying the DAT7 fatal-gate policy (applied here as change 3).

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05; .map contents change for unset-slot objects |
| 2 | Approved | Arthur Georges | 2026-09-05; .ped sex column changes for scalar callers |
| 3 | Approved | Arthur Georges | 2026-09-05; DAT7 policy: fatal accept = "SNP" gate |
| 4 | Approved | Arthur Georges | 2026-09-05; --a2-allele wired, list named from outfile |
| 5 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 6 | Approved | Arthur Georges | 2026-09-05, blanket class approval |
| 7 | Approved | Arthur Georges | 2026-09-05, blanket class approval |

## Outcome (Phase C)

All seven changes applied on branch `review-gl2plink` (base
`upstream/dev` = ddaed27), 2026-09-05. Change 4 is implemented as the
report's first option: the allele list is written as
`<outfile>_a2_alleles.txt` and passed with explicit column numbers
(`--a2-allele <file> 2 1`); still static (no PLINK binary), as in
Phase A.

Verification: baseline characterization test updated -- every diff maps
to an approved finding (F1 unmapped 0/0 coordinates, F2 sex recoding for
scalars, F3 SilicoDArT gate, F5 gated warnings) -- 17 assertions pass.
End-to-end at verbose = 3 on platypus.gl (assigned coordinates) and
testset.gl (unset slots): clean, warnings print as gated. Sibling grep
across the 8 clones: two package-code callers, both verified unaffected --
`gl2vcf` assigns position and chromosome (its own "0" fallbacks) before
calling, so the F1 fallback never engages there, and its SilicoDArT
pass-through previously crashed inside the recoding and now stops at the
clear datatype gate; `gl.report.ld.map` requires mapped data by contract,
and the fabricated-coordinate path it could have inherited is precisely
the F1 defect. Neither breaks.
PR: https://github.com/green-striped-gecko/dartR.base/pull/350

```json
{
  "function": "gl2plink",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "applied", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 7}
  ],
  "coverage_skipped": ["bed.files=TRUE PLINK invocation: no binary", "chr.format=numeric coercion path: not exercised", "DAT6: no FBM fixture"],
  "status": "phase-c-complete",
  "pr": 350
}
```
