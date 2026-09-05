# Review: gl2vcf (dartR.base)

- Family mode: io
- Date: 2026-09-05
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v1.0.0
- Package commit: ed99203
- Datasets: platypus.gl (dartR.data), testset.gs
- Baseline: tests/testthat/test-gl2vcf.R (snapshot captured pre-review; 41
  assertions passing, including an end-to-end run against a real PLINK
  1.9 binary — gated on env var PLINK19_DIR)

## Verdicts

**Standards: Needs work** — the datatype gate admits SilicoDArT (a named
DAT7 offender), the delegated `gl2plink` call and the raw PLINK output
defeat `verbose = 0`, and the external-binary dependency is unguarded.
**Spec: Needs work** — with a clean object and explicit `snp.pos`/`snp.chr`
the emitted file is a valid VCFv4.2 whose CHROM/POS/ALT and GT coding all
verified correct; but the coordinate-precedence logic silently ignores the
documented arguments whenever `@position` is already populated — on the
shipped `platypus.gl` (which still carries the stale tag-offset copy) the
function's own roxygen example emits tag offsets as POS — and REF degrades
to `N` for loci fixed for the alternate allele.

What works well: genotype coding verified end-to-end against `as.matrix` —
dosage 0 -> `0/0`, 1 -> `0/1`, 2 -> `1/1`, NA -> `./.` — with allele
orientation matching `loc.all` whenever the reference allele is observed;
missing `snp.pos`/`snp.chr` fields stop with a clear message; the post-#330
`is.null` test correctly zeroes CHROM/POS instead of sniffing magnitudes.

## Findings

**F1 [HIGH, confidence: high] — a populated @position slot silently overrides the explicit snp.pos/snp.chr arguments (DOC5 proposed; audit C1/PR #330 follow-on)**
`R/gl2vcf.r:105-127,130-156` — the fill logic runs only when the slot is
NULL (or wrong length); when `@position`/`@chromosome` are already set, the
user's `snp.pos`/`snp.chr` are ignored with no message at any verbosity.
Failure scenario: verified — `platypus.gl` as currently shipped in
dartR.data still carries the pre-#330 stale copy (`@position` = tag
offsets, `@chromosome` NULL). Running the function's own `@examples` call
(`gl2vcf(platypus.gl, snp.pos = 'ChromPos_...', snp.chr = 'Chrom_...')`)
produced a VCF with genuine genome CHROM names paired with POS = 36, 13,
54, 65 — the within-tag offsets — instead of 2438118, 32077451, ... . The
file is well formed and silently wrong; any position-aware downstream use
(merging, annotation, LD) is corrupted. Reproduced equally with a manually
pre-set slot plus conflicting `snp.pos` (POS = the slot values).
Proposed change: when `snp.pos`/`snp.chr` are supplied explicitly, use them
(overwriting the slot) or stop on conflict; when the slot is used as found,
say so at `verbose >= 1`. Pairs naturally with running
`gl.compliance.check` (post-#330 it clears provably stale copies —
verified: it nulls the platypus slot).

**F2 [HIGH, confidence: high] — REF collapses to N for loci without an observed reference allele (no catalogue rule: target-format fidelity)**
`R/gl2vcf.r:184-194,215` — `mylist.txt` passes only the second `loc.all`
allele via `--reference-allele` (PLINK's A1); the REF column comes from
PLINK's inferred A2, and the `--a2-allele` route that would pin REF from
`loc.all` sits commented out (218-219).
Failure scenario: verified — locus `45063222-23-C/T`, scored only `2` or NA
in the fixture, emitted `REF=N ALT=T` although `loc.all` records `C/T`.
Every locus fixed for the alternate allele in the exported sample loses its
known reference base; `bcftools norm/merge` against a reference genome then
fails or mismatches on those records.
Proposed change: also write the first `loc.all` allele and pass it with
`--a2-allele` (the commented-out lines), so REF always comes from
`loc.all`.

**F3 [MEDIUM, confidence: high] — verbose = 0 is defeated by gl2plink chatter and raw PLINK output (VRB5)**
`R/gl2vcf.r:178,227-236` — `gl2plink` is called with literal
`verbose = NULL` (so it resolves the global default, 2) and
`system_verbose` wraps PLINK's output in an ungated `message()`.
Failure scenario: verified — at `verbose = 0` a run prints 3 stdout lines
("Starting gl2plink", ...) and ~38 message lines of PLINK banner/log.
Proposed change: pass `verbose = verbose` through, and gate the
`system_verbose` message at `verbose >= 2` (keep the captured output for
error reporting).

**F4 [MEDIUM, confidence: high] — SilicoDArT admitted, fails deep in gl2plink with a cryptic error (DAT7)**
`R/gl2vcf.r:95` — default `accept`; this function is one of the two named
offenders in the DAT7 rule.
Failure scenario: verified — `gl2vcf(testset.gs[1:4, 1:6], ...)` dies with
"number of items to replace is not a multiple of replacement length"
(inside gl2plink's genotype recoding, `loc.all` being absent), after
already resolving paths and starting work.
Proposed change: `accept = "SNP"`.

**F5 [MEDIUM, confidence: high] — a factor-typed snp.pos field exports factor level codes as POS (no catalogue rule; same defect class as gl2fasta F2)**
`R/gl2vcf.r:125` — `as.integer(metrics[[snp.pos]])` on a factor returns
level codes.
Failure scenario: verified — with the position column stored as a factor
(routine when metrics pass through `read.csv`/merges), the emitted POS
values were 1, 3, 4, 5 instead of 2438118, 27775762, ...; no warning, file
well formed.
Proposed change: coerce via `as.integer(as.character(...))` and stop if the
result contains NAs.

**F6 [MEDIUM, confidence: high] — external PLINK dependency unguarded; default path is almost never right (DEP1 analogue)**
`R/gl2vcf.r:70,196-239` — no existence check on the binary;
`plink.bin.path` defaults to `getwd()`; the ad-hoc `system_verbose` wrapper
duplicates the existing `utils.plink.run` helper (used by `gl.read.PLINK`).
Failure scenario: verified — with no PLINK at the path the run fails only
after `gl2plink` has written its temp files, with
`'<path>/plink' not found`; a user following the docs from a default call
gets this instead of an instruction to set `plink.bin.path`. The roxygen
default annotation also carries a typo: "[default getwd())]".
Proposed change: check for `plink`/`plink.exe` at `plink.bin.path` up front
and stop with the download URL (DEP1 wording); prefer `utils.plink.run`;
fix the doc typo.

**F7 [LOW, confidence: high] — wrong-length @position silently zeroed (VRB4 proposed)**
`R/gl2vcf.r:105-108` — `length(x$position) != nLoc(x)` routes into the fill
branch; with no `snp.pos` the user's partial coordinates become all zeros
without a message even at `verbose = 2` (verified).
Proposed change: warn at `verbose >= 1` when discarding a malformed slot.

**F8 [LOW, confidence: high] — outpath pollution and a fixed temp-file prefix (FS7 adjacent)**
`R/gl2vcf.r:168-181` — `gl_plink_temp.map/.ped`, `<outfile>.log` and
`<outfile>.nosex` remain in `outpath` after a successful run (verified
listing); the hard-coded `gl_plink_temp` prefix collides across concurrent
sessions sharing an outpath.
Proposed change: write intermediates to `tempdir()` (as `mylist.txt`
already is), clean up PLINK by-products, or document them.

**F9 [LOW, confidence: high] — roxygen drift (DOC1; DOC2, DOC7 proposed rules)**
`R/gl2vcf.r:31-33,49-50` — outdated `verbose` default clause (DOC2);
`@author` names only a custodian (DOC7); "This function requires to
download" grammar; `@details` does not state that CHROM/POS are written as
0 when no genome coordinates are available.
Proposed change: doc pass + `devtools::document()` (DOC4).

## Proposed changes

1. Make explicit `snp.pos`/`snp.chr` take precedence over (or conflict-check
   against) a populated `@position`/`@chromosome`, and announce which source
   was used at `verbose >= 1` (F1).
   **Consequence: numerical output changes (POS column) for callers passing
   snp.pos on objects with a populated position slot.**
2. Pin REF from `loc.all` via `--a2-allele` alongside the existing
   `--reference-allele` list (F2).
   **Consequence: REF changes from N to the recorded reference base for
   alt-fixed loci.**
3. Pass `verbose` through to `gl2plink` and gate PLINK's captured output
   (F3).
4. Restrict the datatype gate with `accept = "SNP"` (F4).
5. Robust position coercion: `as.integer(as.character(...))` + NA check
   (F5).
6. Up-front PLINK binary guard with informative stop; prefer
   `utils.plink.run`; fix the default-annotation typo (F6).
7. Warn when a malformed `@position` slot is discarded (F7).
8. Move intermediates to `tempdir()` and tidy PLINK by-products (F8).
9. Roxygen conformance pass + `devtools::document()` (F9).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (n/a — no plot), STY — run
- Spec: end-to-end runs against PLINK 1.9 (v1.90b7.2 win64) on platypus.gl
  subsets — header validity, CHROM/POS sourcing (slot vs snp.pos vs
  defaults), REF/ALT vs `loc.all`, GT coding vs dosage matrix, missing
  data; stale-slot, conflicting-slot, factor-field, wrong-length-slot and
  missing-binary paths — run empirically
- VCF validation with an external validator (vcftools/GATK): SKIPPED — not
  installed; header and column checks done manually against the 4.2 spec
- `chr.format = "numeric"` path: SKIPPED — platypus chromosome names are
  non-numeric; behaviour inherited from gl2plink and not exercised
- FBM path (DAT6): SKIPPED — no FBM fixture for converters

## Approval (Phase B)

Blanket approval by Arthur Georges, 2026-09-05, via the formal approval
boxes: all findings at every severity approved, explicitly acknowledging
that converted outputs change where they were wrong (POS precedence, REF
no longer degrading to N); the DAT7 class fix (fatal `accept = "SNP"`
gate wherever SilicoDArT is silently admitted) is approved campaign-wide
and applies here as change 4. The genome-only convention (PR #330)
stands: `@position`/`@chromosome` are genome coordinates only, NULL
until assigned; when NULL and no `snp.pos` is supplied, coordinates are
written as 0 — no tag offsets or fabricated coordinates.

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | Approved | Arthur Georges | 2026-09-05; output-change consequence acknowledged |
| 2 | Approved | Arthur Georges | 2026-09-05; output-change consequence acknowledged |
| 3 | Approved | Arthur Georges | 2026-09-05 |
| 4 | Approved | Arthur Georges | 2026-09-05; DAT7 fatal gate |
| 5 | Approved | Arthur Georges | 2026-09-05 |
| 6 | Approved | Arthur Georges | 2026-09-05 |
| 7 | Approved | Arthur Georges | 2026-09-05 |
| 8 | Approved | Arthur Georges | 2026-09-05 |
| 9 | Approved | Arthur Georges | 2026-09-05 |

## Outcome (Phase C)

Applied 2026-09-05 on branch `review-gl2vcf`; PR #352 (base upstream/dev,
ddaed27 — which predates the PR #330 `is.null` test; the rewritten
coordinate block below embeds the genome-only convention either way).
All nine changes applied:

1. F1 — coordinate block restructured: an explicit `snp.pos`/`snp.chr`
   always wins (overwriting the slot); otherwise a valid slot is used as
   found; otherwise coordinates are zeroed. The source used is announced
   at `verbose >= 1`.
2. F2 — REF pinned from `loc.all`. Mechanism deviation, documented:
   PLINK 1.9 rejects `--reference-allele`/`--a1-allele` together with
   `--a2-allele` in one invocation (verified: "Error: --a2-allele cannot
   be used with --a1-allele." — the reason the lines were commented
   out), so the alleles are pinned in two passes: A1 (= ALT) forced from
   `loc.all` while converting ped -> bed, then A2 (= REF) forced while
   recoding bed -> VCF. Outcome verified: REF/ALT match `loc.all` for
   every locus, including alt-fixed loci that previously emitted REF=N.
3. F3 — `verbose` passed through to `gl2plink`; PLINK's captured log
   prints only at `verbose >= 2`, and PLINK's stderr stream (which
   bypasses R's capture) is suppressed below that level via
   `ignore.stderr`.
4. F4 — `utils.check.datatype(x, accept = "SNP", ...)`.
5. F5 — positions coerced `as.integer(as.character(...))`; any NA in
   the result is a fatal error naming the field.
6. F6 — partial, documented: the up-front binary guard (stop with the
   download URL when neither `plink` nor `plink.exe` exists at
   `plink.bin.path`) and the `[default getwd())]` typo fix are applied;
   `utils.plink.run` was NOT adopted because it streams PLINK output
   ungated via `system()` and `setwd()`s into the data directory, which
   would defeat the F3 gating — consolidation left to the custodian.
7. F7 — a wrong-length `@position`/`@chromosome` slot is discarded with
   a warning at `verbose >= 1`.
8. F8 — ped/map (and the new bed) intermediates are written to
   `tempdir()` with PLINK by-products (.log, .nosex) and are removed on
   success (left in place on failure for diagnosis); `outpath` receives
   only the requested `.vcf`.
9. F9 — roxygen: DOC2 verbose clause, DOC7 author labels, grammar fix,
   `@details` states the coordinate precedence and the CHROM/POS = 0
   convention, house tag order; `devtools::document()` run (gl2vcf.Rd
   only).

Verification: baseline characterization test passes (47 assertions,
including the end-to-end block against PLINK 1.9 v1.90b7.2 win64 via
PLINK19_DIR). Every flipped expectation maps to an approved finding:
explicit-argument precedence and slot-as-found (F1), REF = `loc.all`
reference for all loci (F2), `verbose = 0` silence in stdout and
messages (F3), SilicoDArT fatal (F4), factor labels not level codes and
non-numeric field fatal (F5), missing-binary fatal with URL (F6),
wrong-length-slot warning (F7), intermediates/by-products absent from
outpath (F8). End-to-end run on testset.gl at `verbose = 3` clean: 755
body records = nLoc, outpath contains only the VCF. Sibling grep across
the 8 clones: one caller, `dartR.base/R/gl.impute.r:421` (beagle path)
— it sets a valid-length `@position` (1:nLoc) that the new precedence
uses as found (previously zeroed for objects under 1000 loci) and reads
only the `.vcf` from its outpath, so it does not break; the change
repairs its small-dataset path. NEWS.md entry added.

```json
{
  "function": "gl2vcf",
  "package": "dartR.base",
  "family": "io",
  "skill_version": "1.0.0",
  "commit": "ed99203482f4cc51dd2cef70f07493200e534ee7",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "none:target-format-fidelity", "status": "applied", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "applied", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "none:metadata-type-coercion", "status": "applied", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "applied", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB4", "status": "applied", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "FS7", "status": "applied", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 9}
  ],
  "coverage_skipped": ["external VCF validator: not installed", "chr.format numeric path: no numeric-chromosome fixture", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": 352
}
```
