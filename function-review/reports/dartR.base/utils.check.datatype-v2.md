# Review: utils.check.datatype (dartR.base) — second pass (v2)

## Provenance

- Model: Claude Fable 5 (via dartr-dev agent); Skill: dartr-function-review
  v2.0.0; Package commit reviewed: ddaed27 (upstream/dev), exported to an
  isolated tree via `git archive` because the working checkout
  (integration-local) already carries the PR #296 fix; Date: 2026-09-06.
- Prior review: `utils.check.datatype.md` (v1) — findings F1–F5 applied on
  branch `review-utils.check.datatype`, PR #296 open. This pass reviews the
  same upstream state to the standard of the package's single most-called
  function, verifies the PR #367 lead the v1 review did not record, and
  separates what #296 fixes from what remains.
- Family mode: analysis/utility (internal dispatcher; its verdict steers
  SNP/SilicoDArT dispatch in every caller, including the converter
  `accept =` gates added by the current campaign).
- Datasets: testset.gl, testset.gs, platypus.gl (dartR.data 1.2.5);
  constructed fixtures: all-all-NA dartR object, individual-only all-NA,
  bare/flagless/contradictory genlights, mixed ploidy, zero-locus,
  zero-individual (slot surgery), FBM-backed (`gl.gen2fbm`), fd, glPca,
  data.frame, matrix, dist, list, NULL, atomic vectors, unknown S3 class.
- Baseline: `tests/testthat/test-utils.check.datatype-edges.R` (new file;
  the existing `test-utils.check.datatype.R` is the PR #296 baseline and is
  left untouched). 37 assertions, all passing on integration-local.
- Checks skipped: Google Group / GitHub issue search (not available: no
  browser session — but the PR #367 review's lead is incorporated and
  verified). Everything else run.

## Verdicts (for the reviewed state, upstream/dev = ddaed27)

**Standards: Needs work** — the structure conforms (verbosity resolution,
`stop(error())` idiom, invisible return, read-only), but a message-gated
code path runs a full `gl.filter.allna` filter on every function entry at
default verbosity, one warning prints at `verbose = 0`, and the ploidy
branch carries unreachable dead code.

**Spec: Needs work** — classification is a pure ploidy-slot lookup that
works as documented for well-formed objects, but at `verbose >= 2` the
courtesy all-NA check can kill the caller's run outright on degenerate
objects (the PR #367 lead, confirmed), and the individuals warning carries
the loci wording.

Open PR #296 already fixes the four worst findings (F1, F3, F4, F5 below);
merging it is the single action that moves both verdicts to Ready apart
from the residuals listed at F2 and F6–F9.

## The detection contract (established empirically)

Classification consults **only the `@ploidy` slot** — never genotype
values, dartR flags, or `loc.metrics`:

| Input | Result (verbose 0/1) | Result (verbose >= 2) |
|---|---|---|
| genlight/dartR, all ploidy 1 | `"SilicoDArT"` | same + report line |
| genlight/dartR, all ploidy 2 | `"SNP"` | same + report line |
| mixed/other ploidy (1,2 / 4) | `"SNP"` silently | ddaed27: `"SNP"`, no notice at any verbosity; #296: warns |
| ploidy slot NULL/empty (bare `new("genlight")`, zero-individual object) | fatal "ploidy not set … run gl.compliance.check" | same |
| flagless genlight from a matrix (no dartR metadata) | classified normally | same |
| silico 0/1 values stamped ploidy 2 | `"SNP"` (passes `accept="SNP"` gates) | same |
| 0/1/2 values stamped ploidy 1 | `"SilicoDArT"` | same |
| all loci all-NA (dartR class) | `"SNP"` | ddaed27: **ERROR** "Subsetting resulted in zero loci."; #296: warns |
| zero-locus (plain genlight) | `"SNP"` | ddaed27: **ERROR** "subscript out of bounds"; #296: completes |
| FBM-backed dartR (`@gen` empty, ploidy populated) | `"SNP"` | same, after a full `as.matrix` densification |
| fd (with genlight `$gl`) | `"fd"` | same + report |
| fd with non-genlight `$gl` | fatal "Fixed Difference object expected!" | same |
| dist / matrix / glPca / list | that name | same + report |
| data.frame | `"list"` (via `is.list`) | same |
| NULL, atomic vectors, genind, anything else | `class(x)[1]`, then the accept gate (fatal for the default accept) | same; ddaed27 prints the class warning even at verbose 0 |

`accept` matching is exact `%in%` — case-sensitive ("snp" and "silicodart"
never match), no partial matching. On ddaed27 the default entries
"genlight" and "dartR" are dead (datatype is never literally either), so
`accept = "genlight"` alone rejects every genlight; #296 makes them admit
both genotype datatypes unless a specific datatype is listed (which then
governs). The fatal reads:
`Fatal Error: inappropriate object passed to function, found <datatype>
expecting <accept, " or "-collapsed>` (no trailing newline).

The function is read-only (input `identical()` before/after at verbose 3,
verified) and returns the datatype string invisibly. `verbose = 0` and
`verbose = 1` are byte-identical (no begin/end banner — the house pattern
for this util; `utils.flag.start` is the caller's job).

## Findings

**F1 [HIGH, confidence: high] — verbosity flag changes control flow: the
courtesy all-NA check can crash the caller's preamble (VRB5, DOC5; defect
class: plot/message code coupled to results, cf. PLT3)**
`R/utils.check.datatype.R:87-103` (ddaed27) — at `verbose > 1` the check
calls `gl.filter.allna(x, verbose = 0)`. When every locus is all-NA,
`gl.filter.allna`'s subset `x[, !x$loc.names %in% loc.list]`
(`R/gl.filter.allna.r:117`) hits the dartR `[` guard
(`R/utils.dartR.class.def.r:366`) and the run dies with the contextless
message "Subsetting resulted in zero loci.". A zero-locus input dies
earlier, in `gl.filter.allna`'s `for (i in 1:nL)` loop with "subscript out
of bounds". Failure scenario: default verbosity is 2
(`gl.check.verbosity(NULL)` → 2), so ANY `gl.*` call on such an object —
e.g. a population subset in which no locus remains scored — crashes in the
preamble before the function proper (or its own `accept` gate: verified
identical with `accept = "SNP"`) ever sees the object, while the same call
at `verbose = 0/1` proceeds. Measured cost when it does not crash: ~87 ms
and two full `as.matrix` densifications per function entry on the 250x755
testset. This is the PR #367 lead, confirmed and root-caused.
Proposed change: none new — PR #296's replacement (direct
`colSums`/`rowSums` over `is.na(as.matrix(x))`) removes the filter call and
with it every crash path; verified no-crash on integration-local. The v1
review/#296 did not record the crash, so the fix was accidental and
unpinned — the new baseline `test-utils.check.datatype-edges.R` now pins
it against regression.

**F2 [MEDIUM, confidence: high] — the courtesy scan densifies the full
matrix on every function entry at default verbosity (DAT6 proposed, STY2)**
`R/utils.check.datatype.R:87-103` — the scan is O(nInd x nLoc) with a full
`as.matrix(x)` (twice on ddaed27, once under #296: ~87 ms → ~24 ms on
250x755). FBM-backed objects (`@gen` empty, verified via `gl.gen2fbm`)
classify correctly, but at `verbose >= 2` every function entry
materialises the entire FBM as a dense matrix — the exact pattern DAT6
exists to prevent, paid once per `gl.*` call. Failure scenario: a
50k-loci x 5k-individual FBM dataset densifies (~2 GB) on every call at
default verbosity. Proposed change (residual after #296; echoes the v1
report's own note): consult `loc.metrics.flags$allna` where fresh, or scan
column-wise via FBM accessors; a team-level design decision, recorded not
actioned.

**F3 [MEDIUM, confidence: high] — individuals-all-NA warning carries the
loci wording (DOC5, VRB2)**
`R/utils.check.datatype.R:96-102` (ddaed27) — the `nInd(tmp) < nInd(x)`
branch prints "loci that are scored NA across all individuals". Verified:
an object whose only defect is one all-NA individual warns about loci.
Fixed by #296 (its F1); pinned by both test files.

**F4 [MEDIUM, confidence: high] — `accept = "genlight"`/`"dartR"` alone
rejects every genlight (DAT7)**
`R/utils.check.datatype.R:161` (ddaed27) — datatype is never literally
"genlight"/"dartR", so those accept entries are dead: verified
`accept = "genlight"` → "Fatal Error: … found SNP expecting genlight".
Latent, not live: the family census (below) shows no caller passes them
without a specific datatype. Fixed by #296 (its F3, with the
specific-listing-governs clause — verified `c("genlight","SNP")` remains
SNP-only on integration-local).

**F5 [LOW, confidence: high] — unknown-class warning prints at
`verbose = 0` (VRB5, VRB3)**
`R/utils.check.datatype.R:155` (ddaed27) — the fallback
`cat(warn("Found object of class …"))` is ungated; verified it prints at
verbose 0 (including when the unknown class is accepted and the call
succeeds). Fixed by #296 (its F4, gated at >= 1).

**F6 [LOW, confidence: high, proposed rule] — classification-affecting
ploidy notice absent below verbose 2 (VRB4)**
`R/utils.check.datatype.R:73-85` (ddaed27) — mixed/abnormal ploidy is
silently classified "SNP" with no notice at ANY verbosity, and the
`else { stop("SNP or SilicoDArT coding misspecified") }` at :78-85 is
unreachable (the NULL-ploidy case already stopped at :53), so the
documented misspecification error can never fire. #296 removes the dead
code and warns at `verbose >= 2` — but under VRB4 [proposed] a warning
that a classification decision was made for the user belongs at
`verbose >= 1`. Residual; docs-level severity.

**F7 [LOW, confidence: high] — contract not documented: detection is
ploidy-slot-only; contradictory content passes the campaign's accept gates
(DAT5, DOC5)**
Whole genlight branch — genotype values are never examined: a 0/1-valued
(SilicoDArT-style) object stamped ploidy 2 returns "SNP" and passes
`accept = "SNP"` into dosage-based converter math; 0/1/2 values stamped
ploidy 1 return "SilicoDArT" (both verified). A zero-individual object
exits via the misleading "ploidy not set" fatal (the ploidy accessor
returns NULL for an empty slot). Not a code defect — the speed of a pure
slot lookup is why this function can run on every entry — but the roxygen
("checks … whether the genlight object is a SNP dataset or a SilicoDArT
object") implies inspection it does not do, and the campaign's converter
gates inherit the limitation. Proposed change (docs-only): state in
`@details` that classification is by the ploidy slot alone and content is
not verified (`gl.compliance.check` is the content-level check). Residual
after #296.

**F8 [LOW, confidence: high] — a data.frame is classified "list" (DOC5)**
`R/utils.check.datatype.R:149` — `is.list()` catches data.frames before
any data.frame branch; verified "found list expecting matrix" when a
caller gates on `accept = "matrix"`. Confusing but fail-safe; a behaviour
change would touch the 4 callers that gate on "list" (API1), so the
proposed change is docs-only: note the mapping in `@return`. Residual
after #296; pinned as BUG in the new baseline.

**F9 [INFO, confidence: high, proposed rule] — `@author` names a Custodian
only, no `Author(s):` line (DOC7)**
Header — residual after #296. Proposed change: add
`Author(s): Arthur Georges.` per the DOC7 default.

Notes, not findings: `accept` matching is case-sensitive and exact
(undocumented; no caller misspells — census below); the accept-gate fatal
lacks the house trailing `\n`; `is(x,"dartR") | is(x,"genlight")` uses
vectorised `|` where `||` is meant (harmless on scalars); the ddaed27
header defects (accept default typo, 'mat'/'glPCA', stale
gl.check/fd.check note, `@return` after `@export`) are all fixed by #296
(its F5). Scope-rule notes on OTHER functions: `gl.filter.allna` cannot
return a zero-locus/zero-individual result (`[` guard) and its `1:nLoc`
loop breaks on zero-locus input — surfaced here only through the courtesy
call; recorded, not fixed. `testset.gl` (dartR.data) itself contains 3
all-NA loci, so every `gl.*` call on the standard fixture at default
verbosity prints the all-NA warning.

## Caller census (all 8 clones, grep of R/ trees)

Call-site lines containing `utils.check.datatype` (includes the definition
file and roxygen mentions in dartR.base): dartR.base 148, dartR.captive
35, dartR.popgen 24, dartR.sim 1, dartR.spatial 5, dartR.sexlinked 3,
dartRstartup 1, dartRverse 0 — 217 total.

`accept =` patterns (paren-aware extraction of actual calls):

| Pattern | base | captive | popgen | sim | spatial | sexlinked | startup |
|---|---|---|---|---|---|---|---|
| default (omitted) | 104 | 32 | 18 | 1 | 5 | 0 | 1 |
| `"SNP"` | 18 | 0 | 5 | 0 | 0 | 0 | 0 |
| `c("genlight","SNP")` | 5 | 0 | 0 | 0 | 0 | 0 | 0 |
| `c("SNP","SilicoDArT")` | 4 | 0 | 0 | 0 | 0 | 0 | 0 |
| `c("genlight","SNP","SilicoDArT")` | 3 | 2 | 0 | 0 | 0 | 0 | 0 |
| `"SilicoDArT"` | 2 | 1 | 0 | 0 | 0 | 0 | 0 |
| `c("dartR","genlight","SNP","SilicoDArT")` | 0 | 0 | 0 | 0 | 0 | 3 | 0 |
| `"fd"` / fd-containing vectors | 5 | 0 | 1 | 0 | 0 | 0 | 0 |
| dist/matrix/glPca/list vectors | 5 | 0 | 0 | 0 | 0 | 0 | 0 |

No caller anywhere passes "genlight"/"dartR" without a specific genotype
datatype (F4 is latent, not live); no caller uses a lowercase or variant
spelling; every `accept = "SNP"` gate the converter campaign added behaves
identically before and after #296 (verified on both trees). No caller
depends on the verbosity-coupled crash, so removing it (merging #296)
breaks nothing — callers only gain the ability to survive degenerate
objects at default verbosity long enough to fail in their own, better
error paths.

## Interaction with campaign PRs

- **PR #296 (this function, open)**: fixes F1 (incidentally — the crash
  was never recorded there; now pinned), F3, F4, F5, most doc defects, and
  adds the F6 warning at `verbose >= 2`. Merging it is the main outcome
  this review confirms. Residuals after merge: F2 (FBM densification /
  per-call scan cost), F6 (warning level, proposed rule), F7–F9
  (docs-level).
- **PR #367 (gl.compliance.check)**: its residual note IS F1; resolved by
  merging #296. The compliance-check path that repairs all-NA data can now
  be reached at default verbosity instead of dying in the preamble.
- **Converter accept-gates (gl2* campaign)**: gates verified working on
  both trees; F7 documents their limit — they guarantee the ploidy stamp,
  not value-level datatype consistency.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY (PLT n/a — no plotting) — run
- Spec axis: behaviour vs roxygen, every branch traced and exercised
  empirically on ddaed27 (isolated `git archive` tree) and cross-checked
  on integration-local (#296) — run
- Degenerate-input matrix (all-NA, zero-dim, bare, flagless,
  contradictory, mixed ploidy, NULL, non-genlight, malformed fd) — run
- FBM path (DAT6): run via `gl.gen2fbm(testset.gl)` (empty `@gen`
  confirmed; classification and verbose-2 scan exercised); large-FBM
  densification cost extrapolated, not measured — no large FBM fixture
- Verbosity contract: verbose 0/1/2/3 captured (sink-based, error-safe)
- Read-only contract: `identical()` before/after — run
- Timing: 20-call means, ddaed27 vs #296 — run
- Caller census: all 8 clones — run
- Google Group / GitHub issues: SKIPPED (no browser session; PR #367 lead
  used instead)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| Amend PR #296: record and pin F1 | Approved 2026-09-06 | Arthur Georges | Addendum + edge baseline pushed to the PR branch |
| F2 FBM-safe courtesy scan | Approved 2026-09-06 | Arthur Georges | Applied |
| F6 mixed-ploidy notice at verbose >= 1 | Approved 2026-09-06 | Arthur Georges | Applied |
| F7 content-vs-ploidy consistency check | Approved 2026-09-06 | Arthur Georges | Applied as a behaviour change (caught at the gate), upgraded from the docs-only proposal |
| F8 data.frame classified "data.frame" | Approved 2026-09-06 | Arthur Georges | Applied as a behaviour change, upgraded from the docs-only proposal |
| F9 Author(s) line | Approved 2026-09-06 | Arthur Georges | Applied |
| Notes: gate-fatal trailing newline, scalar `\|\|`, accept case-sensitivity documented | Approved 2026-09-06 | Arthur Georges | Applied with the LOWs |

## Outcome

Phase C complete (2026-09-06). Two actions:

1. **PR #296 amended** (branch `review-utils.check.datatype`, commit
   e42eca3): the edge baseline committed to `tests/testthat/` (passing,
   37 assertions, pinning the F1 crash fix), an addendum appended to the
   v1 report recording that the F1 crash existed on ddaed27 and is
   resolved by that PR, and the PR body updated with an Addendum section.

2. **Residuals applied** on the stacked branch
   `review-utils.check.datatype-2`:
   - F2: the courtesy all-NA scan reads FBM-backed objects in 1024-column
     blocks via the `@fbm` accessor — no full densification. The
     gen-backed path keeps the single `as.matrix` scan (unchanged cost).
     The `loc.metrics.flags$allna` shortcut was NOT taken (stale-flag
     semantics would change behaviour for clean objects).
   - F6: the mixed-ploidy notice prints from `verbose >= 1` (VRB4).
   - F7: content-vs-ploidy consistency at the gate — an object of uniform
     ploidy 2 whose non-missing genotypes are all 0/1, and which carries
     no SNP metadata (empty `loc.all` slot, no `SNP`/`SnpPosition` locus
     metrics), is fatal at any verbosity (not verbosity-coupled, avoiding
     the F1 defect class), with an actionable message
     (`ploidy(gl) <- rep(1, nInd(gl))`). The metadata clause was added
     during Phase C verification: a blanket no-2s fatal broke two
     committed campaign tests that deliberately protect clean SNP subsets
     with the homozygous-alternate class absent
     (`test-gl.report.basics.R` "SNP subset with a genotype class
     absent", `test-gl.smearplot.R` "subset with no homozygote-alternate
     scores still plots") — violating the approval's own clean-objects-
     unchanged constraint. The `loc.all` slot discriminates cleanly:
     populated with allele pairs on `testset.gl` (and preserved by
     subsetting), empty on `testset.gs` even though its `loc.metrics`
     carries SNP-style recalculated columns. The scan is early-exit
     (first genotype of 2), so clean SNP objects pay for one individual
     (or one FBM column block); all-NA and zero-locus objects are exempt
     (no content to contradict — a zero-locus SNPbin decodes to a
     spurious 0, guarded explicitly). The reverse mismatch (0/1/2
     content stamped ploidy 1) remains uncaught: proving the absence of
     2s in clean SilicoDArT data would cost a full scan on every silico
     function entry at every verbosity. A bare genlight built from a
     0/1 matrix and stamped ploidy 2 with no metadata at all is caught
     (pinned).
   - F8: a data.frame is classified `"data.frame"`; caller check: the
     only `accept`-"list" gates in the family are two sites in
     `gl.pcoa.plot.r`, neither a plausible data.frame recipient.
   - F9: `Author(s): Arthur Georges.` added; Rd regenerated.
   - Notes applied: gate fatal ends with `\n`; `|` → `||`; the
     ploidy-slot contract, accept case-sensitivity and the data.frame
     mapping documented in the header.

Verification: #296 baseline 19/19 unchanged; edge baseline 40/40 (37 + 3
assertions added), flips/extensions confined to pins annotated
`# [approved F2/F7/F8]`. Spot runs: testset.gl/testset.gs preamble
output byte-identical at verbose 0 and 2; testset.gs stamped ploidy 2
now fatal at verbose 0 and 2 (list- and FBM-backed); a testset.gl
subset with no homozygous-alternate scores passes an `accept = "SNP"`
gate; data.frame returns its own name; plain list still "list". Timing
(20-call means, 253x255 testset.gl): verbose 0 0.5 ms vs 2.0 ms (#296),
verbose 2 22.3 ms vs 22.2 ms (noise); FBM verbose 2 5.2 ms vs 3.4 ms —
the block scan trades ~2 ms on the small fixture for not materialising
the full matrix on large FBMs (the DAT6 point). Full-tree regression
sweep: identical failure set with and without the residuals (40 failed
assertions, 8 errors, same 33 test blocks; 832 vs 829 passing — the
difference is the three added edge assertions; the failing set is
environmental — the installed dartR.data 1.2.5 vs expectations
recomputed for CRAN 1.2.2 — plus the pre-existing dev failures already
itemised in the v1 report; details in the PR). The
first sweep iteration caught the blanket-F7 breakage of the two
protected no-2s tests, which drove the metadata refinement above.
Caller safety: census of 217 sites re-confirmed; the 23
`accept = "SNP"` sites gain rejection of mislabelled silico input (the
campaign's gate intent); no caller passes a data.frame to a "list"
gate.

Manifest: new row `pr-open` for `review-utils.check.datatype-2`
appended alongside the untouched 296 row.

```json
{
  "function": "utils.check.datatype",
  "package": "dartR.base",
  "family": "utility",
  "review_pass": 2,
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "VRB5/DOC5", "status": "fixed-by-pr-296", "pinned_by": "test-utils.check.datatype-edges.R"},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "applied", "note": "FBM block scan; flags shortcut rejected (stale-flag semantics)"},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5/VRB2", "status": "fixed-by-pr-296"},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT7", "status": "fixed-by-pr-296"},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB5", "status": "fixed-by-pr-296"},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB4", "proposed_rule": true, "status": "applied"},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5/DOC5", "status": "applied", "note": "upgraded at approval: behaviour change, content-vs-ploidy mismatch fatal at the gate"},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "note": "upgraded at approval: data.frame classified 'data.frame'"},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DOC7", "proposed_rule": true, "status": "applied"}
  ],
  "coverage_skipped": ["Google Group/GitHub search: no browser session", "large-FBM densification cost: no large FBM fixture"],
  "datasets": ["testset.gl", "testset.gs", "platypus.gl", "constructed degenerate fixtures", "gl.gen2fbm(testset.gl)"],
  "baseline_test": "tests/testthat/test-utils.check.datatype-edges.R",
  "status": "applied",
  "pr": 296,
  "pr_followup": "pending",
  "branch_followup": "review-utils.check.datatype-2"
}
```
