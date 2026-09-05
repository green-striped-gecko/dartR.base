# Review: gl.pcoa (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code)
- Skill: dartr-function-review v2.0.0 (verify-then-apply pass: Phase A
  verification + Phase C apply)
- Package commit: ddaed27 (upstream/dev)
- Branch: review-gl.pcoa
- Datasets: testset.gl, crafted dist matrices (10-, 4- and 3-entity),
  testset.gl[1:4,] low-rank subset, gl.gen2fbm(testset.gl) FBM fixture
- Family mode: analysis
- **Unusual provenance**: findings F1-F6 originated from an EXTERNAL
  session's code-read (no empirical runs, no baseline test) and were
  handed over by the custodian (Arthur Georges, 2026-09-06) with
  decisions already made: apply 1, 3, 4, 5, 6; finding 2 reserved for the
  custodian's own numerical verification; Tracy-Widom flagged as
  unverified with no action. Because the findings arrived unverified,
  every approved finding was reproduced empirically in this session
  before being applied. All reproduced; all five were applied.

### Base reconciliation

The handoff reviewed the `integration-local` copy of `R/gl.pcoa.r`,
reported to carry two tweaks (FBM imputation frequency -> neighbour;
NaN-zeroing of big_SVD trailing singular values) relative to the
2026-08-16 state. Verified: `git diff upstream/dev integration-local --
R/gl.pcoa.r` is empty, and both tweak commits (23f8edf, 4c95bc5) are
ancestors of upstream/dev ddaed27. The reviewed copy is byte-identical to
upstream/dev, so the branch bases on upstream/dev normally, with the
tweaks present. All handoff line numbers re-anchored cleanly.

## Verdicts

- **Standards: Needs work** — the entry scaffold (verbosity, flag
  start/end, datatype check, wd check) conforms; seven standards-walk
  notes remain open (worst MEDIUM), chiefly plots printed regardless of
  `plot.out`/`verbose` and a full densification of FBM objects.
- **Spec: Needs work** — at baseline, `nfactors` was silently ignored on
  the FBM path and the dist path crashed below `nfactors + 1` entities;
  both fixed here. `@return` still names the wrong class and the FBM
  rescaling remains custodian-reserved.

## Findings (handed over; verified this session)

**F1 [HIGH, confidence: high] — `nfactors` silently ignored on the FBM
path (spec/API)**
`R/gl.pcoa.r:774-794` — `glPca(x, nf = nfactors)` truncates, but the
big_SVD branch always used `k = nInd(x) - 1` and returned untruncated
scores/loadings.
Reproduced: `gl.pcoa(gl.gen2fbm(testset.gl), nfactors = 5)` returned
scores 274 x 273 and loadings 611 x 273.
Applied: scores/loadings truncated to `min(nfactors, ncol(dummy$u))`
columns with a `verbose >= 2` warning when clamped; `$eig` kept full
length (matching glPca, and required by the `pc.select` criteria).
Post-fix: scores 274 x 5, loadings 611 x 5, `$eig` length 273,
eigenvalues unchanged.

**F2 [custodian-reserved — NOT touched] — FBM reconstruction rescaling**
`R/gl.pcoa.r:786` (`scores = u %*% diag(d)/2`, `eig = d^2/(4*nInd(x))`,
`loadings = v*2`, marked `#!# intermediate fbm fix`) — internally
self-consistent but numerically unverified against the glPca path
(observed: FBM eig[1] 2.2597 vs glPca eig[1] 1.3822 on testset.gl).
Per the custodian's decision this triple is his to verify; the line is
untouched and carries no added comments. **Open item.**

**F3 [MEDIUM, confidence: high] — second algorithm undocumented (DOC5)**
The roxygen `@details` described the function only as a wrapper for
glPca/pcoa; the big_SVD/FBM branch (with neighbour imputation) was
absent. Verified by reading the header. Applied: a paragraph added to
`@details`; `man/gl.pcoa.Rd` regenerated.

**F4 [LOW, confidence: high] — dead assignment (STY1)**
`R/gl.pcoa.r:848` — `e <- pca$eig[pca$eig > sum(pca$eig /
length(pca$eig))]` immediately overwritten by `e <- eig.top`. Verified by
reading; deletion cannot change behaviour (confirmed by identical
baseline pins post-fix). Applied: line deleted.

**F5 [MEDIUM, confidence: high] — NA lines in `verbose >= 3` output
(DOC5/VRB)**
SNP branch (`R/gl.pcoa.r:851-855`): no length guards on `e[1]+e[2]` /
`e[1]+e[2]+e[3]`. Reproduced: `gl.pcoa(testset.gl[1:4,], verbose = 3)`
printed "PCA Axis 1-3 combined explain NA % of the total variance".
Dist branch (`:727,733`): guard `if(length(eig.top >= 2))` takes
`length()` of a logical vector — truthy whenever `eig.top` is non-empty
(precedence bug). Reproduced: 3-entity dist at `verbose = 3` printed the
Axis 1-3 NA line although only 1 informative and 2 positive axes existed.
Applied: SNP branch gains `if (length(e) >= 2)` / `>= 3` guards; dist
branch guards test `length(eig.raw.pos.pc) >= 2` / `>= 3` — the vector
actually indexed there. Post-fix: no NA lines; combination lines print
only when enough axes exist.

**F6 [HIGH, confidence: high] — crash on small distance matrices (spec)**
`R/gl.pcoa.r:744,746` — `pco$vectors[, 1:nfactors]` with fewer than
`nfactors + 1` entities; the `nInd(x) < 2` guard covers only genlight
input. Reproduced: 4-entity dist with `nfactors = 5` (default) died with
"subscript out of bounds". Applied: `nf.use <- min(nfactors,
ncol(pco$vectors))` before both subsets, `verbose >= 2` warning when
clamped — same pattern as F1. Post-fix: returns scores 4 x 3 with the
warning. (Note: with `correction = "none"` ape::pcoa returns no
`vectors.cor`, so `$loadings` is NULL on this path — pre-existing,
pinned, unchanged.)

**Tracy-Widom criterion** — flagged in the handoff as unverified
(`pc.select = "Tracy-Widom"`, `twtest`/`tw.statistics`). No action taken;
no verification attempted here. **Open item.**

## Standards-walk notes (this session; report-only, NOT applied)

Scope was the custodian's five approved fixes; these are recorded for a
future pass. None reaches BLOCKER/HIGH.

- **N1 [MEDIUM] (VRB5, PLT3)** — the SNP branch calls
  `bs.statistics(eig.raw.pos, plot = T)`, printing a diagnostic ggplot
  regardless of `plot.out` and at `verbose = 0`; `tw.statistics` likewise
  defaults `plot = TRUE` when `pc.select = "Tracy-Widom"`. (The dist
  branch calls `bs.statistics` with plot off.)
- **N2 [MEDIUM] (DAT6)** — the FBM branch tests
  `is.na(sum(as.matrix(x)))`, densifying the entire file-backed matrix
  merely to detect missing values.
- **N3 [MEDIUM] (DOC5)** — `@return` says "An object of class pcoa"; the
  function returns class `glPca` (and the details text says so).
- **N4 [LOW] (VRB5)** — the default `plot.colors = gl.colors(2)` prints
  "Starting gl.colors / Selected color type 2 / Completed: gl.colors" on
  every call, including `verbose = 0` (root cause in `gl.colors`'s
  default verbosity; observed in every run this session).
- **N5 [LOW] (VRB5)** — on the FBM path at `verbose = 0`, the
  `gl.impute(method = "neighbour")` chain prints "Calculating the
  unscaled distance matrix -- euclidean" (cross-function; observed).
- **N6 [LOW] (VRB2)** — two `stop("Error: incorrect specification...")`
  calls lack the `error()` wrapper; each is followed by an unreachable
  `pc.select <- "Tracy-Widom"`.
- **N7 [INFO] (STY1, DOC6)** — ~50 lines of commented-out superseded
  `bs.statistics`, a duplicated comment line, and a stray `#'` mid-
  sentence in the details text ("will not #' exactly represent").

## Proposed changes

1. Truncate FBM-path scores/loadings to `nfactors`, warn when clamped,
   keep `$eig` full length (F1).
2. Document the big_SVD/FBM algorithm in `@details` (F3).
3. Delete the dead `e <-` assignment (F4).
4. Guard the axis-combination `verbose >= 3` lines in both branches;
   fix the dist-branch precedence bug (F5).
5. Clamp `nfactors` to available axes on the dist path, warn when
   clamped (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (notes N1-N7)
- glPca path on testset.gl (nfactors 2/5/10, pinned eigenvalues and
  score mass) — run
- Low-rank verbose path (testset.gl[1:4,]) — run
- dist path: 10-, 4- and 3-entity matrices, `correction = "none"` — run
- FBM path via gl.gen2fbm(testset.gl) (bigstatsr 1.6.x installed) — run
- dist path with correction on small matrices: `correction = "cailliez"`
  on 4 entities fails inside `ape::pcoa` ("arguments imply differing
  number of rows: 2, 0") before gl.pcoa's subsetting — out of scope,
  not pursued
- Tracy-Widom path: SKIPPED — flagged unverified by the custodian, no
  action mandated
- SilicoDArT and fd input paths: SKIPPED — untouched by all five fixes
  (shared code verified via the SNP-branch tests)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 (F1) | approved | Arthur | handed over 2026-09-06; verified here |
| — (F2) | reserved | Arthur | custodian will verify the rescaling numerically himself; not touched |
| 2 (F3) | approved | Arthur | verified here |
| 3 (F4) | approved | Arthur | verified here |
| 4 (F5) | approved | Arthur | verified here |
| 5 (F6) | approved | Arthur | verified here |
| — (Tracy-Widom) | no action | Arthur | flagged unverified |

## Outcome

- All five approved findings reproduced empirically before application;
  all five applied on branch review-gl.pcoa.
- Characterization test `tests/testthat/test-gl.pcoa.R` (new file): at
  baseline, 7 failures, every one at an assertion annotated
  `[approved F1/F5/F6]`; post-apply, all 74 assertions pass. F3/F4
  produce no behavioural diff (F3 docs-only; F4 dead code).
- glPca-path numerical output identical to baseline (pinned eigenvalues
  1.382214/1.344892/1.114741 and score mass 1186.296 unchanged); FBM
  eigenvalues unchanged (2.259724/2.184270/1.652706), only score/loading
  truncation differs.
- E2E: defaults + `verbose = 3` on testset.gl byte-identical verbose
  output; 4-entity dist now clamps with a warning instead of crashing.
- Caller grep across all 8 clones: `gl.assign.mahal`/`gl.assign.
  mahalanobis`/`gl.assign.pca` (dartR.captive, dartR.popgen) and
  `gl.pcoa.plot` — all genlight-path callers with small `nfactors`;
  no caller reaches the changed edge cases; all clear.
- `man/gl.pcoa.Rd` regenerated (F3); NEWS.md entry added.
- Open custodian items: F2 rescaling triple; Tracy-Widom criterion.
- PR: pending (number added in follow-up commit)

```json
{
  "function": "gl.pcoa",
  "package": "dartR.base",
  "family_mode": "analysis",
  "commit": "ddaed27",
  "skill_version": "2.0.0",
  "provenance": "external code-read handed over by custodian 2026-09-06; verified empirically this session",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec/API1", "status": "approved", "reproduced": true, "applied": true, "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "medium", "rule": "spec", "status": "custodian-reserved", "reproduced": null, "applied": false},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "reproduced": true, "applied": true, "change": 2},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "reproduced": true, "applied": true, "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5/VRB", "status": "approved", "reproduced": true, "applied": true, "change": 4},
    {"id": "F6", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "reproduced": true, "applied": true, "change": 5},
    {"id": "N1", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5/PLT3", "status": "noted"},
    {"id": "N2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT6", "status": "noted"},
    {"id": "N3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "noted"},
    {"id": "N4", "severity": "LOW", "confidence": "high", "rule": "VRB5", "status": "noted"},
    {"id": "N5", "severity": "LOW", "confidence": "high", "rule": "VRB5", "status": "noted"},
    {"id": "N6", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "noted"},
    {"id": "N7", "severity": "INFO", "confidence": "high", "rule": "STY1/DOC6", "status": "noted"}
  ],
  "open_items": ["F2 FBM rescaling triple (custodian)", "Tracy-Widom criterion unverified"],
  "coverage_skipped": ["Tracy-Widom path (custodian-flagged, no action)", "SilicoDArT/fd input (untouched by fixes)", "dist+correction small-matrix ape failure (out of scope)"],
  "datasets": ["testset.gl", "crafted dist 10/4/3", "testset.gl[1:4,]", "gl.gen2fbm(testset.gl)"],
  "baseline_test": "tests/testthat/test-gl.pcoa.R",
  "status": "pr-open",
  "pr": null
}
```
