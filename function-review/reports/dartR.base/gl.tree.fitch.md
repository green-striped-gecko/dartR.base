# Review: gl.tree.fitch (dartR.base)

- Family mode: analysis
- Date: 2026-09-06
- Reviewer: Claude (Claude Fable 5, via dartr-dev agent), dartr-function-review v2.0.0
- Package commit: ddaed27 (upstream/dev; local HEAD ed99203 — `git diff upstream/dev -- R/gl.tree.fitch.r` is empty, so `load_all()` exercised the reviewed code exactly)
- Datasets: testset.gl (dartR.data 1.2.5; 5 largest populations, `gl.filter.callrate(threshold = 1)`, distances via `gl.dist.phylo(subst.model = "F81")`)
- External toolchain: PHYLIP 3.695 Windows binaries at `D:/workspace/R/phylip-3.695/exe` (`fitch.exe`, `consense.exe`)
- Baseline: `tests/testthat/test-gl.tree.fitch.R` (new file; snapshot captured pre-review, all 7 tests pass)

What the function is: despite the campaign's working assumption of Fitch *parsimony*, `gl.tree.fitch` is a wrapper around the PHYLIP `fitch` binary — **Fitch–Margoliash weighted least-squares distance trees** (power 2, non-negative branch lengths) — taking a precomputed `dist` object (typically from `gl.dist.phylo`), with `consense` supplying majority-rule bootstrap support. Genotypes are never touched directly; the only genotype-side work is the bootstrap resampling of `x` via `gl.subsample.loc` and re-computation of distances.

**Empirical method note.** As shipped, the function cannot complete a run on this Windows machine (F1), so downstream behaviour was exercised through a shimmed copy: the function object verbatim, with its lexical scope given a `system()` binding that routes through `shell()` so the `< fitch.cmd` redirection the code writes actually happens. Every empirical claim below other than F1 was observed through that shim; the shim changes nothing but the redirection. The same shim is used in the baseline test file.

## Verdict

**Standards: Needs work** — the FS preamble is present and the datatype guard works, but the plot-parameter bundle is absent, verbosity gating is widely violated (`verbose = 0` prints the full PHYLIP dialogue), and the error idiom raises empty-message errors.

**Spec: Rework** — the core FM tree is numerically correct (independent exhaustive verification: RF 0), but the function as shipped never gets that far on Windows (F1), and the entire bootstrap-support pipeline delivers wrong numbers on the wrong nodes of the wrong tree and then discards them (F2–F4); that pipeline needs redesign as a unit, not a patch.

What works well: the Fitch–Margoliash topology and the menu-scripting protocol (including the append-prompt "A" answers) are correct when the plumbing runs.

## Independent verification (spec axis, analysis mode)

Reproduced the method independently: all 15 unrooted 5-taxon topologies enumerated; each fitted by weighted least squares with FM weights `1/D_ij^2` (power 2) under non-negativity (L-BFGS-B), from the same `dist` object.

- Topology: fitch's tree matches the exhaustive optimum, **RF = 0**; the optimum is distinct (best SSQ 0.3805 vs next 0.4012).
- Branch lengths: agree only to the 6-decimal precision of PHYLIP's `outtree` — see F10; max cophenetic discrepancy 7.8e-06 against a maximum branch of 5e-05.
- Input fidelity: the `infile` written from the `dist` object is faithful; PHYLIP 3.695 parses R's scientific-notation output (`8.111075e-05`) correctly (verified by branch-length scale of the result).

## Findings

**F1 [BLOCKER, confidence: high] — PHYLIP invocation performs no redirection on Windows (platform-portability principle; no catalogue rule covers platform-specific `system()` calls — gap noted; DOC5 for the failed documented workflow)**
`R/gl.tree.fitch.r:150` (also 155, 160; executed at 285, 388, 412) — `system("fitch.exe < fitch.cmd")`: on Windows, `system()` calls CreateProcess without a shell, so `<` is passed as an argument and the command file is never fed to `fitch`.
Failure scenario (reproduced): under `Rscript`/`R CMD check`, `fitch` reads EOF from stdin, aborts with "Made 100 attempts to read input in loop", no `outtree` is written, and the function dies with the opaque `cannot open the connection`. With an open interactive stdin, `fitch` instead blocks indefinitely waiting for keyboard input (also reproduced). The function has no working path on Windows.
Proposed change: replace all three invocations with `system2(prog, stdin = "fitch.cmd", stdout = ...)` — portable across all three platforms and removing the copied-binary/`cwd` dance dependency on shell semantics.

**F2 [BLOCKER, confidence: high] — bootstrap support values are extracted from the wrong edges and drawn on the wrong tree (numerical-correctness principle, analysis mode; DOC5)**
`R/gl.tree.fitch.r:424-433, 443` — the code assumes `tree.bstraps$edge.length` is ordered "tip edges first, then internal edges" and slices `branch_lengths[(num_tips+1):length(branch_lengths)]`. ape's edge order is cladewise, not tips-first.
Failure scenario (reproduced, `bstrap = 5`, seed 42): true internal-edge supports were (1.0, 0.8, 0.6); the code's slice picked three *tip* edges and displayed (1, 1, 1) — every node labelled full support regardless of the data. Additionally the slice yields `Nnode - 1` labels where `nodelabels()` expects `Nnode` (3 vs 4 observed → recycling/misplacement), and the labels are drawn on `tree_1` (the best tree) using node numbers from the *consensus* tree, two objects with unrelated node numbering. Users read fabricated bootstrap support.
Proposed change: compute supports on `tree_1` directly with `ape::prop.clades(tree_1, <replicate trees>)/bstrap` (reading the replicate trees fitch writes before consense), attach them as `tree_1$node.label`, and plot from `node.label`.

**F3 [HIGH, confidence: high] — bootstrap replicate distances ignore the model that built `D` (DOC5)**
`R/gl.tree.fitch.r:326` — replicates call `gl.dist.phylo(tmp, by.pop = TRUE)` with defaults: `subst.model = "F81"`, `pairwise.missing = TRUE`, `min.tag.len = NULL`. The user's `D` may have been built with any other settings (the docs themselves show `subst.model` as a user choice).
Failure scenario: `D` built with `subst.model = "K80"`; the tree is a K80 tree but its "support" values measure the stability of F81 trees — a different statistic, silently.
Proposed change: accept and forward the distance parameters (explicit `subst.model`/`pairwise.missing`/`min.tag.len` arguments or `...`) to the replicate `gl.dist.phylo` calls.

**F4 [HIGH, confidence: high] — bootstrap results exist only as a (mis)drawn plot; nothing is returned (PLT3)**
`R/gl.tree.fitch.r:416-443, 454` — `tree.bstraps` (the consensus) and `node_values` are never attached to the return value; `tree_1$node.label` is `NULL` (verified empirically). The return is bit-identical to a `bstrap = 1` run.
Failure scenario: any user who computes 1000 bootstraps and wants to report support values, re-plot, or write the tree to a file has nothing — the computation is unrecoverable once the plot window closes.
Proposed change: attach supports as `tree_1$node.label` (per F2) so `write.tree()` round-trips them; optionally return the consensus tree as an attribute.

**F5 [HIGH, confidence: high] — working directory left at `tempdir()` after any bootstrap run (FS7; cwd-hygiene precedent from `gl2paup.parsimony`)**
`R/gl.tree.fitch.r:353-354` — `hold <- getwd(); setwd(tempdir())` opens the bootstrap block; the matching restore is commented out at line 413 and never executed on any path.
Failure scenario (reproduced): after `gl.tree.fitch(..., bstrap = 5)` returns, `getwd()` is `C:/Users/.../Temp/Rtmp...` — every subsequent relative-path operation in the user's session lands in a temporary directory that vanishes at session end. Also a CRAN policy issue.
Proposed change: `on.exit(setwd(hold), add = TRUE)` immediately after the first `setwd()`, covering error paths too (which also fixes the stranded cwd when `readLines`/`read.tree`/`plot` fail mid-run).

**F6 [HIGH, confidence: high] — `out.path` is documented but never used (DOC5; silently-ignored parameter)**
`R/gl.tree.fitch.r:13, 79` — `out.path` ("Path to the directory to save files produced by the analysis [default tempdir()]") appears only in the signature; every file (infile, outfile, outtree, intree, command files, copied executables) goes to `tempdir()` unconditionally, and the closing message says so.
Failure scenario (reproduced): `out.path = <dir>` leaves the directory empty; the user's requested copies of the PHYLIP outputs are destroyed at session end.
Proposed change: copy `infile`/`outfile`/`outtree` to `out.path` when it differs from `tempdir()`; otherwise remove the parameter (removal is an API2 change).

**F7 [MEDIUM, confidence: high] — `verbose = 0` is not silent (VRB5, VRB3, VRB2)**
`R/gl.tree.fitch.r:285, 388, 412` — `system()` without stdout suppression echoes the full PHYLIP menus, progress, and "Press enter to quit" at every verbosity (reproduced at `verbose = 0`). `:111, 121, 132, 137` — `cat(warn(...))` ungated. `:356` — raw ungated `cat("    Deleting old outtree file\n")` (also VRB2: not routed through a message helper).
Failure scenario: scripted/pipeline use at `verbose = 0` floods the console with several screens of PHYLIP dialogue per call (hundreds with bootstrapping).
Proposed change: capture/suppress subprocess output unless `verbose >= 2` (a `system2` migration per F1 makes this a `stdout=` argument); gate the warnings at `verbose >= 1` (VRB4 applies — they change results); delete or gate the raw `cat`.

**F8 [MEDIUM, confidence: high] — no plot bundle; the unconditional plot couples plotting to the result (PLT1, PLT3, VRB5)**
`R/gl.tree.fitch.r:441-444` — `plot()` runs on every call; there is no `plot.display`/`plot.file`/`plot.dir`/`plot.theme`, no `verbose == 0` plot gate, and the plot precedes the `return`.
Failure scenario (reproduced via F9): when `plot()` throws, the already-computed tree is lost — the function errors after the tree exists in a local variable. Also `plot.phylo` warnings show `edge.width`/`edge.color` leaking into `plot.window` for degenerate inputs.
Proposed change: add `plot.display` (house preamble: `if (verbose == 0) plot.display <- FALSE`), plot after the result is secured (or wrap in `tryCatch`), and adopt the `plot.file`/`plot.dir` idiom (PLT2) for saving.

**F9 [MEDIUM, confidence: high] — fewer than four taxa fails opaquely (FS5)**
`R/gl.tree.fitch.r` (no guard exists) — PHYLIP `fitch` writes an *empty* `outtree` for 3 taxa (despite claiming success) and segfaults for 2; `read.tree` then returns `NULL` and the unconditional `plot(NULL)` throws `need finite 'xlim' values` (both reproduced).
Failure scenario: a user with 3 populations gets a graphics error pointing nowhere near the actual cause.
Proposed change: `stop(error(...))` when `length(attr(D, "Labels")) < 4` with a message stating the PHYLIP minimum.

**F10 [MEDIUM, confidence: medium] — branch lengths quantised to PHYLIP's 6-decimal output at SNP-scale distances (numerical-correctness principle; same scaling class as gl.tree.nj F7)**
`R/gl.tree.fitch.r:227-238, 299` — F81 distances from SNP data are ~1e-05 here; `outtree` carries 6 decimals, so branch lengths survive with 1–2 significant digits.
Failure scenario (measured): true branch 5.18e-06 came back as 1e-05 (+93%); true 2.14e-06 came back as 0. Relative branch-length error up to ~100%; topology unaffected.
Proposed change: multiply the matrix by a fixed scale factor (e.g. 1e4) before writing `infile` and divide `tree$edge.length` by it after reading — invariant for FM (weights are relative), restoring precision.

**F11 [MEDIUM, confidence: high] — a valid `outgroup` does not root the returned tree (DOC5; same class as gl.tree.nj F5)**
`R/gl.tree.fitch.r:135-143, 259-262` — the outgroup index is passed to PHYLIP's `O` option, which only reorients PHYLIP's own drawing; `ape::is.rooted()` on the return is `FALSE` (reproduced). The `@param` text "default NULL, no outgroup, tree not rooted" promises the converse for a supplied outgroup.
Proposed change: `tree_1 <- ape::root(tree_1, outgroup = <label>, resolve.root = TRUE)` before plotting/returning, or correct the documentation. **Consequence: the returned tree's representation changes for calls that pass `outgroup`.**

**F12 [LOW, confidence: high] — missing-binary errors have empty condition messages (VRB2; style-guide error idiom)**
`R/gl.tree.fitch.r:170-181, 208-219` — `cat(error(...))` followed by bare `stop()`.
Failure scenario (reproduced): `tryCatch`/logging captures `conditionMessage == ""`; in a pipeline the failure is undiagnosable from the condition object.
Proposed change: `stop(error("Fatal Error: cannot find", prog, "in", phylip.path, "\n"))`.

**F13 [LOW, confidence: high] — documentation/default mismatches (DOC1, DOC5, DOC3)**
`R/gl.tree.fitch.r:21` `bstrap "[default 1000]"` vs signature default `1`; `:31` `offset "[default 1.8]"` vs `1.2`; `:74` `@return "The tree file in newick format"` vs an ape `phylo` object; `:57` the `\dontrun` example calls `gl.phylip()`, a name that does not exist in the package; `:35` the `verbose` default clause deviates from the DOC2 canon.
Failure scenario: a user copying the example gets `could not find function "gl.phylip"`; a user relying on the documented default believes bootstrapping is on.
Proposed change: correct the five items; run `devtools::document()` (DOC4).

**F14 [LOW, confidence: high] — `tree.method` is a dead parameter with an inverted, ungated message (DOC5; same class as gl.tree.nj F3)**
`R/gl.tree.fitch.r:129-133` — any string is accepted; FM runs regardless; the "warning" (phrased as information, printed at all verbosities) fires only when the value is *not* FM. The `@details` do say only FM is implemented.
Proposed change: `match.arg(tree.method, c("FM"))` or remove the parameter (removal is API2).

**F15 [LOW, confidence: medium] — bootstrap replicates ignore `global.rearrange`/`randomize`; seeds are hard-coded (DOC5)**
`R/gl.tree.fitch.r:367-374` — the search options are commented out for replicate runs, and seeds 12345/331 are fixed and undocumented.
Failure scenario: the best tree and its replicates are searched with different intensity; users cannot vary the PHYLIP-side seed.
Proposed change: document the behaviour (or honour the options for replicates); note the R-side resampling is still governed by the session RNG.

**F16 [LOW, confidence: medium] — no consistency check between `D` and `x` (FS5)**
`R/gl.tree.fitch.r:109-116, 324-332` — nothing verifies that `x`'s populations match `attr(D, "Labels")` (set or order); the outgroup *index* computed from `D`'s labels (line 140) is applied to replicate matrices derived from `x`.
Failure scenario: `D` from a filtered subset, `x` the full object — replicate trees carry different taxa; consense mixes them or PHYLIP errors mid-loop with the cwd already moved (F5 compounds).
Proposed change: fail fast when `sort(popNames(x)) != sort(attr(D, "Labels"))`.

**F17 [LOW, confidence: medium] — 10-character name truncation can create duplicate taxon labels (FS5)**
`R/gl.tree.fitch.r:230-231, 330-331` — population names are truncated to 10 characters for the PHYLIP format with no collision check (testset names such as `EmmacBrisWiv` are already at 12).
Failure scenario: two populations sharing a 10-character prefix become indistinguishable tips; PHYLIP output and `read.tree` silently misassign.
Proposed change: check `anyDuplicated(substr(labels, 1, 10))` and stop with the offending names.

**F18 [INFO, confidence: high] — verbosity monotonicity oddities (VRB1)**
`R/gl.tree.fitch.r:417-423` — the bootstrap newick echo is gated `verbose == 3` exactly, so `verbose = 5` ("full report") prints *less* than 3; `:290-296` — `readLines` runs at `verbose >= 2` but its output prints only at `> 2` (dead read at exactly 2).
Proposed change: gate both blocks `verbose >= 3`.

**F19 [INFO, confidence: medium] — scripted menu answers depend on PHYLIP prompt sequencing (STY3)**
`R/gl.tree.fitch.r:361-379, 395-402` — the leading "A" in the bootstrap `fitch.cmd` and `consense.cmd` answers the "outfile exists — Replace/Append/File/Quit" prompt, which only appears because the first run's `outfile` is deliberately left in place. Verified working with 3.695; a different PHYLIP version, or any change to the deletion logic above, silently desynchronises every subsequent menu answer.
Proposed change: none required now; a comment documenting the protocol would prevent an innocent "cleanup" from breaking it.

**F20 [LOW, confidence: high, proposed rule] — `@author` names a custodian only, no `Author(s):` line (DOC7)**
`R/gl.tree.fitch.r:48-49`.
Proposed change: `Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}`.

### Notes on other functions (scope rule — one line each, not findings here)

- `gl.subsample.loc` / adegenet `SNPbin[]`: sampling with `replace = TRUE` corrupts NA genotypes on repeated indices (known, previously flagged, unfixed) — this poisons `gl.tree.fitch`'s bootstrap replicates whenever the data contain missing values; the baseline here dodged it only by using a callrate-1 fixture.
- `gl.dist.phylo`: not reviewed; the bootstrap path inherits whatever it does with its defaults (see F3).
- PHYLIP `fitch` 3.695 itself: claims "Tree also written onto file outtree" for 3 taxa while writing an empty file, and segfaults for 2 taxa (upstream, unfixable here; guarded by F9).

## Proposed changes

1. Replace the three `system("... < ...cmd")` invocations with `system2(prog, stdin = <cmd file>, stdout = <gated>)` on all platforms (F1, and the subprocess half of F7). **Consequence: the function changes from always failing on Windows to working; no output changes on platforms where it already ran.**
2. Rewrite bootstrap support handling: derive supports with `ape::prop.clades` on `tree_1` from the replicate trees, attach as `tree_1$node.label`, plot from `node.label` (F2, F4). **Consequence: numerical output changes — displayed support values change (they are currently wrong), and the returned object gains `node.label`.**
3. Forward the distance-model parameters to the replicate `gl.dist.phylo` calls (F3). **Consequence: bootstrap supports change whenever the user's `D` was not built with F81 defaults.**
4. Restore the working directory with `on.exit(setwd(hold), add = TRUE)` on every `setwd` (F5).
5. Honour `out.path` by copying `infile`/`outfile`/`outtree` there, or remove the parameter (F6; removal is API2).
6. Verbosity hygiene: gate the four ungated `cat(warn())` calls and the raw `cat` at line 356; suppress subprocess stdout below `verbose >= 2` (F7, F18).
7. Add `plot.display` to the signature with the house `verbose == 0` gate, and secure the return value before/independent of plotting (F8).
8. Guard `length(attr(D, "Labels")) >= 4` with an informative fatal error (F9).
9. Scale the distance matrix (×1e4) into `infile` and unscale `edge.length` after reading (F10). **Consequence: numerical output changes — branch lengths become more precise; topology unchanged.**
10. Root the returned tree on the supplied outgroup via `ape::root(..., resolve.root = TRUE)`, or amend the docs to state the tree stays unrooted (F11). **Consequence: returned tree representation changes for calls passing `outgroup`.**
11. Convert the two missing-binary exits to `stop(error(...))` (F12).
12. Documentation pass: `bstrap`/`offset` defaults, `@return`, example function name, DOC2 verbose text, hard-coded seeds and 4-taxon minimum documented, `Author(s):` line (F13, F15, F20); `devtools::document()` in the same change.
13. Validate `tree.method` with `match.arg` (F14).
14. Fail fast on `D`/`x` population mismatch and on 10-character name collisions (F16, F17).

## Coverage

- Standards walk (FS, DOC, VRB, DAT, DEP, PLT, STY): run.
- Spec axis, empirical: run — via the shimmed copy described above for everything past the `system()` calls; the unshimmed failure itself pinned separately.
- Independent verification (analysis mode): run — exhaustive 15-topology weighted-least-squares FM fit; RF 0; branch lengths agree within `outtree` precision (F10).
- Bootstrap path: run (`bstrap` 3 and 5, seed 42) — support extraction, return content, cwd checked.
- Edge cases: run — 3 taxa, 2 taxa, all-zero distances (returns an all-zero star tree), bad outgroup, unsupported `tree.method`, missing binary, genlight-as-`D`, even `n.jumble`, `bstrap > 1` without `x`.
- verbose = 0 silence: run (text side; plot side confirmed unconditional by code path and exercised under `pdf(NULL)`).
- SilicoDArT dispatch: PARTIAL — `D` is a `dist`, so datatype dispatch applies only to `x`; a SilicoDArT `x` is admitted for bootstrapping (default `accept`) and handed to `gl.dist.phylo`, which has a silico branch. Not exercised empirically: no silico distance fixture built. (DAT7 not cited as a finding since `gl.dist.phylo` owns the datatype-specific math.)
- Linux/macOS invocation path (`./fitch < fitch.cmd` via a shell): SKIPPED — Windows machine; redirection is expected to work there via `sh`, untested.
- FBM path (DAT6): SKIPPED — no FBM fixture; `x` is only touched via `gl.subsample.loc`/`gl.dist.phylo`.
- dartR Google Group / GitHub issues: searched 2026-09-06; no `gl.tree.fitch`-specific complaints found.
- `bstrap = 1000`-scale run: SKIPPED — runtime; protocol identical to the `bstrap = 5` run by construction (`M` answer is `bstrap`).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1–14 | approved | Arthur Georges (2026-09-06) | Via the formal approval boxes, each stated consequence acknowledged. F2+F3+F4 approved as a unit redesign of the bootstrap pipeline (all bootstrap outputs change; previous displayed values were fabricated). F1: function starts working on Windows, unchanged elsewhere. F11: returned tree representation changes for calls passing `outgroup`. F10: branch lengths change (more precise); topology unchanged. F19 (INFO): no action. |

## Outcome

All 14 approved changes applied on branch `review-gl.tree.fitch`
(base `ddaed27`, upstream/dev) and submitted as
[PR #371](https://github.com/green-striped-gecko/dartR.base/pull/371),
covering findings F1–F18 and F20; F19 is INFO/no action. Verification on this Windows machine with PHYLIP 3.695:

- (a) F1: unshimmed end-to-end run completes; RF = 0 to the Phase A pinned
  topology.
- (b) F2: `node.label` equals an independent recomputation of the clade
  frequencies from the replicate trees fitch wrote (canonicalised
  bipartition content, `bstrap = 5`, seed 42: supports 1.0, 0.4, 0.8) —
  exact match, tolerance 0.
- (c) F3: with `subst.model = "K80"`, the first replicate matrix written to
  `infile` equals a hand-computed same-seed
  `gl.dist.phylo(subst.model = "K80")` matrix to 3.8e-15, and differs from
  the old hard-wired F81 default by 1.1e-04.
- (d) F4: supports returned as `node.label`; PHYLIP majority-rule consensus
  attached as attribute `consensus.tree`; documented in `@return`.
- (e) F5: `getwd()` identical before/after a `bstrap = 5` run.
- (f) F6: `out.path` receives `infile`/`outfile`/`outtree` (and `intree`
  for bootstrap runs).
- (g) MEDIUMs: 3-taxa and 2-taxa inputs stop with "requires at least 4
  taxa" before any PHYLIP or plot work; a `verbose = 0` bootstrap run
  emits zero lines (PHYLIP dialogue suppressed) and renders no plot even
  with `plot.display = TRUE`; `outgroup` gives `is.rooted() == TRUE`; an
  invalid `plot.type` fails inside `tryCatch` and the tree is still
  returned; the distance matrix is scaled x10,000 into `infile` and branch
  lengths rescaled on read (documented), max |cophenetic - D| = 5.1e-05 on
  the fixture.
- (h) Updated characterization file: 28 assertions pass, 0 failures; every
  flipped expectation carries an `# [approved Fn]` tag; no unexplained
  diff.
- Caller grep across the 8 dartRverse clones: no callers of
  `gl.tree.fitch` outside its own file and tests — the signature additions
  (`subst.model`, `pairwise.missing`, `min.tag.len`, `plot.display`, all
  defaulted) break no sibling.

```json
{
  "function": "gl.tree.fitch",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "ddaed27",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "platform-portability (no catalogue rule; gap noted)/DOC5", "status": "applied", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "numerical-correctness/DOC5", "status": "applied", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "applied", "change": 2},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "applied", "change": 4},
    {"id": "F6", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 5},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "VRB5", "status": "applied", "change": 6},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "PLT3", "status": "applied", "change": 7},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "applied", "change": 8},
    {"id": "F10", "severity": "MEDIUM", "confidence": "medium", "rule": "numerical-correctness", "status": "applied", "change": 9},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 10},
    {"id": "F12", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "applied", "change": 11},
    {"id": "F13", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "applied", "change": 12},
    {"id": "F14", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "applied", "change": 13},
    {"id": "F15", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "applied", "change": 12},
    {"id": "F16", "severity": "LOW", "confidence": "medium", "rule": "FS5", "status": "applied", "change": 14},
    {"id": "F17", "severity": "LOW", "confidence": "medium", "rule": "FS5", "status": "applied", "change": 14},
    {"id": "F18", "severity": "INFO", "confidence": "high", "rule": "VRB1", "status": "applied", "change": 6},
    {"id": "F19", "severity": "INFO", "confidence": "medium", "rule": "STY3", "status": "no-action", "change": null},
    {"id": "F20", "severity": "LOW", "confidence": "high", "rule": "DOC7 (proposed rule)", "status": "applied", "change": 12}
  ],
  "coverage_skipped": [
    "Linux/macOS invocation path: Windows machine",
    "DAT6/FBM: no fixture",
    "SilicoDArT bootstrap x: no silico distance fixture",
    "bstrap=1000 scale run: runtime; protocol identical by construction"
  ],
  "status": "pr-open",
  "pr": 371
}
```
