# Review: gl.run.faststructure + gl.plot.faststructure (dartR.popgen)
- Family mode: analysis (fastStructure wrapper and its bar plot)
- Scope: one review, one report, one PR for the pair (agreed with Luis,
  2026-09-23); two manifest rows
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: a16da26 (origin/dev)
- Datasets: testset.gl, three populations (31 individuals, 17 loci);
  fastStructure (Structure_threader macOS binary, `~/programs/fastStructure`)
  and PLINK 1.9 (`~/programs/plink`); hand-built run objects for the plot
- Baseline: tests/testthat/test-gl.faststructure.R (plot tests run
  anywhere; run tests need `FASTSTRUCTURE_EXEC`, `PLINK_DIR` and a
  dartR.base whose `gl2plink` passes `--make-bed`)

## Verdict

**Standards: Rework** — `gl.run.faststructure` writes into the working
directory by default, checks neither binary, prints at every verbosity and
lacks the plot arguments; `gl.plot.faststructure` is a copy of
`gl.plot.structure` from before #94, with no `verbose`, checks or save
option.
**Spec: Rework** — `gl.run.faststructure` works only when `k.range` starts
at 2 and has no gaps, in an empty folder: `k.range = 3:4` and
`c(2, 4)` stop after all runs with "subscript out of bounds", and
`k.range = 1:2` drops K = 1 and returns an empty "3". A fixed `seed` makes
every replicate identical. `gl.plot.faststructure` carries the
`gl.plot.structure` defects fixed in #94 (modes shown as one replicate,
K = 1 duplicated, dendrogram on distances of distances).

What works well: for `k.range = 2:n` in a fresh folder, q-matrices are read
in `indNames(x)` order (the PLINK `.fam` order matches) with the right
marginal likelihoods, and the plot output feeds `gl.map.structure`.

**Environment note:** with the released dartR.base 1.2.3, `gl2plink` calls
PLINK without `--make-bed`, PLINK stops ("Basic file conversions do not
support regular filtering operations") and every fastStructure run fails.
dartR.base `dev` already passes `--make-bed`; all runs in this review
loaded dartR.base `dev` (commit on origin/dev, via a scratch worktree).
Users on CRAN dartR.base cannot run this function until dartR.base is
released.

## Findings

**F1 [HIGH, confidence: high] — results are indexed by K value in a list sized by the number of K (DOC5)**
`R/gl.run.faststructure.r:185-252` — likelihoods and q-matrices are stored
at position `K` of containers with `length(k.range)` (+1) rows, then the
first element is dropped.
Failure scenario: `k.range = 3:4` or `c(2, 4)`: all fastStructure runs
finish, then the call stops with "subscript out of bounds" and the runs
are lost. `k.range = 1:2`: the K = 1 result is dropped, `q_list[["2"]]`
holds K = 2 and `q_list[["3"]]` is NA.
Proposed change: key likelihoods and q-matrices by K (names = K values),
for any `k.range`, including K = 1 and gaps. **Consequence: calls with
`k.range` not starting at 2 now return results; `k.range = 1:n` now
includes K = 1 and has no empty element.**

**F2 [HIGH, confidence: high] — every earlier run in `output` is read back; default writes to the working directory (FS7, PLT2)**
`R/gl.run.faststructure.r:88, 179-183, 229-230` — all
`genotypes_output*` files in `output` are parsed, not only this call's,
and `output = getwd()`.
Failure scenario: running `k.range = 2` in a folder that holds an earlier
`k.range = 2:3, num.k.rep = 2` run stops with "subscript out of bounds";
other combinations can mix old and new runs silently. A default call
leaves 11+ files (PLINK and fastStructure) in the user's working
directory.
Proposed change: run in a new time-stamped folder created inside `output`,
with `output` defaulting to `tempdir()`, and read only that folder.
**Consequence: files are no longer written to the working directory by
default, and each call writes to its own subfolder.**

**F3 [MEDIUM, confidence: high] — a fixed seed gives identical replicates (DOC5)**
`R/gl.run.faststructure.r:151-174` — every replicate gets `--seed=seed`.
Failure scenario: `seed = 5, num.k.rep = 2` returns two identical
q-matrices per K, so the replicates add nothing and
`gl.plot.faststructure` sees one mode by construction.
Proposed change: replicate r uses `seed + r - 1`: reproducible and
different. **Consequence: results with a `seed` change for replicates
after the first.**

**F4 [MEDIUM, confidence: high] — no checks on binaries, folder or run success (FS5, VRB2)**
`R/gl.run.faststructure.r:106-123, 129-175, 194-199` — the executable,
PLINK and `output` are not checked; a failed fastStructure run is noticed
only when its log is read; Windows is not refused; the `gsubfn` guard
prints and returns -1 instead of stopping.
Failure scenario: a wrong `exec` or a missing `output` folder prints
fastStructure tracebacks for every run and then stops with "cannot open
the connection" (with a missing folder, PLINK writes its files to
`tempdir()` and fastStructure looks in `output`).
Proposed change: before running, stop with `error()` on Windows, when
`exec` or the PLINK binary is missing (with the download URLs), and create
`output` if needed; after each run, stop naming K and replicate if its
`.meanQ` file is missing; read the likelihood with base R and drop
`gsubfn`.

**F5 [LOW, confidence: high] — output ignores `verbose` (VRB1, VRB3)**
`R/gl.run.faststructure.r:128, 130-175` — "Running K = ... Replicate = ..."
is printed with `print()` at every verbosity, and PLINK and fastStructure
output always reaches the console.
Failure scenario: `verbose = 0` prints one line per run plus all program
output.
Proposed change: one progress line per run at `verbose >= 2`
(`cat(report())`), program output only at `verbose >= 3`, as in
`gl.run.structure`.

**F6 [LOW, confidence: high] — likelihood plot outside the dartR plot idiom (PLT1, PLT2)**
`R/gl.run.faststructure.r:208-221` — `theme_bw()`, x breaks fixed at 1:10
(K above 10 unlabelled), always printed, no `plot.out`, `plot.file` or
`plot.dir`; the data are passed as loose vectors, so the ggplot has no
data frame.
Failure scenario: `k.range = 2:12` has no labels for K = 11, 12; a script
cannot suppress or save the plot.
Proposed change: plot a data frame of mean marginal likelihood per K with
`plot.theme = theme_dartR()`, breaks at the K run, `plot.out`, `plot.dir`,
`plot.file` (`utils.plot.save()`).
**Consequence: the plot looks different (dartR theme).**

**F7 [HIGH, confidence: high] — gl.plot.faststructure is the pre-#94 gl.plot.structure (DOC5, PLT1, FS2, FS5)**
`R/gl.plot.faststructure.r:106-352` — the body repeats `gl.plot.structure`
as it was before #94, so it keeps the defects fixed there.
Failure scenario (from the baseline): six K = 2 replicates in two modes of
three return two tables that are each a single replicate, not the mode
means; `k.range = c(2, 1)` returns panels `2.1`, `2.2`, `1.1` (K = 2
twice); with `den = TRUE` the dendrogram clusters `dist()` of a distance
matrix; `aes_()` is deprecated; `k.range` has no default although
documented as NULL; there is no `verbose`, `plot.out`, `plot.file` or
argument check.
Proposed change: convert the fastStructure run object to the
`structure.result` layout (id, pct.miss = 0, orig.pop, Group.1..K per
replicate) and call `gl.plot.structure`, keeping the current arguments
(`colors_clusters` passed as `color_clusters`; `k.range` default NULL for
all K) and adding `dis.mat`, `plot.out`, `plot.dir`, `plot.file`,
`verbose`. `label.size` is kept by adding an optional `label.size = 12`
argument to `gl.plot.structure` (its current fixed value).
**Consequence: returned q-matrices change when a K has more than one mode;
`k.range` with K = 1 after another K no longer duplicates panels; the
returned tables are data.tables (as from `gl.plot.structure`); the
dendrogram order changes with `den = TRUE`.**

**F8 [LOW, confidence: high] — documentation (DOC1, DOC2, DOC5, DOC7, FS3)**
Both files — `gl.run.faststructure`: `exec.plink` described as "working
directory" but is a folder path; `@return` does not describe `q_list`
(K → replicate → data frame of id, orig.pop, V1..VK) or `plot`; Windows
limitation only in details; `utils.flag.start(build = "Jody")`.
`gl.plot.faststructure`: `k.range` default, `den` refers to
`gl.run.structure`, `@return` says "List of Q-matrices" without columns;
both lack `@family` and the Author(s)/Custodian structure (DOC7, proposed
rule).
Failure scenario: a user cannot tell from the help how to reach the
replicate q-matrices or that `exec.plink` is a folder.
Proposed change: rewrite both headers; regenerate Rd.

## Proposed changes

1. Results keyed by K for any `k.range` (F1). **Consequence: `k.range`
   not starting at 2 works; K = 1 is kept.**
2. Each call runs in a new subfolder of `output`, default `tempdir()`;
   only that folder is read (F2). **Consequence: no files in the working
   directory by default; old runs are never mixed in.**
3. `seed + replicate - 1` per replicate (F3). **Consequence: seeded
   replicates now differ.**
4. Binary, OS and folder checks; stop on a failed run; drop `gsubfn` (F4).
5. `verbose` gating of progress and program output (F5).
6. Likelihood plot with `plot.theme`, `plot.out`, `plot.dir`, `plot.file`
   (F6). **Consequence: plot look changes.**
7. `gl.plot.faststructure` delegates to `gl.plot.structure`; optional
   `label.size` added to `gl.plot.structure` (F7). **Consequence: output
   changes as in #94 for modes, K = 1 and dendrograms; returned tables are
   data.tables.**
8. Documentation for both, Rd regenerated, NEWS entry (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both.
- Spec: real fastStructure + PLINK runs (k.range 2:3, 3:4, c(2, 4), 1:2,
  reused folder, seed, default output, wrong exec) — run; plot on
  hand-built run objects (one mode, two modes, K = 1) — run.
- `cv > 0` and `prior = "logistic"`: SKIPPED — not exercised; parsing of
  their logs unchecked.
- Windows: SKIPPED — no Windows machine (fastStructure has no Windows
  binary).
- Released dartR.base 1.2.3: run fails in `gl2plink` (see environment
  note); all other runs used dartR.base `dev`.
- dartR Google Group / GitHub issues search: not run (no search access).
- Downstream callers: none in other `dartR.*` packages. dartr2shiny runs
  both (`input_Tabs.csv`, `variables_matrix.csv`: the run result is passed
  to the plot as `Myfaststructure`); every proposed change keeps the
  returned `list(q_list, plot)` shape and adds arguments only after the
  existing ones.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |
| 8 | approved | Luis | |

## Outcome

- Changes 1–8 applied in commit 13ae6ac on `review-faststructure`, PR #101 to `dev`.
- Characterization test: 10 tests, 37 expectations pass with both binaries;
  every diff from baseline is tagged `[approved n]`.
- Old vs new, seeded single replicate K = 2:3: identical `q_list` and
  marginal likelihoods; plot tables for single-mode runs identical.
- `verbose = 0` silent in the terminal; needed suppressing `gl2plink`'s
  `message()` of PLINK output (a dartR.base issue: it ignores `verbose`).
- A first rewrite dropped `@importFrom stats ...` and `@import ggdendro`
  from the plot header, which other functions rely on; restored before
  commit, `NAMESPACE` unchanged.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; no new NOTE.

## Machine block

```json
{
  "function": "gl.run.faststructure",
  "also_covers": ["gl.plot.faststructure"],
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "a16da26",
  "verdict_standards": "rework",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS7", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "PLT1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 8}
  ],
  "coverage_skipped": [
    "cv > 0 and logistic prior: not exercised",
    "Windows: no machine",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 101
}
```
