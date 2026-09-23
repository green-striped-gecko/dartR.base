# Review: gl.read.structure (dartR.popgen)
- Family mode: io (reads STRUCTURE output files into a `structure.result`)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 063c0c7 (origin/dev)
- Datasets: STRUCTURE 2.3.4 output written by `gl.run.structure(delete.files
  = FALSE)` for three populations of testset.gl (31 individuals, 18 loci):
  K = 1:3 without priors, K = 2:3 with `pop.prior = "usepopinfo"`; a real
  hand-run folder (80 files, K = 1:8 x 10 replicates, 251 individuals,
  `~/dartR.popgen/str/`, not committed)
- Baseline: tests/testthat/test-gl.read.structure.R with fixture files in
  tests/testthat/fixtures/structure/ (72 KB; runs without the STRUCTURE
  binary)

## Verdict

**Standards: Needs work** — the function has most of the house anatomy but
uses the obsolete `build =` flag argument, guards packages that are
Imports, prints messages with doubled spaces, and has no examples, author
or details.
**Spec: Rework** — for runs with `usepopinfo`, every individual with a
population prior is read with wrong membership (0.5/0.5 instead of
0.042/0.958 in the fixture) and an all-NA `prior.anc`; with `x`, files
whose ids are not full individual names get `orig.pop = NA` for every
individual without a message. The parser was rewritten from the one in
`utils.structure.run`, which reads the same files correctly.

What works well: runs without priors are read exactly — every q-matrix and
summary matches `gl.run.structure`'s own parse of the same files, and the
real 80-file folder (251 individuals) reads in 7 s with correct LnP values.

## Findings

**F1 [HIGH, confidence: high] — usepopinfo individuals read with wrong q and empty prior ancestry (DOC5)**
`R/gl.read.structure.r:222-262` — for a line such as
`2  10  (5)  2 :  0.958 | Pop 1: 0.002 0.012 0.028 |`, the first token
after `|` is the word "Pop", so `as.integer(b[1])` is NA and every ancestry
row is skipped. `rowSums(anc, na.rm = TRUE)` is then 0 for all
populations, the fallback copies `Group.1` (0.958) into every group, and
the row is normalised to equal shares. `df$Group.1 <- NULL` also discards
the probability of the individual's own population.
Failure scenario: `usepopinfo` fixture, K = 2, individual 10 (population
2, POPFLAG = 1): the correct q is 0.042 / 0.958; `gl.read.structure`
returns 0.5 / 0.5, and `prior.anc` is NA for every flagged individual. At
K = 3 the largest error is 0.66. `gl.run.structure` parses the same file
correctly.
Proposed change: parse the ancestry blocks as `utils.structure.run` does
(drop the "Pop" token; the own population's row takes the value before
`|`, other populations the sum over generations), without `na.rm`.

**F2 [MEDIUM, confidence: high] — `x` matched by exact id only; unmatched individuals get NA silently (DAT5)**
`R/gl.read.structure.r:344-362` — `orig.pop` is replaced by a `merge()` on
`id` with `indNames(x)`; nothing checks how many matched.
Failure scenario: files written by `gl.run.structure` label individuals
by index (1..n), so `gl.read.structure(folder, x = x)` returns
`orig.pop = NA` for all 31 individuals, overwriting the population numbers
the file had. `gl.map.structure` then fails with "attempt to select less
than one element". STRUCTURE also truncates labels longer than 11
characters, so hand-run files lose the match for those individuals.
Proposed change: match by name; when the ids are exactly `1..nInd(x)`
(the `gl.run.structure` layout) match by index and restore the names in
`id`, as `gl.run.structure` does; stop with an error naming the unmatched
ids when some ids match neither way.

**F3 [MEDIUM, confidence: high] — a STRUCTURE run folder cannot be read without `pattern` (DOC5)**
`R/gl.read.structure.r:53-63, 330-339` — every file in the folder is
parsed, and the first file that is not STRUCTURE output stops the call.
Failure scenario: `gl.read.structure(folder)` on a folder kept by
`gl.run.structure(delete.files = FALSE)` stops with "Could not locate
Q-matrix block in file: k1.r1_log". The same happens with any README,
data or params file next to the outputs.
Proposed change: keep only files that contain "Estimated Ln Prob of Data"
(STRUCTURE output), skip the rest with a note at `verbose >= 2`, and stop
only when no output file is left.

**F4 [LOW, confidence: high] — replicate numbers and labels do not follow the files (DOC5)**
`R/gl.read.structure.r:303-315` — files are sorted alphabetically within K,
so `rep10` sorts before `rep2`, and labels are built from the longest
common file-name prefix.
Failure scenario: the real folder `str_K1_rep1_f ... str_K1_rep10_f` gives
`str_K1_rep10_f` the label `str_K.k1.r2` and `str_K1_rep2_f` the label
`str_K.k1.r3`; the fixture folder gives `k.k1.r1`. A user matching a label
back to a file picks the wrong run. `gl.run.structure` names runs
`k1.r1`.
Proposed change: order replicates by the numbers in the file name (natural
order) and name runs `k<K>.r<replicate>` as `gl.run.structure` does, with
`prefix` prepended only when the user sets it. **Consequence: run names
change (e.g. `k.k2.r1` → `k2.r1`).**

**F5 [LOW, confidence: medium] — `rename_files = TRUE` can overwrite existing files (principle: do not destroy user data)**
`R/gl.read.structure.r:318-324` — `file.rename()` replaces an existing file
of the target name without asking.
Failure scenario: a folder that already holds `<label>_out` (for example
from an earlier renaming) loses that file when another output is renamed
onto the same name.
Proposed change: stop before renaming when any target file exists, naming
it.

**F6 [LOW, confidence: high] — house structure and messages (FS3, FS5, VRB2, DEP1)**
`R/gl.read.structure.r:39-47, 319, 328, 366` — `utils.flag.start()` gets
the obsolete `build = "Jody"`; `purrr` and `dplyr` are guarded although
they are Imports; messages built as `report("Processing ", n, " STRUCTURE
...")` print doubled spaces; the end flag starts with a blank line; `x` is
checked only after all files are parsed.
Failure scenario: `verbose = 2` prints "Processing  6  STRUCTURE output
files."; a wrong `x` is reported after the whole folder has been read.
Proposed change: drop `build =` and the Imports guards; single-spaced
messages; check `x` before reading.

**F7 [LOW, confidence: high] — documentation gaps (DOC1, DOC3, DOC7)**
`R/gl.read.structure.r:1-26` — no `@examples`, `@author`, `@family` or
`@details`; the help does not say which files are read, how K is
detected, that `orig.pop` is a population number unless `x` is given, or
how `x` is matched; `verbose` text is not the DOC2 wording.
Failure scenario: a user passes the output of `gl.read.structure(folder)`
(no `x`) to `gl.map.structure` and gets "No coordinates in x for these
populations of qmat: 1, 2, 3" without knowing why.
Proposed change: complete the header, with a `\dontrun{}` example reading
a folder kept by `gl.run.structure`; regenerate Rd.

**F8 [INFO, confidence: high] — two STRUCTURE output parsers (STY3)**
`R/gl.read.structure.r:102-288` and `R/utils.structure.run.r:112-447` — the
same file format is parsed by two separate implementations; F1 is a
divergence between them.
Failure scenario: a future fix to one parser does not reach the other.
Proposed change: none now; a shared internal reader used by both is a
refactor across two functions, better done as its own change.

**Out of scope, noted for gl.map.structure:** a q-matrix whose `orig.pop`
is all NA makes `gl.map.structure` stop with "attempt to select less than
one element in integerOneIndex" instead of a clear message (seen through
F2). Not changed here.

## Proposed changes

1. Parse usepopinfo ancestry lines correctly (F1). **Consequence: q-matrix
   values and `prior.anc` change for every individual with a population
   prior.**
2. Match `x` by name, or by index when ids are `1..nInd(x)` (restoring the
   names in `id`); error naming unmatched ids (F2). **Consequence: with
   `x`, files written by `gl.run.structure` now carry population names and
   individual names instead of NA and index numbers; files with
   unmatched ids now stop with an error instead of returning NA.**
3. Read only STRUCTURE output files from the folder; skip others with a
   note (F3). **Consequence: folder reads that errored now succeed.**
4. Natural replicate order and `k<K>.r<replicate>` run names, `prefix`
   only when set (F4). **Consequence: run names and the replicate number
   given to each file change.**
5. Refuse to rename onto an existing file (F5).
6. Drop `build =` and Imports guards, single-spaced messages, check `x`
   before reading (F6).
7. Complete the documentation with an example; regenerate Rd; NEWS entry
   (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run.
- Spec: every run read back compared with `gl.run.structure`'s parse of the
  same files (independent parser) — run, no-prior and usepopinfo; real
  hand-run folder (80 files) — run for timing, K detection, LnP and row
  sums.
- Pipeline: result passed to `gl.plot.structure`, `gl.evanno`,
  `gl.map.structure` — run.
- POPDATA = 0 files (no population column) and LOCPRIOR outputs: SKIPPED —
  no fixture; `gl.run.structure` always writes POPDATA = 1.
- Files from Windows (CRLF line endings): SKIPPED — no fixture.
- `rename_files = TRUE`: exercised by reading the code only; not run on
  user data.
- dartR Google Group / GitHub issues search: not run (no search access).
- Downstream callers: none in `dartR.*` packages or dartr2shiny (grepped).

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

## Outcome

- Changes 1–7 applied in commit a495713 on `review-gl.read.structure`, PR #98 to `dev`.
- Characterization test: 42 expectations pass; every diff from baseline is
  tagged `[approved n]` (1–6).
- Independent check: with `x`, fixture folders read back equal
  `gl.run.structure`'s parse of the same files (q max diff 0; ids,
  `orig.pop`, `prior.anc`, summaries identical), no prior and usepopinfo.
- Real 80-file folder: q-matrices and summaries identical to the old
  reader; 6.8 s vs 6.5 s.
- Output-file detection needed both "Estimated Ln Prob of Data" and
  "Estimated Allele Frequencies": the STRUCTURE run log repeats the
  likelihood line (found while applying change 3; within its scope).
- `devtools::document()` also added the missing `@family` cross-links to
  `gl.evanno.Rd`, `gl.map.structure.Rd`, `gl.plot.structure.Rd`.
- `R CMD check`: unrelated failure in `test-gl.ld.haplotype.R`; no new NOTE.

## Machine block

```json
{
  "function": "gl.read.structure",
  "package": "dartR.popgen",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "063c0c7",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "medium", "rule": "principle: do not destroy user data", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "STY3", "status": "no-change", "change": null}
  ],
  "coverage_skipped": [
    "POPDATA = 0 and LOCPRIOR outputs: no fixture",
    "CRLF files: no fixture",
    "rename_files on user data: code reading only",
    "Google Group / issues search: no access"
  ],
  "status": "pr-open",
  "pr": 98
}
```
