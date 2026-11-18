# Review: gl.blast (dartR.popgen)

- Family mode: analysis (external-tool wrapper that returns an annotated
  genlight, so the modify-family checks DAT2 and FS8 were applied as well)
- Date: 2026-09-18
- Reviewer: Claude (claude-fable-5-1), dartr-function-review v2.0.0
- Package commit: 31fb40c (dev_luis, synced with origin/dev)
- Datasets: testset.gl (dartR.data, 255 loci, TrimmedSequence 20 to 69 nt)
  against a synthetic three-contig reference built by the test from the
  object's own sequence tags (40 tags embedded exactly, 30 with two
  substitutions, one random contig); the same tags written to a query fasta
  for the file-input path. BLAST+ 2.16.0 (blastn/makeblastdb on PATH),
  macOS.
- Baseline: tests/testthat/test-gl.blast.R (dartR.popgen), 7 tests, all
  passing on the reviewed state. Three of them pin defects (F1, F3, F4).
- Manifest note: the campaign manifest lives in dartR.base. That checkout's
  dev_luis could not be merged with origin/dev (uncommitted local changes
  in man/ block the merge), so the claim row was added to the local file
  only; the dartR.popgen checkout is synced.

**Standards: Needs work** - the house structure is followed, but the two
shell calls ignore their exit status, paths are never quoted, and console
output does not honour verbose.
**Spec: Needs work** - the core promise holds (one best hit per tag, merged
into loc.metrics with locus order and row count intact, thresholds
applied as documented), but every failure of the external tools returns the
previous run's result as if it were this run's.

What works well: the filter and tie-break match the documentation exactly,
the merge keeps 255 rows in the original locus order, ploidy is untouched,
and the gzipped-reference path gives identical hits to the plain fasta.

## Findings

**F1 [BLOCKER, confidence: high] - a failed run returns the previous run's
hits (FS5 fail fast; DOC5)**
`R/gl.blast.r:182-187`, `:236-244`, `:284-305`, `:307` - the return values
of `file.copy()`, `system(makeblastdb ...)` and `system(blastn ...)` are
discarded. Every intermediate file (`fasta.input`, `db_blast.*`,
`output_blast.txt`) has a fixed name in `tempdir()` and is never removed, so
whatever the failed step did not overwrite is read as this run's result.
`task` is not validated, and neither `ref_genome` nor a path `x` is checked
for existence.
Failure scenario (all reproduced, in the baseline test): after one
successful call in the session, (a) `task = "blastx"` makes blastn exit
with an error and the function returns the 59 hits of the previous call
with "Completed: gl.blast"; (b) a query fasta path that does not exist is
silently not copied and the previous query (the genlight's numbered tags)
is BLASTed and returned as a data frame; (c) a reference-genome path
makeblastdb cannot read (see F2) leaves the previous database in place and
blastn searches it. On a fresh session the same mistakes stop with
"missing value where TRUE/FALSE needed" from `file.info()$size`.
Proposed change: validate up front (`x` is a genlight or an existing file,
`ref_genome` exists, `task` via `match.arg`, thresholds numeric), unlink
the three intermediate files before each run, and stop with the tool's own
message when makeblastdb or blastn exits non-zero.

**F2 [HIGH, confidence: medium] - paths are passed to the shell unquoted
(FS5; platform-specific calls lens)**
`R/gl.blast.r:236-244`, `:284-305` - `ref_genome`, `tempdir()` and the
executable paths are pasted into the command line bare. Lines 203-210 and
260-267 work around this on Windows by refusing any executable path that
contains a space, which is where the NCBI Windows installer puts BLAST by
default (under Program Files). A `tempdir()` under a Windows user name with
a space breaks `-out`/`-query` on every call.
Failure scenario (reproduced on macOS): `ref_genome = "dir with
space/genome.fasta"` makes makeblastdb exit with "Too many positional
arguments (1), the offending value: with"; combined with F1 the call then
returns the previous run's 59 hits. Windows behaviour is inferred from the
code, not run (see Coverage).
Proposed change: build both commands with `system2()` and a quoted
argument vector (`shQuote`, type `cmd` on Windows), and drop the
spaces-in-path refusal, which then has no purpose. This also removes the
newline embedded in the `-outfmt` string at lines 302-303, which works
through `sh` on macOS but is untested on Windows.

**F3 [MEDIUM, confidence: high] - a second run on an annotated object
suffixes the BLAST columns .x/.y (DOC5; DAT2 principle that loc.metrics
describe the loci they sit against)**
`R/gl.blast.r:383-393` - the merge adds the 20 BLAST columns without
removing those left by an earlier run.
Failure scenario (reproduced): `gl.blast(gl.blast(x, ref), ref)` returns
loc.metrics with 60 columns, `sacc.x`/`sacc.y` and no `sacc`, so
`gl.filter.locmetric(metric = "evalue")` on the result fails to find the
metric. Rerunning with a different task or threshold is the normal way to
use the function. The 2023 Google Group thread "blast filtering" reports
users who could not find the BLAST metrics in loc.metrics after running
the function; this and F4c are the two code paths that produce that
symptom.
Proposed change: drop any existing BLAST columns from loc.metrics before
merging, so the latest run replaces the previous one.
**Consequence: a second gl.blast on the same object overwrites the earlier
run's BLAST columns instead of keeping both sets under .x/.y.**

**F4 [MEDIUM, confidence: high] - the no-hit and no-survivor exits do not
match the documentation or the verbosity contract (DOC5, VRB3, FS8, FS9,
VRB4 proposed rule)**
`R/gl.blast.r:307-325`, `:381`, `:395-399` - three related gaps.
(a) With a fasta path as `x` and no alignment, the function returns the
input path string; `@return` promises a data frame.
(b) "No sequences were aligned" prints at `verbose = 0`, and that exit
skips both the history entry and the "Completed" line that every other exit
prints.
(c) When BLAST finds hits but none pass `Percentage_overlap`/`bitscore`,
`plyr::rbind.fill()` of an empty list is NULL, the message reads
"  sequences were aligned after filtering" with a blank count, no column is
added, and the run ends with "Completed". Nothing tells the user the
thresholds removed everything.
Failure scenario (all reproduced): `bitscore = 1e9` on testset.gl returns
the object unchanged with the blank-count message; the forum symptom in F3
is the same outcome seen by a user.
Proposed change: one exit path that reports "k of N sequences aligned after
filtering" at `verbose >= 1` (results-affecting, VRB4), keeps history and
"Completed" on every path, gates the no-alignment message, and for fasta
input with no surviving hit returns a zero-row data frame with the 20 BLAST
columns.
**Consequence: for a fasta-file query with no hit, the return value changes
from the input path string to an empty data frame.**

**F5 [MEDIUM, confidence: high] - external tool output ignores verbose
(VRB1)**
`R/gl.blast.r:236-244`, `:279-281`, `:327-329` - makeblastdb's build log
("Building a new DB ...", eight lines) reaches the console on every call
including `verbose = 0`, because `system()` inherits stdout. "Starting
BLASTing" and "Starting filtering" are progress messages printed at
`verbose >= 1`, a level reserved for begin and end.
Failure scenario (reproduced): `capture.output(gl.blast(x, ref, verbose =
0))` captures nothing while eight lines of makeblastdb output appear on the
terminal; a Shiny or batch log fills with tool chatter.
Proposed change: pass `stdout = FALSE` to the tool calls below
`verbose >= 3` (errors surface through F1's status check), and move the
two progress lines to `verbose >= 2`.

**F6 [LOW, confidence: high] - the NOTE points to functions that do not
exist (DOC5)**
`R/gl.blast.r:418-425` - "Retrieve output files from tempdir using
gl.list.reports() and gl.print.reports()" names functions that are exported
by no dartRverse package (grep of all six NAMESPACE files). The three RDS
tables are written under random `tempfile()` names, so the unfiltered and
filtered tables are unreachable in practice.
Failure scenario: a user following the note gets "could not find function".
Proposed change: print the three file paths at `verbose >= 2` instead of
the note.

**F7 [LOW, confidence: high] - documentation drift (DOC2; DOC5; DOC6 and
DOC7, proposed rules)**
`R/gl.blast.r:43-45` - the verbose text is not the standard DOC2 wording.
`:93-95` - "saved to the working directory (plot.dir tempdir if not set)":
there is no `plot.dir` argument and nothing is written outside `tempdir()`.
`:105-107` - `@return` does not cover the no-hit case (F4a).
`:109-110` - `@author` names authors but no Custodian label (DOC7).
`:132` - `@seealso gl.print.history` is unrelated; `gl.filter.locmetric`
is the function users need next (the forum thread asks exactly this).
`:24-28`, `:33`, `:57-72`, `:124`, `:129` - curly quotes, an en dash and an
accented character in roxygen (DOC6). `man/gl.blast.Rd` is in sync with the
header (DOC4 checked).
Proposed change: docs-only edits listed above, then `devtools::document()`.

**F8 [INFO, confidence: high] - style and idiom (FS3, FS11 proposed, STY3)**
`R/gl.blast.r:157`, `:168`, `:383`, `:428`, `:438` - `class(x)[1] ==
"genlight" | class(x)[1] == "dartR"` five times where `is(x, "genlight")`
covers both (FS11). An input that is neither (a genind, reproduced) falls
through to `file.copy()` and stops with "invalid 'file' argument".
`:152-154` - `build = "Jody"` is the outdated `utils.flag.start` argument
(FS3). `:378` - `decreasing = T`. `:316` - `read.table(quote = "\"")`
strips double quotes inside subject titles (reproduced: `plasmid "quoted"
name` comes back as `plasmid quoted name`); `quote = ""` keeps titles
verbatim.
Proposed change: the input-type check joins change 1; the rest is a
tidy-up.

## Proposed changes

1. Fail fast and never reuse a previous run's files: validate `x`,
   `ref_genome` and `task`, unlink the three intermediate files before each
   run, stop with the tool's message on a non-zero exit (F1, F8 input
   check).
2. Build both commands with `system2()` and quoted arguments; drop the
   spaces-in-path refusal and the embedded newline in `-outfmt` (F2).
   Consequence: a BLAST install under a path with spaces (the Windows
   default) now runs instead of stopping. No numerical change.
3. Replace, not suffix, the BLAST columns when the object already carries
   them (F3). **Consequence: a second gl.blast on the same object
   overwrites the earlier run's BLAST columns.**
4. One exit path with the "k of N" count at `verbose >= 1`, history and
   "Completed" on every path, gated no-alignment message, progress lines at
   `verbose >= 2`, tool stdout hidden below `verbose >= 3` (F4b, F4c, F5).
5. Fasta-file query with no surviving hit returns a zero-row data frame
   with the BLAST columns (F4a). **Consequence: return type changes from
   the input path string to an empty data frame in that case.**
6. Print the paths of the three saved tables instead of the NOTE naming
   `gl.list.reports()`/`gl.print.reports()` (F6).
7. Documentation fixes and `devtools::document()` (F7).
8. Style tidy-up: `is(x, "genlight")`, `TRUE`, drop `build =`,
   `quote = ""` in `read.table` (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT (not applicable, no plot),
  STY, API - run
- Spec: behaviour vs roxygen on testset.gl against the synthetic reference,
  15 probes (genlight, fasta path, gz, no hit, nothing surviving the filter,
  rerun, bad task, missing paths, path with a space, genind input, verbose
  0 capture, dc-megablast on the mutated contig) - run
- DOC4 Rd/NAMESPACE sync: roxygenise in a scratch copy, no diff - run
- Callers: dartr2shiny `Fun_gl.blast.R` calls with all seven named
  arguments; no sibling dartR.* package calls gl.blast - run
- Google Group: thread "blast filtering" (groups.google.com/g/dartr/c/
  q2WCRJCYyqY, Aug to Sep 2023) - read; symptoms consistent with F3/F4c
- Windows path (`where`, cmd quoting, `system()` with an embedded newline):
  SKIPPED - no Windows machine; F2 confidence set to medium for that reason
- FBM path (DAT6): not exercised - the function reads `nLoc()` and
  loc.metrics only and never touches genotypes
- Large reference genome runtime and memory: SKIPPED - out of scope for a
  correctness review

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | rejected | Luis | keep both runs under .x/.y; overwrite not wanted |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |
| 8 | approved | Luis |  |
| F9 addendum | approved | Luis | tempdir() with a space refused with a TMPDIR message |

## Outcome

- Change 1 (F1): inputs validated before any work (`x` genlight or existing
  file, `ref_genome` exists, `task` via `match.arg`, numeric ranges);
  `fasta.input`, `db_blast.*` and `output_blast.txt` unlinked before each
  run; `system2()` exit status checked and a non-zero exit stops with the
  tool's name and status. Evidence: test "invalid inputs and a failed tool
  stop instead of returning the previous run's hits", six error
  expectations including makeblastdb exit 1 on an empty genome file.
- Change 2 (F2): both commands built with `system2()` and `shQuote()`;
  spaces-in-path refusal and the embedded newline removed. Applying it
  showed that makeblastdb splits the values of `-in` and `-out` on
  whitespace and blastn does the same for `-db` (verified in the shell:
  `-in "dir with space/genome.fasta"` fails, stdin input works, blastn
  `-query`/`-out` with spaces work, `-db` with spaces fails), so the genome
  is now fed on stdin with `-title ref_genome`. Evidence: test "paths with
  spaces are quoted for the shell", 59 hits from a genome and a query fasta
  under `dir with space/`.
- Addendum F9 (found while applying change 2, no approval yet): because
  BLAST cannot open a database under a path with spaces, a `tempdir()`
  containing a space is refused up front with a message naming
  `TMPDIR`/`TMP`. Before, such a session failed inside blastn with
  "Database memory map file error" (and, with F1, returned stale hits).
  Not a change in results. Asked for approval in the pre-push message.
- Change 4 (F4b, F4c, F5): one exit path; "k of N sequences aligned after
  filtering" at `verbose >= 1` (as `warn()` when k = 0); history and
  "Completed" on every exit; progress lines at `verbose >= 2`; tool stdout
  hidden below `verbose >= 3`. Evidence: `expect_silent()` at `verbose = 0`
  in three tests; `expect_output()` of "59 of 255" and "0 of 255".
- Change 5 (F4a): fasta query with no surviving hit returns a zero-row data
  frame with the 20 BLAST columns. Evidence: test "no alignment".
- Change 6 (F6): the three RDS paths are printed at `verbose >= 2`.
  Evidence: verbose 3 run in the Outcome log below.
- Change 7 (F7): docs edited as listed; `devtools::document()` regenerated
  `man/gl.blast.Rd` (62 lines changed). The accented author name in the
  1997 reference was kept as is (Encoding: UTF-8 is declared).
- Change 8 (F8): `is(x, "genlight")`, `TRUE`, `build =` dropped,
  `quote = ""` in `read.table`, `pmin()` for the overlap denominator.
- Change 3 rejected: the merge is untouched; the rerun test still pins the
  `.x`/`.y` suffix behaviour.
- Snapshot: 5 of the 29 baseline expectations changed, each mapped to an
  approved change: ungated message x2 and missing history entry (change 4),
  path string returned (change 5), stale hits after a bad `task` (change
  1). The 59 hits, 41 columns, 36 + 23 contig counts, locus order and
  ploidy are unchanged. Final test file: 8 tests, 47 expectations, all
  passing.
- End-to-end at `verbose = 3` on testset.gl against the synthetic genome:
  59 of 255 aligned, database build log shown, three RDS paths printed,
  "Completed: gl.blast".
- `devtools::check()`: 0 errors; a `setNames` NOTE raised against gl.blast
  was fixed before committing. The remaining 2 warnings and 3 notes come
  from untracked local files in the checkout and ade4/ggplot2/dplyr built
  under a newer R.
- Addendum F9 approved by Luis with the pre-push OK.
- PR: green-striped-gecko/dartR.popgen#90, commit 8eea717 on dev_luis.

```json
{
  "function": "gl.blast",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "31fb40c",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "medium", "rule": "FS5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "rejected", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": [4, 5]},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "FS11", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["Windows path: no Windows machine", "DAT6: not applicable", "large-genome runtime: out of scope"],
  "status": "pr-open",
  "pr": 90
}
```
