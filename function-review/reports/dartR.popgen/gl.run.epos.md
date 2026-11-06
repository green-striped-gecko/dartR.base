# Review: gl.run.epos (dartR.popgen)
- Family mode: analysis (wrapper around EPOS, epos2plot and bootSfs)
- Custodian: Bernd Gruber (STY5: this report is the discussion record;
  changes approved by Luis)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 1e90a1e (origin/dev)
- Datasets: possums.gl[1:30, 1:200] after gl.filter.monomorphs (30
  individuals, 176 loci), and without the filter (24 monomorphic loci) for
  the zero-class check; testset.gs[1:20, 1:60]
- External software: EPOS, epos2plot (`~/programs`) and bootSfs from the
  dartRverse `binaries/epos_mac.zip` (same `epos` binary); EPOS example data
  from github.com/EvolBioInf/epos (`testF.dat`, `fig2bFolded.sfs`,
  `fig2bUnfolded.sfs`) for the input format; macOS
- Baseline: tests/testthat/test-gl.run.epos.R (11 tests, 28 expectations;
  9 tests need EPOS_DIR, the bootstrap test also needs bootSfs)

## Verdict

**Standards: Needs work** — failures stop with an empty message or an
unrelated R error, SilicoDArT data is accepted, "Output written" prints at
`verbose = 0`, and there is no way to set a seed.
**Spec: Rework** — every setting other than the default `minbinsize = 1,
folded = TRUE` sends EPOS a wrongly labelled SFS: `folded = FALSE` runs as
a folded SFS of twice the sample size, `minbinsize = 2` shifts every class
down by one, and `minbinsize = 0` (the documented alternative to `L`)
always fails. `upper` and `lower` are never passed to epos2plot.

What works well: with the defaults the SFS matches the EPOS example format
(classes 1..n/2 for n sequences) and the diagnostics (sites, likelihood,
d^2, observed/expected SFS) are parsed correctly; the run folder is a
`tempfile()` and `cleanup` removes only that folder.

## Findings

**F1 [HIGH, confidence: high] — `folded = FALSE` is run as a folded SFS with twice the sample size (DOC5)**
`R/gl.run.epos.R:167-178` — the unfolded SFS from `gl.sfs` (classes
1..2n) is sent without `-U`, the option the documentation names.
Failure scenario (baseline test 7): 30 individuals give classes 1..60;
EPOS reads them as a folded SFS of 120 sequences and returns Ne 15,000.
Sent correctly (classes 1..59, `-U`), EPOS returns Ne 9,440.
Proposed change: add `-U` and send classes 1..(2n - 1) (the class 2n holds
loci fixed for the derived allele).
**Consequence: estimates change for every call with `folded = FALSE`.**

**F2 [HIGH, confidence: high] — `minbinsize >= 2` shifts the SFS by one class per dropped class (DOC5)**
`R/gl.run.epos.R:168` — the SFS is written as `r = 1:length(sfs)`, so
after `gl.sfs` drops classes below `minbinsize` the counts are relabelled
from 1.
Failure scenario (baseline test 6): with `minbinsize = 2` doubletons are
sent as singletons (EPOS sees classes 1..29): Ne 29,000 from 77,300
generations. Sending all classes and excluding singletons with EPOS's own
`-x 1` gives Ne 31,100 from 83,100 generations.
Proposed change: write each count with its true class number and exclude
classes below `minbinsize` with `-x`.
**Consequence: estimates change for every call with `minbinsize >= 2`.**

**F3 [HIGH, confidence: high] — `minbinsize = 0` (zero class instead of `L`) always fails (DOC5)**
`R/gl.run.epos.R:168,174` — the zero class is written as class 1 and `-l`
is dropped, so EPOS finds neither a zero class nor a sequence length.
Failure scenario (baseline test 5): possums.gl[1:30, 1:200] (24
monomorphic loci) stops with "ERROR[epos]: Please include either the
zero-clas in the SFS or enter the sequence length", then R's "'names'
attribute [4] must be the same length as the vector [1]". Written as
class 0, EPOS runs. Note: the zero class counts only the monomorphic loci
kept in the genlight (24 here against L = 1e5), so it rarely stands in for
`L` with DArT data; the documentation should say so.
Proposed change: write the zero class as class 0 (fixed together with F2
by writing true class numbers); document when it is meaningful.
**Consequence: `minbinsize = 0` returns results instead of an error.**

**F4 [MEDIUM, confidence: high] — `upper` and `lower` are ignored (DOC5)**
`R/gl.run.epos.R:192` — `epos2plot` is called without `-u`/`-l`.
Failure scenario: `boot = 100, upper = 0.75, lower = 0.25` returns the
95% limits (epos2plot defaults 0.975/0.025), labelled by the user as 50%
limits. dartR Shiny exposes both arguments.
Proposed change: pass `-u upper -l lower` to epos2plot.
**Consequence: `low`/`high` change for calls with non-default `upper` or
`lower`.**

**F5 [MEDIUM, confidence: high] — command built with empty options (FS5, DOC5)**
`R/gl.run.epos.R:174-178` — `L = NULL` gives `-l ` with no value,
`u = NULL` gives `-u ` with no value, and `other.options` is pasted with
no leading space.
Failure scenario (baseline test 8): each of the three stops with EPOS
"efopen(...) failed", then R's "'names' attribute [4] must be the same
length as the vector [1]". The documentation says that without `u` the
EPOS default (5e-9) is used; the example option style ("-m, -x") fails.
Proposed change: omit `-u` when `u` is NULL; stop with a message when `L`
is NULL and the SFS has no zero class; add a space before `other.options`.

**F6 [MEDIUM, confidence: high] — failures stop with empty or unrelated messages (FS5, VRB2)**
`R/gl.run.epos.R:146-161,186-200` — missing binaries print with `cat()`
then `stop()` with no message; exit statuses of bootSfs/epos/epos2plot
are never checked, so any EPOS error surfaces as the `names` error above.
`chmod 777` runs on Linux only.
Failure scenario (baseline test 1): a wrong `epos.path` gives `Error:`
with an empty message (the reason is on stdout, lost in dartR Shiny).
Proposed change: `stop(error(...))` naming the missing files and
`gl.download.binary("epos")`; check each exit status and stop with EPOS's
message; `Sys.chmod(0755)` on Linux and macOS.

**F7 [MEDIUM, confidence: high] — SilicoDArT data is accepted (FS4)**
`R/gl.run.epos.R:128` — `datatype` is computed and never used.
Failure scenario (baseline test 9): `testset.gs` runs and returns a
history built from presence/absence scores.
Proposed change: stop on SilicoDArT input.
**Consequence: SilicoDArT input errors.**

**F8 [MEDIUM, confidence: high] — an unnamed user `sfs` fails (DOC5)**
`R/gl.run.epos.R:252-253` — class numbers are read from names like `d1`.
Failure scenario (baseline test 10): `sfs = c(6, 6, 3, ...)` (a plain
vector, as from another program) stops with "arguments imply differing
number of rows: 0, 30" after EPOS has run.
Proposed change: take class numbers from `d<k>` names when present,
otherwise from `minbinsize` upwards; document the format.

**F9 [LOW, confidence: high] — runs cannot be made reproducible (principle: reproducibility)**
`R/gl.run.epos.R:172,178` — bootSfs and epos seed from the system clock.
Failure scenario (probe): two identical calls with `boot = 5` return
different histories.
Proposed change: add `seed = NULL`; when set, pass `-s seed` to epos and
bootSfs.
**Consequence: new argument; default behaviour unchanged.**

**F10 [LOW, confidence: high] — messaging ignores verbosity; SFS plot ignores `plot.theme` (FS3, VRB1, VRB3, PLT1)**
`R/gl.run.epos.R:122-124,247,261` — `utils.flag.start` is commented out,
"Output written to" prints at every verbosity, EPOS warnings are shown
unconditionally, `sfs_plot` uses `theme_bw()`.
Failure scenario (baseline test 4): `verbose = 0` prints the output path.
Proposed change: `utils.flag.start`; gate messages at `verbose >= 2`;
`sfs_plot` uses `plot.theme`.

**F11 [LOW, confidence: high] — documentation gaps and mismatches (DOC1, DOC5, DOC7 (proposed rule))**
`R/gl.run.epos.R:1-86` — no `@name`/`@title`/`@family`; `@author` has no
Author(s) part; `outfile` default documented as 'genepop.gen' (it is
'epos.out'); `@return` lists history as "generation, median, low and high"
(the order is generation, low, median, high); the EPOS description
paragraph is repeated; `-l` and the zero class are not explained;
typos (programm, paramter, "provide via the sfs parameter (see below)").
Failure scenario: a user looks for `genepop.gen` or reads the columns in
the documented order.
Proposed change: fix the items above and add the zero-class and
`other.options` notes.

**F12 [INFO, confidence: high] — EPOS doubles the middle class of a folded SFS on input**
EPOS reports 180 polymorphic sites for 176 loci: the class n/2 (4 loci)
is counted twice. EPOS's own example `testF.dat` shows the same (24,148
sites in the file, 24,384 reported). This is EPOS's convention, not a
wrapper defect. No change.

## Proposed changes

1. Send EPOS the SFS with its true class numbers: class 0 when kept,
   `-x` for classes below `minbinsize`, `-U` and classes 1..(2n - 1) for
   `folded = FALSE` (F1, F2, F3).
   **Consequence: estimates change for `folded = FALSE` and for
   `minbinsize >= 2`; `minbinsize = 0` returns results instead of an
   error. Default calls are unchanged.**
2. Pass `upper`/`lower` to epos2plot (F4).
   **Consequence: `low`/`high` change for non-default `upper`/`lower`.**
3. Command and failure handling: omit `-u` when NULL, require `L` without
   a zero class, space before `other.options`, clear errors for missing
   binaries and EPOS failures, `Sys.chmod(0755)` (F5, F6).
4. Reject SilicoDArT input (F7).
   **Consequence: SilicoDArT input errors.**
5. Accept an unnamed user `sfs` (F8).
6. New `seed` argument passed to epos and bootSfs (F9).
   **Consequence: new argument; default behaviour unchanged.**
7. Verbosity and `sfs_plot` theme (F10).
8. Documentation fixes (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour vs roxygen with the real binaries (defaults,
  `minbinsize` 0 and 2, `folded = FALSE`, `L`/`u` missing,
  `other.options`, `boot = 5`, `upper`/`lower`, SilicoDArT, unnamed
  `sfs`, wrong path, unknown method, `verbose = 0`, cleanup) — run
- Correct-input comparison: EPOS run directly on correctly labelled SFS
  files (`-x 1`; `-U` with classes 1..59; class 0) — run
- Callers (API3): dartR Shiny (`dartr2shiny/shiny_fun/Fun_gl.run.epos.R`)
  passes `minbinsize`, `folded`, `upper`, `lower`, `method`, `depth`,
  `outpath`, `outfile`, `plot.theme`; no sibling dartR.* package calls it
  (dartRverse only documents it in `gl.download.binary`) — run
- Windows and Linux: SKIPPED — not run; Windows paths (`.exe`, dlls) read
  only
- `method = "exhaustive"`: not run (command construction read only)
- DAT1-DAT4, FS8: not applicable — no genlight returned
- DAT6 (FBM): SKIPPED — no FBM fixture; `gl.sfs` densifies, to be covered
  in the `gl.sfs` review
- Missing data in the SFS: same issue as gl.run.stairway2 F12, left to the
  `gl.sfs` review
- Google Group / GitHub issues: not searched in this session

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

Branch `review-epos` (from origin/dev 1e90a1e). Evidence:
`tests/testthat/test-gl.run.epos.R`, 11 tests, 37 expectations, all passing
with EPOS, epos2plot and bootSfs (dartRverse `epos_mac.zip`, macOS);
without EPOS_DIR 3 tests run and 8 skip. Every changed assertion is tagged
with its approved change.

- 1 applied: the SFS is written with true classes (class 0 when
  `minbinsize = 0`; all classes from 1 with `-x` below `minbinsize`; `-U`
  and classes 1..2n-1 when unfolded). Each result equals EPOS run directly
  on a correctly labelled file (seed 1): `minbinsize = 2` Ne 31,100 /
  83,100 generations (was 29,000 / 77,300); `folded = FALSE` Ne 9,440 /
  37,100 (was 15,000); `minbinsize = 0` Ne 1.5e7 with 24 monomorphic
  sites (was an error). Side effect: the returned `sfs` no longer carries
  `d1..dn` row names (values identical).
- 2 applied: `-l lower -u upper` passed to epos2plot; with `boot = 5,
  seed = 3`, 0.25/0.75 limits lie inside the 0.025/0.975 limits (test 11).
- 3 applied: `-u` omitted when NULL; `L` required unless `minbinsize = 0`;
  programs called by full path with `system2()`, arguments as a vector (so
  `other.options` needs no leading space); a non-zero exit status stops
  with the program's message ("epos failed (exit status 139)" for
  `-L abc`); missing binaries named with a `gl.download.binary` pointer;
  `Sys.chmod(0755)`.
- 4 applied: SilicoDArT errors (test 9).
- 5 applied: an unnamed `sfs` gives the same history as the computed one
  (test 10).
- 6 applied: `seed` passed as `-s` to epos and bootSfs; two seeded
  bootstrap runs identical (test 11).
- 7 applied: `utils.flag.start`; messages at `verbose >= 2`, EPOS stderr
  at `verbose >= 3`; `sfs_plot` uses `plot.theme`. `verbose = 0` prints
  nothing (test 4).
- 8 applied: roxygen rewritten; `devtools::document()` also rewrote seven
  unrelated Rd files (stale "population structure" family links already
  on dev); those were reverted, only `gl.run.epos.Rd` and
  `gl.run.stairway2.Rd` (new shared family link) are committed.
- Default calls unchanged: with seed 1 the history and diagnostics are
  identical to the pre-change run.
- `devtools::check()`: 0 errors; install WARNING from packages built under
  R 4.4.3; 3 NOTEs (hidden `.git`, time check, NEWS parse) present before
  this change.
- NEWS entry added. Callers: dartR Shiny passes `minbinsize`, `folded`,
  `upper`, `lower` (now effective); it does not read `sfs` row names.
- Windows and Linux not run.
- PR: #105

## Machine block

```json
{
  "function": "gl.run.epos",
  "package": "dartR.popgen",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "1e90a1e",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 4},
    {"id": "F8", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "principle: reproducibility", "status": "approved", "change": 6},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 7},
    {"id": "F11", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 8},
    {"id": "F12", "severity": "INFO", "confidence": "high", "rule": "none (EPOS convention)", "status": "no-change", "change": null}
  ],
  "coverage_skipped": ["Windows/Linux not run", "method exhaustive not run", "DAT6: no FBM fixture", "Google Group / issues not searched"],
  "status": "pr-open",
  "pr": 105
}
```
