# Review: utils.plink.run (dartR.base)

## Provenance

- Model: claude-fable-5 (Claude Code); Skill: dartr-function-review
  v2.0.0; Base: dev at ddaed27; Branch: review-utils.plink.run. Reviewed in the
  infrastructure wave (seven files, one approval round, per-function
  PRs), under the standing member directive that utility functions
  are not for end users.
- Datasets: testset.gl subsets; constructed matrices; ggplot fixtures;
  composed-command inspection (no PLINK binary required).
- Family mode: analysis/infrastructure utility.
- Checks skipped: Google Group not searched (not available: no
  browser session); dartr2shiny not present in the local workspace.

## Verdicts

- **Standards: FAIL** — dead commented blocks; the on.exit/setwd
  workaround documented only by an inline aside.
- **Spec: FAIL** — the composed command is malformed twice (verified:
  "path/nonexistent_exe_999 --file hapmap1--out hapmap1"): (1) the
  default plink.path = "path" — documented as "plink is on the PATH,
  no path needed" — is pasted literally, producing a bogus "path/"
  prefix; (2) there is no space between syntax and "--out", gluing
  the last flag ("hapmap1--out"). Any call that does not both supply
  an explicit plink.path and end syntax with a trailing space runs a
  broken command. One caller in the family (gl.read.PLINK).

## Findings

### I11 — Command composition (MEDIUM) [escalation: composed command changes]

plink.path == "path" now means bare plink.cmd (PATH lookup as
documented); explicit paths joined with file.path; a space guaranteed
before --out. Callers that already worked around the gluing (trailing
space in syntax) produce a harmless double space.

### I12 — Tidy (LOW)

Dead commented blocks removed; @keywords internal (STAYS exported).

## Coverage

test-utils.plink.run.R — 2 assertions on the composed command
(baseline both defects). All pass pre-fix.

## Approval

All findings approved via the approval boxes (2026-09-01).

## Outcome

All findings applied (I11 command composition, I12 tidy). Suite: 2/2; the composed command is now well formed under the default path. PR #326.

## Addendum (2026-09-23): exit status and quoting

Follow-up to the `gl2plink` addendum (PR #418), which found that a failed
PLINK run was reported as success and that unquoted paths broke on spaces.
The same gaps were in `utils.plink.run` (and its caller `gl.read.PLINK`)
and in `gl2vcf`. Requested by Luis ("fix the same gaps elsewhere"),
2026-09-23. Reviewer: Claude (claude-opus-5-5), dartr-function-review
v2.0.0. Reproduced on `origin/dev` at def8e82 with PLINK 1.9.

**B1 [HIGH] failed PLINK run not detected (principle: fail loudly; FS5)**
`utils.plink.run` ran `system()` and ignored the exit status: with a stub
PLINK that exits with status 2 it returned normally, and `gl.read.PLINK`
then stopped with "PLINK did not produce a .bed file. Check that PLINK is
installed", hiding PLINK's reason. `gl2vcf` gave only an R warning and
wrote no VCF. Applied: both stop with an error quoting the last lines of
PLINK's output (stderr captured with `2>&1`); a command that cannot start
(wrong `plink.path`) is reported the same way.

**B2 [MEDIUM] unquoted paths (principle: platform-safe shell calls; DAT5)**
`gl.read.PLINK` on `my data/my file.ped` failed ("--out only accepts 1
parameter"); `gl2vcf` with `outpath` "my out" wrote no VCF. Applied:
`shQuote()` on the executable and `--out` in `utils.plink.run`, on the
`--file` name in `gl.read.PLINK`, and on every path in `gl2vcf`.

Also: `utils.plink.run` prints PLINK's output only at `verbose >= 3`
(stderr previously leaked at `verbose = 0`); outdated `build =` dropped.
`gl2vcf` keeps its reviewed choice of showing the log at `verbose >= 2`.

Evidence: all four reproductions pass after the change (stub failure errors
in both functions; `gl.read.PLINK` reads `my file.ped`; `gl2vcf` writes the
VCF under "my out"). Tests: test-utils.plink.run.R 8/8 (stub failure,
spaces, verbose gating; the missing-executable baseline flipped from a
returned command to an error), test-gl2vcf.R 49/49 (new stub-failure test;
new real-PLINK space test gated on PLINK19_DIR), test-gl2plink.R 17/17.
test-gl.read.PLINK.R 25/26: the one failure ("silent at verbose = 0 on the
.bed path") fails identically on unmodified origin/dev and is not touched
by this change.

```json
{"function": "utils.plink.run", "package": "dartR.base", "family_mode": "analysis",
 "commit": "ddaed27", "skill_version": "2.0.0",
 "verdict_standards": "FAIL", "verdict_spec": "FAIL",
 "findings": [{"id": "I11", "severity": "MEDIUM", "rules": ["spec"], "loc": "R/utils.plink.run.r command", "status": "applied"},
  {"id": "I12", "severity": "LOW", "rules": ["STY"], "loc": "R/utils.plink.run.r", "status": "applied"}],
 "datasets": ["testset.gl", "constructed"],
 "baseline_test": "tests/testthat/test-utils.plink.run.R",
 "pr": 326}
```
