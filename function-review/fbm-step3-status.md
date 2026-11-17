# FBM step 3: where we are (hand-over, 16 Nov 2026, late)

Branch `dev_fbm` of dartR.base (dartR.sim also has a `dev_fbm` branch, with
the WP0 helpers and `gl.sim.mutate`). Step 2 (WP0-WP6) is complete; see
`fbm-step2-plan.md` for its plan, measurements and decisions. Step 3 items
come from Luis's summary (14 Nov 2026), "Proposals, step 3".

To resume on another computer: `git pull` on `dev_fbm` in dartR.base, then
continue with "Next" below.

## Step 3 items

| # | Item | Status |
|---|---|---|
| 1 | Count from disk (reports and filters) | **In progress**, see below |
| 2 | Read straight to disk (PLINK, VCF, DArT, CSV) | Not started. PLINK `.bed` is easy via bigsnpr |
| 3 | Save/load without converting | Not started; low value (gen form is 4x smaller, conversion now fast) |
| 4 | Filter without copying (class change) | Not started; needs Bernd's decision |
| 5 | Kinship and distances from disk | Not started |
| 6 | Read counts and likelihoods on disk | Not started |
| 7 | Examples in both storage types in CI | Runner exists (`run_examples_fbm.R`, WP6); not in CI |
| 8 | New bigsnpr analyses | Not started |

## Count from disk: method

Each function keeps its original per-locus / per-individual expressions;
they are applied to blocks of loci read from the backing file, so the
results are identical by construction (no rederived formulas). Helpers in
`R/utils.fbm.r`:

- `.gl_colwise(x, FUN, rows = NULL)`: per-locus statistic, blocks
  concatenated (or `cbind` if FUN returns a matrix with loci in columns).
  `rows` = one population's individuals, instead of `seppop()` copies.
- `.gl_colblocks(x, FUN, rows = NULL)`: the list of per-block results, for
  callers that combine them themselves (lists of matrices, added counts).
- `.gl_rowwise(x, FUN)`: per-individual sums over loci, added over blocks.
- `.gl_n_na(x)`: number of missing genotypes.
- `.gl_pop_genotype_counts(x)`: counts of 0/1/2 per population and locus.
- Block size: `.fbm_block_size()` = bigstatsr `block_size()` / 16, set by
  `options(bigstatsr.block.sizeGB)`. In-memory objects take the old path
  (`FUN(as.matrix(x))`), unchanged.

Tests: `tests/testthat/test-fbm-counts.R`. Every converted function is one
entry of the `calls` list; the test checks identical value and console
output for gen and FBM (block size 1 Gb and one locus per block), and that
no call decodes the whole FBM (`.fbm_as_matrix` mocked to fail). For
statistical changes also compare the new code with the old version from git
(`git show <commit>:R/<file>`, `sys.source()` into an env whose parent is
the namespace), including seeded bootstraps.

Pitfall: in a test with `local_mocked_bindings()`, every `on.exit()` must use
`add = TRUE`, or the mock leaks into later test files.

## Count from disk: done (all pushed)

Commits 69d68f6 and the HWE commit after it:

- Per locus / per individual: `gl.allele.freq`, `gl.filter.callrate`
  (ind), `gl.filter.heterozygosity`, `gl.report.allna`,
  `gl.report.callrate`, `gl.report.maf`, `gl.report.monomorphs`,
  `gl.report.rdepth`, `gl.report.reproducibility`, `gl.report.taglength`.
- Per population: `gl.report.heterozygosity` (pop and ind;
  `pop.het_fun()` split into `.het_byloc()` + `.het_summary()`),
  `gl.report.basics`, `gl.tree.nj`, `utils.basic.stats` (split into
  `.basic_pop_stats()` + `.basic_stats_combine()`), `gl.report.fstat`
  (pairs and bootstrap from per-locus values), `gl.report.hwe`,
  `gl.filter.hwe`, `gl.hwe.pop`.
- `gl.fst.pop`: Weir & Cockerham theta computed by dartR (identical to
  `StAMPP::stamppFst`, 30-200x faster, FBM no longer converted). bigsnpr's
  `snp_fst()` was rejected: its heterozygosity term assumes HWE (expected,
  not observed heterozygosity), which would change results. `nclusters` is
  ignored now.
- `gl.fixed.diff` reads from disk through `gl.allele.freq`, except its
  `test = TRUE` simulation (`gl.fdsim` uses `gl.keep.pop`).
- Earlier (step 2): compliance check, `utils.recalc.*`, `gl.alf`,
  `gl.filter.allna`, `gl.filter.monomorphs`, `utils.check.datatype`.

Speed (2,000 x 20,000): the block path is 2-3x faster than decoding the
FBM and up to 12x faster than the in-memory path (e.g. `gl.allele.freq`
1.9 s vs 22.5 s); memory about 200-240 Mb, independent of data size.

## Next (group C, remaining)

In this order (decodes / `seppop` copies to replace):

1. Done (17 Nov, f15888b): `gl.report.pa`, `gl.filter.pa`,
   `utils.pa.Chao` (core moved to `.pa_chao()`, built from per-locus
   terms). Also fixed: with `x2`, `test.asym` split the pooled sample at
   the first population of `x`.
2. Done (17 Nov): `gl.report.diversity`. One block pass stores per
   population and locus the mean dosage and genotype counts; pooled pair
   frequencies are (sum of dosages) / (number of genotypes) of both
   populations. Equal to `colMeans(rbind(m1, m2))` bit for bit for pooled
   samples up to 400 individuals (checked exhaustively); `colMeans` divides
   in long double, so above ~1,000 individuals a 1-ulp difference is
   theoretically possible.
3. Done (17 Nov): `gl.report.allelerich` (counts from
   `.gl_pop_genotype_counts()`, now one pass for all populations; the
   bootstrap decodes one population at a time, in blocks for FBM).
4. Done (17 Nov): `gl.test.heterozygosity` (bootstrap decodes one
   population at a time, in blocks for FBM; identical).
5. Done (17 Nov): `gl.report.shannon`, from per-individual counts of
   dosages 1/2 and per-population histograms of pooled dosages. Not
   bit-identical to the old version: the Hill sums are formed per
   distinct value instead of locus by locus (long-double summation order),
   max relative difference 5e-14; gen and FBM identical.
6. Optional: `gl.fdsim` (used by `gl.fixed.diff(test = TRUE)`).

Out of scope for item 1: group D (distances, kinship, replicates:
`gl.dist.ind`, `gl.dist.phylo`, `gl.report.replicates(.KING)`,
`gl.report.contamination`) belongs to item 5; group E (about 25 `gl2*`
exports and `gl.write.csv`) needs block-wise writing. `gl.report.hwe` and
`gl.report.heterozygosity` with `subsample.pop = TRUE` / bootstrap still
decode one population at a time by design.

## Notes

- Not run yet: `R CMD check` (deferred by Bernd; run once later). Last
  check (after WP3) had only environmental issues (ggtern not installed:
  `gl.report.excess.het` example; packages built under R 4.3.3).
- Full testthat suite passed before each push (about 20 minutes).
- `tests/testthat/Rplots.pdf` changes on every test run; never commit it.
- Downstream: dartR.popgen and dartR.spatial call `gl.fbm2gen` for Fst,
  AMOVA and IBD; not changed yet. Run their tests after changes to shared
  helpers (`R_LIBS` with a scratch library if the installed versions are old).
