# FBM step 3: where we are (hand-over, 17 Nov 2026)

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
| 2 | Read straight to disk (PLINK, VCF, DArT, CSV) | Redefined (17 Nov): audit readers for dense intermediates, see "Storage model" |
| 3 | Save/load without converting | **Dropped** (17 Nov), see "Storage model" |
| 4 | Filter without copying (class change) | Not started; needs Bernd's decision |
| 5 | Kinship and distances from disk | **Done** (17 Nov), see below |
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
6. Checked (17 Nov), no change needed: `gl.fdsim` and
   `gl.fixed.diff(test = TRUE)` already read FBM in blocks through
   `gl.keep.pop`, `gl.filter.monomorphs` and `gl.allele.freq` (no full
   decode, results identical for both storage types); `gl.fdsim` added to
   `test-fbm-counts.R`.
7. Done (17 Nov, not in the original list): `gl.filter.maf(by.pop = TRUE)`
   and its per-population plots (`seppop()` copies replaced by one pass,
   `.pop_maf()`; returned object and plot data identical).

Group C is complete. Remaining `seppop`/decodes by design: per-population
bootstraps, `gl.report.ld.map` (LD needs each population's matrix),
`gl.smearplot` (plots every genotype); polyploid functions (FBM refuses
polyploid data). Item 5 (group D) followed; see its section below.

Out of scope for item 1: group D (distances, kinship, replicates:
`gl.dist.ind`, `gl.dist.phylo`, `gl.report.replicates(.KING)`,
`gl.report.contamination`) belonged to item 5 (done); group E (about 25 `gl2*`
exports and `gl.write.csv`) needs block-wise writing. `gl.report.hwe` and
`gl.report.heterozygosity` with `subsample.pop = TRUE` / bootstrap still
decode one population at a time by design.

## Storage model (decided 17 Nov 2026, Bernd)

The genlight (gen-backed dartR object, saved with `gl.save()`) is the
storage and exchange format; the FBM is the working copy. This is how the
code already works: `gl.save()` converts an FBM object with `gl.fbm2gen()`
(block-wise, WP2) and saves the genlight; `gl.load(fbm = TRUE)` and the
readers' `fbm = TRUE` convert with `gl.gen2fbm()` (block-wise, WP1);
`gl.load()` refuses an FBM object saved with `saveRDS()`.

Reasons not to store or pass on FBM objects: the `.bk` file holds only
genotypes (names, populations, metrics and history are in the R object,
so two files must stay together); the object refers to the backing file
by absolute path, and temporary backing files are deleted automatically
(WP4); the FBM takes 1 byte per genotype against 2 bits (plus NA
positions) in the genlight, about 4x larger.

A save cannot skip the in-memory genlight (an `.rds` serialises an R
object), but that is the packed form (about 250 Mb at 2,000 x 200,000)
built block by block. Hence item 3 is dropped, and item 2 becomes: check
that the readers (`gl.read.dart`, `gl.read.vcf`, `gl.read.plink`, CSV)
build no dense intermediate (8+ bytes per genotype) before the genlight;
write blocks straight to an FBM only where the genlight itself would not
fit.

## Item 5: distances and kinship from disk (done, 17 Nov)

- `gl.dist.ind` (SNP, all methods, via `utils.dist.ind.snp`) and
  `gl.report.replicates`: `.gl_pair_counts()` (utils.fbm.r) adds per pair
  the loci scored in both (L) and those differing by 1 (n1) and 2 (n2)
  allele copies, in a compiled loop over blocks of loci (dist-order integer
  vectors, 12 bytes per pair). All distances and the replicate proportions
  follow from these counts with the old expressions: identical results.
  2,000 x 20,000: 94 s FBM vs 330 s in memory (old loop).
- `gl.report.replicates.KING` (calls input): missing data and allele
  frequencies counted per block; `.king_calls9()` reads the FBM block by
  block. Likelihood and count inputs use `other$geno.metrics` (in memory,
  item 6).
- `gl.report.contamination`: `.gl_reader()` (blocks of loci, single rows,
  selected loci); tier 1 and the kinship cross-product in two passes over
  blocks; tiers 3 and 4 read only suspects' rows. Identical to the old
  version with reference BLAS (five data sets) and between storage types.
- `.gl_count_mean()`: count means exactly as `rowMeans()`/`mean()` (long
  double division); a plain division differs by 1 ulp for about 1 in 6,000
  counts above 2,048.
- `gl.dist.pop`: already from disk via `gl.allele.freq`; now in the tests.
- Not converted: `gl.dist.phylo` (DNA alignment for ape, via `gl2fasta`,
  group E); SilicoDArT distances (`utils.dist.binary`; FBM is SNP only).
- Tests: entries in `test-fbm-counts.R` (dist.ind x5, dist.pop x2,
  replicates, KING, contamination), unit tests of `.gl_count_mean` and
  `.gl_pair_counts`, and a gen/FBM test with read counts and the depth
  pattern in `test-gl.report.contamination.R` (its `depth.host()` fixture
  now builds a dartR object).
- Benchmarks: keep them small from now on (Bernd, 17 Nov).

Next: item 2 (reader audit) or group E (`gl2*` exports written block by
block); also downstream dartR.popgen/dartR.spatial `gl.fbm2gen` calls.

## Downstream audit (17 Nov 2026)

About 120 `as.matrix()`/`seppop()`/`gl.fbm2gen()` calls in dartR.popgen,
.spatial, .captive, .sexlinked and .sim; about half act on distance
matrices or data frames, not genotypes. The genotype decodes:

A. Easy, identical by construction (per-locus or per-sex/pop counts):
   dartR.sim `gl.He`, `gl.Ho`, `gl.report.nall`, `gl.sim.WF.run` (freqs),
   `utils.sims.r` (639-661, 3347); dartR.sexlinked `gl.filter.sexlinked`,
   `gl.keep.sexlinked`, `gl.report.sexlinked` (`t(as.matrix(x))` split by
   sex), `gl.infer.sex`; dartR.popgen `gl.sfs`, `gl.select.panel`
   (allele freqs per pop; `gl.fbm2gen` + `seppop` for dapc: decode only the
   pair of populations), `gl.check.panel`, `gl.check.future.panel`.
B. Distances, from base helpers: dartR.spatial `gl.kosman` (sum |a - b| and
   loci in common = `.gl_pair_counts` n1 + 2 n2 and L: exact), `gl.ibd`
   (`gl.fbm2gen` for StAMPP Fst/NeisD between pops: Fst from
   `gl.fst.pop`, identical to stamppFst; Nei D from allele freqs;
   euclidean `dist(as.matrix(x))` from pair counts with dist()'s NA
   scaling), `gl.propShared` (already via `gl.dist.ind`);
   dartR.popgen `gl.amova` (`gl.fbm2gen` for `stamppNeisD` between
   individuals: check its formula; products of 0/0.5/1 are exact).
C. Full matrix needed by an external routine: `gl.grm` (captive) and
   `gl.grm2` (spatial) via `rrBLUP::A.mat`; `gl.relatedness`,
   `gl.diagnostics.relatedness`, `gl.parentage.power`,
   `gl.report.errorrate`, `gl.report.parent.offspring` (captive);
   `gl.opt.ind.all`, `gl.ld.haplotype` (per pop). Decide case by case;
   where a decode stays, make it explicit (message with the memory
   estimate) instead of a silent `as.matrix()`.
D. Writers for external programs (like group E): `gl2colony`,
   `gl.run.EMIBD9`, `gl.run.eems`, `gl.run.bayescan`, `gl.run.mongrail`,
   `gl.nhybrids`, `gl.run.assignpop`, `gl.outflank`, structure utils.
   Leave for now.
E. One individual only (fine): `gl.assign.*` (`as.matrix(x[unknown, ])`).

Shared helpers: dartR.spatial and dartR.captive (zzz.r), and two popgen
functions, keep local copies of `.fbm_or_null`. Plan (agreed): export the
base helpers as `utils.*` functions (`utils.fbm.has`, `utils.fbm.colwise`,
`utils.fbm.rowwise`, `utils.fbm.reader`, `utils.fbm.pair.counts`,
`utils.fbm.count.mean`) and use them downstream instead of copies or
`:::`. Downstream DESCRIPTIONs then need a minimum dartR.base version
(dependency change: confirm with Bernd before bumping).

### Decisions and progress (17 Nov 2026)

Bernd: create `dev_fbm` branches in popgen, spatial, captive, sexlinked
(done, from `dev`); minimum `dartR.base (>= 1.2.7.9000)`; do A and B now,
C later together with Bernd, skip D and E.

Done (A and B), all identical to the gen-backed results, tested in each
package's `tests/testthat/test-fbm.R` (or `test-fbm-counts.R` in sim)
for both storage types, one-locus blocks and no full decode
(`.fbm_as_matrix` mocked in dartR.base):

- dartR.base: exported `utils.fbm.has`, `.colwise`, `.colblocks`,
  `.rowwise`, `.count.mean`, `.reader`, `.blocks`, `.pair.counts`
  (`R/utils.fbm.export.r`).
- dartR.sim: `gl.He`, `gl.Ho`, `gl.report.nall`, `utils.sim.alf` (its
  big_counts path replaced: it could differ from colMeans by 1 ulp),
  `fst_hudson`, `inbreeding_real`, `gl.sim.WF.run(real_freq)`;
  `utils.sim.is.fbm` delegates to `utils.fbm.has`. Note: helpers used on
  parallel workers (`gl.report.nall`) must not call new internal functions;
  workers load the installed package.
- dartR.sexlinked: report/filter/keep (counts by sex in one pass),
  `gl.infer.sex` (per-individual counts; gametolog loci by reader).
- dartR.popgen: `gl.sfs`, `gl.select.panel` (dapc converts only the pair
  of populations), `gl.amova` (`.neis_d_ind_fbm`: StAMPP Nei distance
  between individuals rebuilt from exact pair counts; terms are multiples
  of 0.25).
- dartR.spatial: `gl.kosman` (pair counts), `gl.ibd` (Fst via
  `gl.fst.pop`, Nei D via StAMPP on a per-population frequency table,
  Euclidean from pair counts with dist()'s NA scaling); local
  `.fbm_or_null` removed from zzz.r.

Group C, GRM (decided 17 Nov with Bernd: block-wise for A.mat's default
only, explicit conversion with a warning for the other options; EJ
shrinkage is not the default and no dartR function asks for it): done.
`utils.fbm.grm` (A.mat default in two passes; equal to A.mat up to the
rounding of the summed cross-product, bit for bit with reference BLAS) and
`utils.fbm.amat` (dispatcher) in dartR.base; used by `gl.grm` (captive)
and `gl.grm2` (spatial). `gl.check.panel(parameter = "relatedness")` uses
its own VanRaden GRM, now block-wise for the full data (identical to the
previous version). Gen-backed objects still use rrBLUP::A.mat unchanged.

Group C, rest (decided 17 Nov with Bernd: option 2 = decode with a
message and a memory limit; convert the exact cases): done.
- `utils.fbm.decode()` (dartR.base): as.matrix for gen-backed; for FBM
  reports the memory (verbose >= 2) and stops above
  `getOption("dartR_fbm_max_decode_gb", 4)` Gb (default 4 chosen by
  Claude; Bernd may change it). Also guards `utils.fbm.amat`'s conversion.
  Used by `gl.relatedness` (engine, 4 bytes/genotype),
  `gl.parentage.power`, `gl.opt.ind.all` (2 copies).
- Exact: `gl.report.parent.offspring` (pair counts n2, L),
  `gl.report.errorrate` (replicate rows, counts per block), all-NA
  pre-filters of `gl.relatedness` and `gl.diagnostics.relatedness`.
- `gl.ld.haplotype` needed nothing: it works per population through base
  statistics (block-wise) and `gl2plink` (group D, skipped).
- Not testable here: `gl.relatedness` and `gl.diagnostics.relatedness`
  need dartR.coancestry (not installed); their changes are the pre-filter
  counts and the decode call only.
- dartR.captive: local `.fbm_or_null` removed from zzz.r
  (`gl.filter.parent.offspring` uses `dartR.base::utils.fbm.has`). `gl.check.panel` / `gl.check.future.panel` (popgen, group A
list) not converted yet: they compare a panel with the full data set and
decode both; next to do.

## Notes

- Not run yet: `R CMD check` (deferred by Bernd; run once later). Last
  check (after WP3) had only environmental issues (ggtern not installed:
  `gl.report.excess.het` example; packages built under R 4.3.3).
- Full testthat suite passed before each push (about 20 minutes).
- `tests/testthat/Rplots.pdf` changes on every test run; never commit it.
- Downstream: dartR.popgen and dartR.spatial call `gl.fbm2gen` for Fst,
  AMOVA and IBD; not changed yet. Run their tests after changes to shared
  helpers (`R_LIBS` with a scratch library if the installed versions are old).

## Item 2 (reader audit), in progress (17 Nov 2026, end of session)

Agreed with Bernd: start with gl.read.PLINK, then gl.read.dart; with
fbm = TRUE read straight into an FBM, with fbm = FALSE build the FBM
block-wise and convert with gl.fbm2gen at the end. Audit: all readers hold
dense intermediates (PLINK: write.SnpMatrix text + fread + as.matrix; DArT:
csv as text + two dense matrices in utils.dart2genlight + two decodes for
OneRatio; VCF: vcfR text; CSV: read.csv + t()).

gl.read.PLINK: done (17 Nov). `.plink_read_meta()` (bim/fam read exactly
as snpStats::read.plink) and `.plink_bed_to_fbm()` (SNP-major .bed read in
blocks into a CODE_012 FBM, count of allele.2, per-individual called/het
counts) in R/utils.fbm.r; the .bed route builds the dartR object with
new("dartR", fbm = ...), converts with gl.fbm2gen when fbm = FALSE.
Identical to the old reader (genotypes, metadata, row summary; both block
sizes); tests in test-gl.read.PLINK.R. Found and fixed: ploidy inferred
per individual from the largest genotype (1 or NA) in the old route.

Remaining steps (1 and 2 done):
1. In gl.read.PLINK, .bed route (also after plink.flags conversion): meta +
   .plink_bed_to_fbm; ind.metrics row summary as snpStats::row.summary
   (Call.rate = called/p, Certain.calls, Heterozygosity = het/called;
   check snpStats' values for an individual with no calls); build with
   new("dartR", fbm = G, ...) and run the rest of the pipeline; if
   !fbm, gl.fbm2gen at the end. .ped route unchanged.
2. Compare with the old function (git show HEAD:R/gl.read.PLINK.r) on a
   fixture with missing data, an all-NA individual, a fixed locus and
   chromosome 0 (snpStats::write.plink), both fbm settings; run
   tests/testthat/test-gl.read.PLINK.R.
3. gl.read.dart: done (17 Nov). Locus-metric columns read with
   colClasses (genotype columns skipped); `.dart_read_genotypes()` reads
   the genotype columns in chunks of rows from an open connection into a
   raw-score FBM (with the largest score for the 1-/2-row detection);
   utils.dart2genlight pairs/recodes in blocks of loci into the genotype
   FBM; OneRatio via .gl_colwise; fbm = FALSE converts with gl.fbm2gen.
   Identical to the old reader (2-row with/without metafile, 1-row, both
   fbm, chunk sizes 1 Gb and one row). Tests in test-gl.read.dart.R.
4. Next: gl.read.vcf (vcfR holds the whole VCF as text; options: chunks
   of variants via vcfR's nrows/skip, or another reader), then
   gl.read.csv (low priority).
