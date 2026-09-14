# Characterization tests for utils.plink.run
# Baseline snapshotted before review (infrastructure wave, dev at
# ddaed27). Assertions marked [approved diff] were flipped in Phase C.
# The composed command string is the contract; no PLINK binary needed.

test_that("the composed command is well formed", {
  td <- tempdir()
  o <- capture.output(cmd <- suppressWarnings(utils.plink.run(dir.in = td,
        plink.cmd = "nonexistent_exe_999", syntax = "--file hapmap1",
        verbose = 0)))
  # [approved diff I11] baseline: the default plink.path produced a
  # literal "path/" prefix, and a missing space glued syntax onto --out:
  # "path/nonexistent_exe_999 --file hapmap1--out hapmap1".
  expect_false(grepl("hapmap1--out", cmd))  # [approved diff I11]
  expect_false(grepl("^path/", cmd))       # [approved diff I11]
})
