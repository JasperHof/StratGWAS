# End-to-end tests on the example data shipped in inst/extdata.
# These take ~30 seconds, so they are skipped on CRAN.

test_that("compute_gencov returns complete, symmetric estimates", {
  skip_on_cran()
  w <- example_workflow()
  g <- w$gencov
  K_tot <- ncol(w$strata$multi) - 2

  expect_equal(dim(g$gencov), c(K_tot, K_tot))
  expect_false(anyNA(g$gencov))
  expect_equal(g$gencov, t(g$gencov))
  expect_equal(unname(diag(g$gencor))[g$hers > 0], rep(1, sum(g$hers > 0)))
})

test_that("jackknife estimates are filled in when jack = TRUE", {
  skip_on_cran()
  g <- example_workflow()$gencov
  expect_length(g$jack_ests, 10)
  expect_false(any(vapply(g$jack_ests, anyNA, logical(1))))
})

test_that("compute_gencov aligns SNPs that differ between strata and LD scores", {
  skip_on_cran()
  w <- example_workflow()
  K_tot <- ncol(w$strata$multi) - 2

  ss_list <- lapply(seq_len(K_tot), function(k)
    read.table(paste0(w$outfile, ".pheno", k), header = TRUE))
  lds <- read.table(paste0(w$outfile, ".ldscores"), header = TRUE)

  # Drop different SNPs from different strata, and give the LD scores other column names
  ss_list[[2]] <- ss_list[[2]][-(1:100), ]
  ss_list[[3]] <- ss_list[[3]][-(201:300), ]
  ss_list[[3]] <- ss_list[[3]][sample(nrow(ss_list[[3]])), ]   # scrambled order
  colnames(lds) <- c("SNP", "LDscore")

  g <- quietly(compute_gencov(w$strata, example_bed(), outfile = tempfile(),
                              ss_list = ss_list, lds = lds))
  expect_false(anyNA(g$gencov))
  # same estimates for the untouched pair as in the full run, within noise
  expect_equal(g$gencov[1, 1], w$gencov$gencov[1, 1], tolerance = 0.1)
})

test_that("transform keeps missing phenotypes missing and gives sensible weights", {
  skip_on_cran()
  w <- example_workflow()
  tp <- w$trans$transformed_pheno

  expect_true(all(is.na(tp[tp$IID %in% w$na_ids, 3])))
  controls <- w$pheno[!is.na(w$pheno[, 3]) & w$pheno[, 3] == 0, 2]
  expect_true(all(tp[tp$IID %in% controls, 3] == 0))
  expect_false(anyNA(w$trans$weights))
  expect_equal(sum(w$trans$weights^2), 1, tolerance = 1e-8)
})

test_that("linear runs the final GWAS and writes one row per SNP", {
  skip_on_cran()
  w <- example_workflow()
  expect_no_error(quietly(linear(w$trans, example_bed(), w$outfile)))

  res <- read.table(paste0(w$outfile, ".assoc"), header = TRUE)
  n_snps <- nrow(read.table(paste0(example_bed(), ".bim")))
  expect_equal(nrow(res), n_snps)
  expect_true(all(res$P >= 0 & res$P <= 1))
})
