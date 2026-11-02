data("pheno", package = "StratGWAS")
data("strat_cont", package = "StratGWAS")
data("strat_cat", package = "StratGWAS")

test_that("stratify rejects invalid input", {
  bad <- pheno; bad[1, 3] <- 2
  expect_error(stratify(bad, strat_cont = strat_cont), "only 0")

  dup <- pheno; dup[2, 2] <- dup[1, 2]
  expect_error(stratify(dup, strat_cont = strat_cont), "Duplicate")

  expect_error(stratify(pheno), "stratification variable")
  expect_error(stratify(pheno, strat_cont = strat_cont, K = 1), "at least 2")
})

test_that("continuous stratification gives K equally sized, mutually exclusive strata", {
  K <- 4
  s <- quietly(stratify(pheno, strat_cont = strat_cont, K = K))
  multi <- apply(s$multi[, paste0("continuous_1_", 1:K)], 2, as.numeric)

  is_case    <- pheno[, 3] == 1
  is_control <- pheno[, 3] == 0

  # controls are 0 in every stratum
  expect_true(all(multi[is_control, ] == 0))
  # every case with a stratification value is in exactly one stratum
  expect_true(all(rowSums(multi[is_case, ] == 1, na.rm = TRUE) <= 1))
  # strata are roughly equal in size
  sizes <- colSums(multi == 1, na.rm = TRUE)
  expect_lt(max(sizes) - min(sizes), 0.05 * sum(sizes))
})

test_that("categorical stratification creates one stratum per category", {
  s <- quietly(stratify(pheno, strat_cat = strat_cat))
  n_cat <- length(unique(strat_cat[pheno[, 3] == 1, 3]))
  expect_equal(s$strat_details[[1]]$K, n_cat)
  expect_true(all(paste0("categorical_1_", 1:n_cat) %in% colnames(s$multi)))
})

test_that("individuals with a missing phenotype are missing in all strata", {
  p <- pheno; p[1:20, 3] <- NA
  s <- quietly(stratify(p, strat_cont = strat_cont, K = 3))
  strata_cols <- grep("^continuous_1_", colnames(s$multi), value = TRUE)
  expect_true(all(is.na(s$multi[1:20, strata_cols])))
})
