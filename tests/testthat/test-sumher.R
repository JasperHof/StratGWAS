test_that("sumher recovers simulated heritability", {
  for (h2 in c(0, 0.1, 0.3)) {
    d <- simulate_sumstats(h2 = h2)
    res <- sumher(d$ss, d$ld)
    expect_true(res$converged)
    expect_lt(abs(res$h2_snp - h2), 0.03)
  }
})

test_that("sumher_cov of a trait with itself equals its heritability (rg = 1)", {
  d <- simulate_sumstats(h2 = 0.3)
  res <- sumher_cov(d$ss, d$ss, d$ld)
  expect_lt(abs(res$h2_AB - 0.3), 0.03)
  expect_equal(res$rg, 1, tolerance = 1e-6)
})

test_that("sumher functions refuse misaligned summary statistics and LD scores", {
  d <- simulate_sumstats(M = 2000)
  expect_error(sumher(d$ss[-1, ], d$ld), "differ in length")
  expect_error(sumher_cov(d$ss, d$ss[-1, ], d$ld), "do not match")
})

test_that("sumher_cov survives a singular design matrix (regression test)", {
  d <- simulate_sumstats(M = 2000)
  d$ss$MAF <- 0.3          # constant MAF and LD make the two regressors collinear
  d$ld[] <- 10
  expect_no_error(suppressWarnings(sumher_cov(d$ss, d$ss, d$ld)))
})

test_that("jackknife returns B estimates and a positive standard error", {
  d <- simulate_sumstats(h2 = 0.3, M = 5000)
  j <- sumher_jack(d$ss, d$ld, B = 20)
  expect_length(j$ests, 20)
  expect_gt(j$se, 0)
  expect_false(anyNA(j$ests))

  jc <- sumher_cov_jack(d$ss, d$ss, d$ld, B = 20)
  expect_length(jc$ests, 20)
  expect_gt(jc$se, 0)
})
