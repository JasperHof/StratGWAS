# Shared fixtures for the StratGWAS tests.
# Files named helper-*.R are sourced automatically by testthat before the tests run.

# Evaluate expr without printing progress output, messages or warnings
quietly <- function(expr) {
  utils::capture.output(res <- suppressWarnings(suppressMessages(expr)))
  res
}

example_bed <- function() {
  sub("\\.bed$", "", system.file("extdata", "data.bed", package = "StratGWAS"))
}

# Small simulated summary statistics that follow the SumHer model,
# E[chisq_j] = 1 + N * h2 * q_j * l_j, with q_j proportional to (MAF(1-MAF))^(1+alpha)
simulate_sumstats <- function(M = 20000, N = 10000, h2 = 0.3, alpha = -0.25, seed = 1) {
  withr_seed <- function(seed, code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    code
  }
  withr_seed(seed, {
    ld  <- runif(M, 1, 60)
    maf <- runif(M, 0.05, 0.5)
    q   <- (maf * (1 - maf))^(1 + alpha)
    q   <- q / sum(q)
    mu  <- 1 + N * h2 * q * ld
    z   <- rnorm(M, sqrt(mu - 1), 1) * sample(c(-1, 1), M, replace = TRUE)
    list(
      ss = data.frame(Predictor = paste0("snp", seq_len(M)), N = N, Chisq = z^2,
                      MAF = maf, BETA = z / sqrt(N), SE = 1 / sqrt(N)),
      ld = ld
    )
  })
}

# Run stratify() + compute_gencov() + transform() once on the example data and
# cache the result, so the (slow) GWAS and LD-score steps only run once per test session.
.stratgwas_cache <- new.env()

example_workflow <- function() {
  if (!is.null(.stratgwas_cache$res)) return(.stratgwas_cache$res)

  RcppParallel::setThreadOptions(numThreads = 2)

  pheno      <- get(utils::data("pheno", package = "StratGWAS", envir = environment()))
  strat_cont <- get(utils::data("strat_cont", package = "StratGWAS", envir = environment()))
  strat_cat  <- get(utils::data("strat_cat", package = "StratGWAS", envir = environment()))

  # Give some individuals a missing phenotype, to check they stay missing downstream
  na_ids <- pheno[1:50, 2]
  pheno[1:50, 3] <- NA

  outfile <- tempfile("stratgwas_test")
  strata <- quietly(stratify(pheno, strat_cont = strat_cont, strat_cat = strat_cat, K = 3))
  gencov <- quietly(compute_gencov(strata, example_bed(), outfile = outfile, jack = TRUE, B = 10))
  trans  <- quietly(transform(strata, gencov, outfile))

  .stratgwas_cache$res <- list(pheno = pheno, na_ids = na_ids, strata = strata,
                               gencov = gencov, trans = trans, outfile = outfile)
  .stratgwas_cache$res
}
