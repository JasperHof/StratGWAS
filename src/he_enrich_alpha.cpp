// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

using namespace Rcpp;

// Profiled version of he_multi_part_enrich: evaluates a grid of alpha values
// controlling the MAF-dependence of per-SNP heritability. For each alpha, each
// standardized genotype column is reweighted by [2f(1-f)]^{(1+alpha)/2}, the
// partitioned HE regression is run exactly as before, and a genome-wide fit
// score is accumulated so the best-fitting alpha can be chosen.
//
// alpha = -1 reproduces he_multi_part_enrich exactly (weight = 1 for all SNPs).

// [[Rcpp::export]]
Rcpp::List he_multi_part_enrich_alpha(const std::string& filename,
                                       const SEXP pheno_mat,
                                       const IntegerMatrix snp_cat,
                                       const CharacterVector cat_names,
                                       int block_size,
                                       const std::string& outfile,
                                       Rcpp::NumericVector alpha_grid =
                                           Rcpp::NumericVector::create(-1.0, -0.75, -0.5, -0.25, 0.0)) {

    int n_cat = cat_names.size();
    int n_alpha = alpha_grid.size();

    // Get dimensions of .bim and .fam file, and read genotype ids
    int n_snps = count_lines(filename + ".bim");
    List fam = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    int n_total_inds = geno_iid.size();

    if (snp_cat.nrow() != n_snps)
        stop("snp_cat must have one row per SNP in the .bim file");
    if (snp_cat.ncol() != n_cat)
        stop("snp_cat must have one column per entry in cat_names");

    Rcpp::Rcout << "Number of individuals (from .fam): " << n_total_inds << std::endl;
    Rcpp::Rcout << "Number of SNPs (from .bim): " << n_snps << std::endl;
    Rcpp::Rcout << "Number of categories: " << n_cat << std::endl;
    Rcpp::Rcout << "Number of alpha values: " << n_alpha << std::endl;

    // Read in phenotype
    Rcpp::NumericMatrix pheno;
    CharacterVector pheno_ids;

    if (Rf_isMatrix(pheno_mat) && !Rf_isNull(rownames(pheno_mat))) {
        pheno = as<NumericMatrix>(pheno_mat);
        pheno_ids = rownames(pheno_mat);
    } else {
        stop("Phenotype must be a numerical matrix with IDs as rownames");
    }

    int n_pheno = pheno.cols();

    CharacterVector pheno_names;
    if (Rf_isNull(colnames(pheno_mat))) {
        pheno_names = CharacterVector(n_pheno);
        for (int i = 0; i < n_pheno; ++i) pheno_names[i] = "V" + std::to_string(i + 1);
    } else {
        pheno_names = colnames(pheno_mat);
    }

    // Match phenotype ids to genotype ids
    IntegerVector match_idx = match(geno_iid, pheno_ids);

    std::vector<int> geno_keep;
    std::vector<int> pheno_keep;

    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            geno_keep.push_back(i);
            pheno_keep.push_back(match_idx[i] - 1);
        }
    }

    if (geno_keep.empty())
        stop("No overlapping individuals between genotype and phenotype");

    NumericMatrix pheno_sub(geno_keep.size(), n_pheno);
    for (size_t i = 0; i < geno_keep.size(); ++i) {
        for (int j = 0; j < n_pheno; ++j) {
            pheno_sub(i, j) = pheno(pheno_keep[i], j);
        }
    }

    pheno = pheno_sub;
    int n_inds = pheno.nrow();

    Rcout << "Number of individuals with complete data: " << n_inds << "\n";

    // Accumulators across all blocks, now indexed by alpha.
    // h2_total[a] is (n_pheno x n_cat) for alpha value a; likewise se2_total[a].
    std::vector<Rcpp::NumericMatrix> h2_total(n_alpha);
    std::vector<Rcpp::NumericMatrix> se2_total(n_alpha);
    for (int a = 0; a < n_alpha; ++a) {
        h2_total[a] = Rcpp::NumericMatrix(n_pheno, n_cat);
        se2_total[a] = Rcpp::NumericMatrix(n_pheno, n_cat);
        colnames(h2_total[a]) = cat_names;
        rownames(h2_total[a]) = pheno_names;
    }

    // Genome-wide h2 and its variance, per alpha x phenotype
    std::vector<std::vector<double>> genome_h2_total(n_alpha, std::vector<double>(n_pheno, 0.0));
    std::vector<std::vector<double>> genome_se2_total(n_alpha, std::vector<double>(n_pheno, 0.0));

    // Genome-wide HE fit score per alpha x phenotype.
    // The HE objective minimizes || y y' - sum_c sigma2_c K_c ||_F^2. For the
    // fitted components the minimized residual reduces (per block) to
    //   ||y||^4 - sol' B_fit,  where B_fit are the fitted RHS moments.
    // We accumulate the per-block reduction in fit so that, summed over the
    // genome, a smaller residual (equivalently a larger explained moment) marks
    // the better-fitting alpha. We store the explained part; larger = better.
    std::vector<std::vector<double>> fit_score(n_alpha, std::vector<double>(n_pheno, 0.0));

    // Output file for block-level results (now carries an alpha column)
    std::ofstream fout(outfile);
    if (!fout.is_open())
        stop("Could not open outfile for writing");
    fout << "alpha\tblock\tblock_start\tblock_end\tphenotype\tcategory\tm_c\th2\tse\tenrichment\n";

    int nr_blocks = (n_snps + block_size - 1) / block_size;
    int nmcmc = 20;

    for (int b = 0; b < nr_blocks; ++b) {

        if (b % 100 == 0) {
            Rcpp::Rcout << "Processing genotype block " << b + 1 << "/" << nr_blocks << "\n";
        }

        int block_start = b * block_size;
        int block_end = std::min(n_snps - 1, (b + 1) * block_size - 1);
        int n_snps_block = block_end - block_start + 1;

        IntegerMatrix geno_block_full = readBedBlock(
            filename + ".bed",
            n_total_inds, n_snps,
            0, n_total_inds - 1,
            block_start, block_end
        );

        IntegerMatrix geno_block(n_inds, n_snps_block);
        for (int i = 0; i < n_inds; ++i) {
            for (int j = 0; j < n_snps_block; ++j) {
                geno_block(i, j) = geno_block_full(geno_keep[i], j);
            }
        }

        Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(geno_block);

        // MAF per column, captured during standardization so we can reweight
        // the standardized columns by alpha later without re-reading genotypes.
        Eigen::VectorXd maf_block(n_snps_block);

        // Standardize columns (alpha = -1 baseline: unit-variance columns)
        for (int col = 0; col < X.cols(); ++col) {
            Eigen::VectorXd v = X.col(col);

            double sum = 0.0;
            int n_valid = 0;
            for (int i = 0; i < v.size(); ++i) {
                if (v[i] != -1) {
                    sum += v[i];
                    n_valid++;
                }
            }

            if (n_valid == 0) {
                X.col(col).setZero();
                maf_block(col) = 0.0;
                continue;
            }

            double mean = sum / n_valid;

            // minor allele frequency for a 0/1/2 genotype column
            double maf = mean / 2.0;
            if (maf > 0.5) maf = 1.0 - maf;
            maf_block(col) = maf;

            double sq_sum = 0.0;
            for (int i = 0; i < v.size(); ++i) {
                if (v[i] == -1) v[i] = mean;
                sq_sum += (v[i] - mean) * (v[i] - mean);
            }

            double stddev = std::sqrt(sq_sum / (v.size() - 1));

            if (stddev > 1e-10) {
                X.col(col) = (v.array() - mean) / stddev;
            } else {
                X.col(col).setZero();
            }
        }

        // Group columns of this block by category (categories may overlap,
        // so a column can be added to more than one list). This grouping is
        // independent of alpha, so we do it once.
        std::vector<std::vector<int>> cat_col_idx(n_cat);
        for (int j = 0; j < n_snps_block; ++j) {
            int snp_row = block_start + j;
            for (int c = 0; c < n_cat; ++c) {
                if (snp_cat(snp_row, c) == 1) {
                    cat_col_idx[c].push_back(j);
                }
            }
        }

        std::vector<int> present_cats;
        for (int c = 0; c < n_cat; ++c) {
            if (!cat_col_idx[c].empty()) present_cats.push_back(c);
        }
        int n_present = present_cats.size();

        if (n_present == 0) continue; // no categorized SNPs in this block

        // Precompute Y once (shared across alpha)
        Eigen::MatrixXd Y(n_inds, n_pheno);
        for (int i = 0; i < n_pheno; ++i) {
            Y.col(i) = Eigen::Map<Eigen::VectorXd>(pheno.column(i).begin(), n_inds);
        }

        // ---- Loop over alpha values -------------------------------------
        for (int a = 0; a < n_alpha; ++a) {
            double alpha = alpha_grid[a];

            // Build per-category submatrices, reweighting each standardized
            // column by w = [2 f (1 - f)]^{(1+alpha)/2}. At alpha = -1 the
            // exponent is 0 so w = 1 and this is identical to the original.
            std::vector<Eigen::MatrixXd> Xc(n_present);
            std::vector<int> m_c(n_present);
            int m_block_total = 0;

            for (int ci = 0; ci < n_present; ++ci) {
                int c = present_cats[ci];
                int mc = cat_col_idx[c].size();
                m_c[ci] = mc;
                m_block_total += mc;

                Eigen::MatrixXd Xsub(n_inds, mc);
                for (int k = 0; k < mc; ++k) {
                    int jcol = cat_col_idx[c][k];
                    double f = maf_block(jcol);
                    double w = 1.0;
                    if (f > 0.0 && f < 1.0) {
                        w = std::pow(2.0 * f * (1.0 - f), (1.0 + alpha) / 2.0);
                    }
                    Xsub.col(k) = X.col(jcol) * w;
                }
                Xc[ci] = Xsub;
            }

            // Trace estimation: tr(K_c K_d) for all category pairs via Hutchinson
            Eigen::MatrixXd tr_KK = Eigen::MatrixXd::Zero(n_present, n_present);

            for (int r = 0; r < nmcmc; ++r) {
                NumericVector z_r = rnorm(n_inds, 0.0, 1.0);
                Eigen::VectorXd z = Eigen::Map<Eigen::VectorXd>(z_r.begin(), n_inds);

                std::vector<Eigen::VectorXd> Kz(n_present);
                for (int ci = 0; ci < n_present; ++ci) {
                    Eigen::VectorXd Xtz = Xc[ci].transpose() * z;
                    Kz[ci] = (Xc[ci] * Xtz) / m_c[ci];
                }

                for (int ci = 0; ci < n_present; ++ci) {
                    for (int cj = ci; cj < n_present; ++cj) {
                        double val = Kz[ci].dot(Kz[cj]) / nmcmc;
                        tr_KK(ci, cj) += val;
                        if (ci != cj) tr_KK(cj, ci) += val;
                    }
                }
            }

            // Exact tr(K_c) per category. Under alpha-weighting a column has
            // squared norm w_v^2 * (N-1), so tr(K_c) = (N-1) * mean(w_v^2),
            // which differs from n_inds for alpha != -1 and differs BETWEEN
            // categories (their MAF distributions differ). Using n_inds here
            // would bias the category-vs-residual coupling in an alpha- and
            // category-dependent way, contaminating exactly the per-category
            // alpha comparison. So compute it exactly from the weighted columns.
            // At alpha = -1 every w_v = 1, so tr_Kc = (N-1) ~ n_inds, recovering
            // the original approximation.
            std::vector<double> tr_Kc(n_present);
            for (int ci = 0; ci < n_present; ++ci) {
                double t = 0.0;
                for (int k = 0; k < m_c[ci]; ++k) {
                    t += Xc[ci].col(k).squaredNorm();
                }
                tr_Kc[ci] = t / m_c[ci];
            }

            // Per-category X'y (for all phenotypes at once)
            std::vector<Eigen::MatrixXd> XtY_c(n_present);
            for (int ci = 0; ci < n_present; ++ci) {
                XtY_c[ci] = Xc[ci].transpose() * Y; // m_c x n_pheno
            }

            // Build the (n_present+1) x (n_present+1) system matrix A (categories
            // + residual). Identical across phenotypes within a block+alpha, so
            // invert once.
            int p = n_present + 1;
            Eigen::MatrixXd A(p, p);
            for (int ci = 0; ci < n_present; ++ci) {
                for (int cj = 0; cj < n_present; ++cj) {
                    A(ci, cj) = tr_KK(ci, cj);
                }
                // off-diagonal coupling to the residual/identity component is
                // tr(K_c * I) = tr(K_c), now exact and category-specific
                A(ci, n_present) = tr_Kc[ci];
                A(n_present, ci) = tr_Kc[ci];
            }
            A(n_present, n_present) = n_inds;

            Eigen::MatrixXd A_inv = A.completeOrthogonalDecomposition().pseudoInverse();

            for (int i = 0; i < n_pheno; ++i) {

                Eigen::VectorXd B(p);
                for (int ci = 0; ci < n_present; ++ci) {
                    // y_i' K_c y_i = ||X_c' y_i||^2 / m_c
                    B(ci) = XtY_c[ci].col(i).squaredNorm() / m_c[ci];
                }
                B(n_present) = Y.col(i).squaredNorm();

                Eigen::VectorXd sol = A_inv * B;

                double h2_block = 0.0;
                double h2_block_var = 0.0;
                for (int ci = 0; ci < n_present; ++ci) {
                    h2_block += sol[ci];
                    for (int cj = 0; cj < n_present; ++cj) {
                        h2_block_var += A_inv(ci, cj);
                    }
                }
                h2_block_var = 2.0 * std::max(0.0, h2_block_var);
                double h2_block_se = std::sqrt(h2_block_var);

                // HE fit contribution for this block: sol' B is the moment
                // explained by the fitted components (larger = better fit). We
                // exclude the residual/identity term (last entry) so the score
                // reflects only variance explained by the genetic kernels.
                double fit_contrib = 0.0;
                for (int ci = 0; ci < n_present; ++ci) {
                    fit_contrib += sol[ci] * B(ci);
                }
                fit_score[a][i] += fit_contrib;

                for (int ci = 0; ci < n_present; ++ci) {
                    int c = present_cats[ci];
                    double h2_c = sol[ci];
                    double var_c = 2.0 * std::max(0.0, A_inv(ci, ci));
                    double se_c = std::sqrt(var_c);

                    h2_total[a](i, c) += h2_c;
                    se2_total[a](i, c) += var_c;

                    double expected_share = (double)m_c[ci] / m_block_total;
                    double h2_share = (h2_block != 0.0) ? h2_c / h2_block : NA_REAL;
                    double enrichment = (h2_block != 0.0) ? h2_share / expected_share : NA_REAL;

                    fout << alpha << "\t" << b << "\t" << block_start << "\t" << block_end << "\t"
                         << pheno_names[i] << "\t" << cat_names[c] << "\t"
                         << m_c[ci] << "\t" << h2_c << "\t" << se_c << "\t" << enrichment << "\n";
                }

                // Block-level total across categories
                fout << alpha << "\t" << b << "\t" << block_start << "\t" << block_end << "\t"
                     << pheno_names[i] << "\t" << "TOTAL" << "\t"
                     << m_block_total << "\t" << h2_block << "\t" << h2_block_se << "\t" << "NA" << "\n";

                genome_h2_total[a][i] += h2_block;
                genome_se2_total[a][i] += h2_block_var;
            }
        } // end alpha loop
    } // end block loop

    // Genome-wide summary rows per alpha (blocks treated as independent)
    for (int a = 0; a < n_alpha; ++a) {
        double alpha = alpha_grid[a];
        for (int i = 0; i < n_pheno; ++i) {
            for (int c = 0; c < n_cat; ++c) {
                fout << alpha << "\tALL\tNA\tNA\t" << pheno_names[i] << "\t" << cat_names[c] << "\t"
                     << "NA\t" << h2_total[a](i, c) << "\t" << std::sqrt(se2_total[a](i, c)) << "\tNA\n";
            }
            fout << alpha << "\tALL\tNA\tNA\t" << pheno_names[i] << "\tTOTAL\tNA\t"
                 << genome_h2_total[a][i] << "\t" << std::sqrt(genome_se2_total[a][i]) << "\tNA\n";
        }
    }

    fout.close();

    // Assemble return objects. For each alpha we return h2 and se matrices; we
    // also return the genome-wide h2, its se, and the fit score, as
    // (n_alpha x n_pheno) matrices so the best alpha can be read off directly.
    Rcpp::List h2_list(n_alpha);
    Rcpp::List se_list(n_alpha);
    for (int a = 0; a < n_alpha; ++a) {
        Rcpp::NumericMatrix se_mat(n_pheno, n_cat);
        colnames(se_mat) = cat_names;
        rownames(se_mat) = pheno_names;
        for (int i = 0; i < n_pheno; ++i)
            for (int c = 0; c < n_cat; ++c)
                se_mat(i, c) = std::sqrt(se2_total[a](i, c));
        h2_list[a] = h2_total[a];
        se_list[a] = se_mat;
    }

    Rcpp::NumericMatrix genome_h2_out(n_alpha, n_pheno);
    Rcpp::NumericMatrix genome_se_out(n_alpha, n_pheno);
    Rcpp::NumericMatrix fit_out(n_alpha, n_pheno);
    rownames(genome_h2_out) = as<CharacterVector>(wrap(alpha_grid));
    rownames(genome_se_out) = as<CharacterVector>(wrap(alpha_grid));
    rownames(fit_out)       = as<CharacterVector>(wrap(alpha_grid));
    colnames(genome_h2_out) = pheno_names;
    colnames(genome_se_out) = pheno_names;
    colnames(fit_out)       = pheno_names;
    for (int a = 0; a < n_alpha; ++a) {
        for (int i = 0; i < n_pheno; ++i) {
            genome_h2_out(a, i) = genome_h2_total[a][i];
            genome_se_out(a, i) = std::sqrt(genome_se2_total[a][i]);
            fit_out(a, i)       = fit_score[a][i];
        }
    }

    // Best alpha per phenotype = the one maximizing the genome-wide fit score
    Rcpp::NumericVector best_alpha(n_pheno);
    for (int i = 0; i < n_pheno; ++i) {
        int best_a = 0;
        double best_val = fit_out(0, i);
        for (int a = 1; a < n_alpha; ++a) {
            if (fit_out(a, i) > best_val) { best_val = fit_out(a, i); best_a = a; }
        }
        best_alpha[i] = alpha_grid[best_a];
    }
    best_alpha.names() = pheno_names;

    return Rcpp::List::create(
        Rcpp::Named("alpha_grid") = alpha_grid,
        Rcpp::Named("h2") = h2_list,
        Rcpp::Named("se") = se_list,
        Rcpp::Named("genome_h2") = genome_h2_out,
        Rcpp::Named("genome_se") = genome_se_out,
        Rcpp::Named("fit_score") = fit_out,
        Rcpp::Named("best_alpha") = best_alpha
    );
}
