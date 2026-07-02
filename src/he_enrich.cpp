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

// [[Rcpp::export]]
Rcpp::List he_multi_part_enrich(const std::string& filename,
                                 const SEXP pheno_mat,
                                 const IntegerMatrix snp_cat,
                                 const CharacterVector cat_names,
                                 int block_size,
                                 const std::string& outfile) {

    int n_cat = cat_names.size();

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

    // Accumulators across all blocks (h2 and its variance, per phenotype x category)
    Rcpp::NumericMatrix h2_total(n_pheno, n_cat);
    Rcpp::NumericMatrix se2_total(n_pheno, n_cat); // variance, not SE yet
    colnames(h2_total) = cat_names;
    rownames(h2_total) = pheno_names;

    std::vector<double> genome_h2_total(n_pheno, 0.0);
    std::vector<double> genome_se2_total(n_pheno, 0.0);

    // Output file for block-level results
    std::ofstream fout(outfile);
    if (!fout.is_open())
        stop("Could not open outfile for writing");
    fout << "block\tblock_start\tblock_end\tphenotype\tcategory\tm_c\th2\tse\tenrichment\n";

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

        // Standardize columns
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
                continue;
            }

            double mean = sum / n_valid;

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
        // so a column can be added to more than one list)
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

        // Build per-category submatrices
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
                Xsub.col(k) = X.col(cat_col_idx[c][k]);
            }
            Xc[ci] = Xsub;
        }

        // Trace estimation: tr(K_c K_d) for all category pairs via Hutchinson estimator
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

        // tr(K_c) approx n_inds, same approximation as the unpartitioned model
        double tr_K_approx = n_inds;

        // Per-category X'y (for all phenotypes at once)
        Eigen::MatrixXd Y(n_inds, n_pheno);
        for (int i = 0; i < n_pheno; ++i) {
            Y.col(i) = Eigen::Map<Eigen::VectorXd>(pheno.column(i).begin(), n_inds);
        }

        std::vector<Eigen::MatrixXd> XtY_c(n_present);
        for (int ci = 0; ci < n_present; ++ci) {
            XtY_c[ci] = Xc[ci].transpose() * Y; // m_c x n_pheno
        }

        // Build the (n_present+1) x (n_present+1) system matrix T (categories + residual).
        // This is identical across phenotypes within a block, so invert it once.
        int p = n_present + 1;
        Eigen::MatrixXd A(p, p);
        for (int ci = 0; ci < n_present; ++ci) {
            for (int cj = 0; cj < n_present; ++cj) {
                A(ci, cj) = tr_KK(ci, cj);
            }
            A(ci, n_present) = tr_K_approx;
            A(n_present, ci) = tr_K_approx;
        }
        A(n_present, n_present) = n_inds;

        // A_inv is reused both for the point estimates (A_inv * B) and for the
        // SE formula Cov[sigma^2] = 2 * A_inv, per the Method-of-Moments result
        // Cov[q]_ij = 2*tr(Ki Kj) = 2*A_ij, so Cov[sigma^2] = A^-1 (2A) A^-1 = 2*A^-1.
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

            for (int ci = 0; ci < n_present; ++ci) {
                int c = present_cats[ci];
                double h2_c = sol[ci];
                double var_c = 2.0 * std::max(0.0, A_inv(ci, ci));
                double se_c = std::sqrt(var_c);

                h2_total(i, c) += h2_c;
                se2_total(i, c) += var_c;

                double expected_share = (double)m_c[ci] / m_block_total;
                double h2_share = (h2_block != 0.0) ? h2_c / h2_block : NA_REAL;
                double enrichment = (h2_block != 0.0) ? h2_share / expected_share : NA_REAL;

                fout << b << "\t" << block_start << "\t" << block_end << "\t"
                     << pheno_names[i] << "\t" << cat_names[c] << "\t"
                     << m_c[ci] << "\t" << h2_c << "\t" << se_c << "\t" << enrichment << "\n";
            }

            // Block-level total across categories (categories may overlap, so this
            // is the sum of partitioned components, not a disjoint total)
            fout << b << "\t" << block_start << "\t" << block_end << "\t"
                 << pheno_names[i] << "\t" << "TOTAL" << "\t"
                 << m_block_total << "\t" << h2_block << "\t" << h2_block_se << "\t" << "NA" << "\n";

            genome_h2_total[i] += h2_block;
            genome_se2_total[i] += h2_block_var;
        }
    }

    // Genome-wide summary rows (blocks treated as independent, so variances add)
    for (int i = 0; i < n_pheno; ++i) {
        for (int c = 0; c < n_cat; ++c) {
            fout << "ALL\tNA\tNA\t" << pheno_names[i] << "\t" << cat_names[c] << "\t"
                 << "NA\t" << h2_total(i, c) << "\t" << std::sqrt(se2_total(i, c)) << "\tNA\n";
        }
        fout << "ALL\tNA\tNA\t" << pheno_names[i] << "\tTOTAL\tNA\t"
             << genome_h2_total[i] << "\t" << std::sqrt(genome_se2_total[i]) << "\tNA\n";
    }

    fout.close();

    Rcpp::NumericMatrix se_total(n_pheno, n_cat);
    colnames(se_total) = cat_names;
    rownames(se_total) = pheno_names;
    for (int i = 0; i < n_pheno; ++i)
        for (int c = 0; c < n_cat; ++c)
            se_total(i, c) = std::sqrt(se2_total(i, c));

    Rcpp::NumericVector genome_h2_out(genome_h2_total.begin(), genome_h2_total.end());
    Rcpp::NumericVector genome_se_out(n_pheno);
    for (int i = 0; i < n_pheno; ++i) genome_se_out[i] = std::sqrt(genome_se2_total[i]);
    genome_h2_out.names() = pheno_names;
    genome_se_out.names() = pheno_names;

    return Rcpp::List::create(
        Rcpp::Named("h2") = h2_total,
        Rcpp::Named("se") = se_total,
        Rcpp::Named("genome_h2") = genome_h2_out,
        Rcpp::Named("genome_se") = genome_se_out
    );
}
