// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────
// Reads a burden file written by compute_gene_burden(): no header, one row
// per gene, tab-separated allele counts per individual. Returns an
// (n_genes x n_inds) matrix.
//
// Assumes rows are in genomic order (true for compute_gene_burden's output,
// since genes are processed in .bim-sorted order) — this is what makes the
// "adjacent rows" = "adjacent genomically" assumption valid downstream.
// ─────────────────────────────────────────────────────────────────────────
static MatrixXd read_burden_file(const std::string& path, int n_inds_hint = -1) {
    std::ifstream in(path);
    if (!in.is_open()) Rcpp::stop("Could not open burden file: " + path);

    std::vector<std::vector<double>> rows;
    std::string line;
    int n_inds = -1;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row;
        if (n_inds_hint > 0) row.reserve(n_inds_hint);
        std::string tok;
        while (std::getline(iss, tok, '\t')) {
            row.push_back(std::stod(tok));
        }
        if (n_inds == -1) n_inds = (int)row.size();
        else if ((int)row.size() != n_inds)
            Rcpp::stop("Inconsistent number of columns in burden file at row " +
                       std::to_string(rows.size() + 1));
        rows.push_back(std::move(row));
    }
    in.close();

    int n_genes = (int)rows.size();
    if (n_genes == 0) Rcpp::stop("Burden file is empty: " + path);

    MatrixXd B(n_genes, n_inds);
    for (int g = 0; g < n_genes; ++g)
        for (int i = 0; i < n_inds; ++i)
            B(g, i) = rows[g][i];

    return B;
}

// ─────────────────────────────────────────────────────────────────────────
// Standardizes each row (gene/burden) to mean 0, sd 1 across individuals.
// Burdens with zero variance (e.g. monomorphic / all-zero burden) are left
// as exact zeros after centering and flagged, so they contribute r = 0 to
// every pairwise correlation rather than producing NaN/Inf.
// ─────────────────────────────────────────────────────────────────────────
static void standardize_rows(MatrixXd& B, std::vector<bool>& is_zero_var) {
    int n_genes = B.rows();
    int n_inds  = B.cols();
    is_zero_var.assign(n_genes, false);

    for (int g = 0; g < n_genes; ++g) {
        double mean = B.row(g).mean();
        VectorXd centered = B.row(g).array() - mean;
        double ss = centered.squaredNorm();
        if (ss < 1e-12) {
            B.row(g).setZero();
            is_zero_var[g] = true;
        } else {
            double sd = std::sqrt(ss / (n_inds - 1));
            B.row(g) = centered / sd;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// compute_burden_weights
//
// Computes pairwise correlations between each burden (gene) and its
// genomically-adjacent neighbours (within a local window of at most
// max_neighbors genes on EACH side, so up to 2*max_neighbors total
// comparisons per burden), then derives LDAK-style weights:
//
//     w_j = 1 / (1 + sum_{k in window, k != j} r_jk^2)
//
// Burdens with many highly-correlated neighbours (redundant signal, e.g.
// overlapping window sizes or genes in high LD) get down-weighted; burdens
// with little local correlation keep weight close to 1.
//
// Arguments:
//   burden_file   - path to a burden file from compute_gene_burden()
//                    (no header, tab-separated, genes in genomic row order)
//   max_neighbors - max number of adjacent genes to compare against on
//                    EACH side (so total window size <= 2*max_neighbors + 1)
//   n_inds_hint   - optional hint for row-reading buffer sizing (cosmetic)
//
// Returns a List with:
//   weights    - NumericVector, one weight per gene, same order as input rows
//   m_eff      - sum(weights): the "effective number of independent burdens"
//   mean_r2    - mean squared correlation across all computed neighbour pairs
//                (diagnostic: how much redundancy is in the data overall)
//   pair_chr_idx_i / pair_chr_idx_j / pair_r (optional diagnostics, see below)
// [[Rcpp::export]]
List compute_burden_weights(
    const std::string& burden_file,
    int max_neighbors = 100,
    int n_inds_hint = -1,
    bool return_pairs = false   // if true, also returns the full pairwise list (can be large)
) {
    if (max_neighbors < 1) Rcpp::stop("max_neighbors must be >= 1");

    MatrixXd B = read_burden_file(burden_file, n_inds_hint);
    int n_genes = B.rows();
    int n_inds  = B.cols();

    Rcout << "Read " << n_genes << " burdens x " << n_inds << " individuals\n";

    std::vector<bool> is_zero_var;
    standardize_rows(B, is_zero_var);   // rows now standardized -> dot product = correlation

    int n_zero_var = 0;
    for (bool z : is_zero_var) if (z) n_zero_var++;
    if (n_zero_var > 0)
        Rcout << "Warning: " << n_zero_var
              << " burdens have zero variance and are treated as uncorrelated (r=0) with everything\n";

    VectorXd sum_r2 = VectorXd::Zero(n_genes);   // sum_{k in window} r_jk^2, per gene j
    double total_r2_sum = 0.0;
    long   total_pairs   = 0;

    // Diagnostic pairwise output (optional, can get large for genome-wide runs)
    std::vector<int>    pair_i, pair_j;
    std::vector<double> pair_r;

    int denom = n_inds - 1;   // standardized rows -> correlation = dot/(n-1)

    for (int j = 0; j < n_genes; ++j) {
        int lo = std::max(0, j - max_neighbors);
        int hi = std::min(n_genes - 1, j + max_neighbors);

        for (int k = lo; k <= hi; ++k) {
            if (k == j) continue;
            // Only compute each pair once (k > j), but accumulate into BOTH
            // j's and k's running sum_r2, since correlation is symmetric.
            if (k < j) continue;

            double r = B.row(j).dot(B.row(k)) / (double)denom;
            // Numerical safety: standardized rows can give |r| fractionally > 1
            if (r > 1.0) r = 1.0;
            if (r < -1.0) r = -1.0;

            double r2 = r * r;
            sum_r2(j) += r2;
            sum_r2(k) += r2;
            total_r2_sum += r2;
            total_pairs++;

            if (return_pairs) {
                pair_i.push_back(j);
                pair_j.push_back(k);
                pair_r.push_back(r);
            }
        }

        if ((j + 1) % 5000 == 0 || j + 1 == n_genes)
            Rcout << "Processed correlations for " << j + 1 << "/" << n_genes << " burdens\n";
    }

    // ── Derive weights ──────────────────────────────────────────────────
    NumericVector weights(n_genes);
    for (int j = 0; j < n_genes; ++j) {
        weights[j] = 1.0 / (1.0 + sum_r2(j));
    }

    double m_eff = 0.0;
    for (int j = 0; j < n_genes; ++j) m_eff += weights[j];

    double mean_r2 = (total_pairs > 0) ? (total_r2_sum / (double)total_pairs) : 0.0;

    Rcout << "Done. m_raw = " << n_genes
          << ", m_eff = " << m_eff
          << ", mean pairwise r^2 (within window) = " << mean_r2 << "\n";

    List result = List::create(
        Named("weights")  = weights,
        Named("m_raw")    = n_genes,
        Named("m_eff")    = m_eff,
        Named("mean_r2")  = mean_r2,
        Named("n_zero_var_burdens") = n_zero_var
    );

    if (return_pairs) {
        result["pair_i"] = wrap(pair_i);   // 0-indexed row of first burden in pair
        result["pair_j"] = wrap(pair_j);   // 0-indexed row of second burden in pair
        result["pair_r"] = wrap(pair_r);   // correlation between them
    }

    return result;
}