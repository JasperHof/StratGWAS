// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────
// Lightweight index of byte offsets for each line (gene/burden row) in the
// burden file, built with a single fast pass over the file. This is what
// lets us "readBedBlock"-style seek directly to a block of rows instead of
// loading the whole file — the text-file analogue of random access into a
// .bed file.
// ─────────────────────────────────────────────────────────────────────────
struct BurdenFileIndex {
    std::vector<std::streampos> line_offsets;  // byte offset of the start of each row
    int n_genes;
    int n_inds;

    BurdenFileIndex(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in.is_open()) Rcpp::stop("Could not open burden file: " + path);

        n_inds = -1;
        std::string line;
        std::streampos pos = in.tellg();

        while (true) {
            pos = in.tellg();
            if (!std::getline(in, line)) break;
            if (line.empty()) continue;

            if (n_inds == -1) {
                // Count tab-separated fields on the first non-empty line
                int tabs = 0;
                for (char c : line) if (c == '\t') tabs++;
                n_inds = tabs + 1;
            }
            line_offsets.push_back(pos);
        }
        n_genes = (int)line_offsets.size();
        in.close();

        if (n_genes == 0) Rcpp::stop("Burden file is empty: " + path);
        Rcout << "Indexed burden file: " << n_genes << " genes x " << n_inds << " individuals\n";
    }
};

// ─────────────────────────────────────────────────────────────────────────
// Reads rows [row_start, row_end] (inclusive, 0-indexed) from the burden
// file into an Eigen matrix, parsing directly into pre-sized storage (no
// vector<vector<double>> intermediate — fixes the doubling issue from the
// whole-file version).
//
// Uses strtod for parsing, which is meaningfully faster than std::stod
// (no exceptions, no std::string construction per token) — relevant given
// how many numeric conversions a 40GB file requires.
// ─────────────────────────────────────────────────────────────────────────
static MatrixXd read_burden_block(
    const std::string& path,
    const BurdenFileIndex& idx,
    int row_start,
    int row_end
) {
    int n_block = row_end - row_start + 1;
    int n_inds  = idx.n_inds;

    MatrixXd block(n_block, n_inds);

    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) Rcpp::stop("Could not reopen burden file: " + path);

    std::string line;
    for (int r = 0; r < n_block; ++r) {
        in.seekg(idx.line_offsets[row_start + r]);
        std::getline(in, line);

        const char* p = line.c_str();
        char* end;
        for (int c = 0; c < n_inds; ++c) {
            block(r, c) = std::strtod(p, &end);
            p = end;
            if (*p == '\t') ++p;
        }
    }
    in.close();
    return block;
}

// ─────────────────────────────────────────────────────────────────────────
// Standardizes each row to mean 0, sd 1. Zero-variance rows are zeroed out
// (contribute r = 0 to all pairs) rather than producing NaN.
// ─────────────────────────────────────────────────────────────────────────
static void standardize_rows_inplace(MatrixXd& B, std::vector<bool>& is_zero_var) {
    int n_rows = B.rows();
    int n_inds = B.cols();
    is_zero_var.assign(n_rows, false);

    for (int r = 0; r < n_rows; ++r) {
        double mean = B.row(r).mean();
        VectorXd centered = B.row(r).array() - mean;
        double ss = centered.squaredNorm();
        if (ss < 1e-12) {
            B.row(r).setZero();
            is_zero_var[r] = true;
        } else {
            double sd = std::sqrt(ss / (n_inds - 1));
            B.row(r) = centered / sd;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// compute_burden_weights_blockwise
//
// Sliding-window version of compute_burden_weights(): mirrors
// computeLDscoresFromBED()'s block-read strategy. Reads at most
// (2 * max_neighbors) burdens into memory at a time (a "padded" block),
// computes correlations only for the central "core" rows of that block
// against the full padded neighbourhood, then slides forward by
// max_neighbors — so peak memory is O(max_neighbors * n_inds) rather than
// O(n_genes * n_inds), exactly analogous to how the LD-score code keeps
// only ~2 * max_block SNPs in memory rather than all n_snp.
//
// Arguments:
//   burden_file   - path to a burden file from compute_gene_burden()
//   max_neighbors - number of adjacent genes considered on EACH side when
//                    computing a burden's correlation sum (so each burden's
//                    weight reflects up to 2*max_neighbors comparisons)
//
// Returns a List with weights, m_eff, mean_r2 (same fields as the
// whole-file version), plus block-size diagnostics.
// [[Rcpp::export]]
List compute_burden_weights_blockwise(
    const std::string& burden_file,
    int max_neighbors = 100
) {
    if (max_neighbors < 1) Rcpp::stop("max_neighbors must be >= 1");

    BurdenFileIndex idx(burden_file);
    int n_genes = idx.n_genes;
    int n_inds  = idx.n_inds;

    int read_block_size = 2 * max_neighbors;   // mirrors `2 * max_block` in the LD-score code

    VectorXd sum_r2 = VectorXd::Zero(n_genes);
    std::vector<bool> is_zero_var_global(n_genes, false);
    double total_r2_sum = 0.0;
    long   total_pairs   = 0;

    int start = 0;
    int last_printed = -10000;

    while (start < n_genes) {

        if (start / 5000 > last_printed / 5000) {
            Rcout << "Processing gene " << (start / 5000) * 5000 << " / " << n_genes << "\n";
            last_printed = (start / 5000) * 5000;
        }

        int end = std::min(start + read_block_size - 1, n_genes - 1);
        int n_block = end - start + 1;

        // ── Determine which rows within this block are the "core" rows we
        //    actually compute final sum_r2 for (avoids double-counting,
        //    exactly mirroring compute_start/compute_end in the LD-score code) ──
        int core_start, core_end;   // indices relative to block start (0-based)

        if (start == 0 && end == n_genes - 1) {
            core_start = 0;
            core_end   = n_block - 1;
        } else if (start == 0) {
            core_start = 0;
            core_end   = std::min(start + (3 * max_neighbors) / 2 - 1, n_genes - 1) - start;
        } else if (end == n_genes - 1) {
            core_start = max_neighbors / 2;
            core_end   = end - start;
        } else {
            core_start = max_neighbors / 2;
            core_end   = std::min((3 * max_neighbors) / 2 - 1, end - start);
        }

        // ── Read and standardize this block only (not the whole file) ──────
        MatrixXd B = read_burden_block(burden_file, idx, start, end);
        std::vector<bool> is_zero_var_block;
        standardize_rows_inplace(B, is_zero_var_block);
        for (int r = 0; r < n_block; ++r)
            is_zero_var_global[start + r] = is_zero_var_block[r];

        int denom = n_inds - 1;

        // ── For each core row, correlate against all OTHER rows in the
        //    padded block (its full local neighbourhood), same as the
        //    LD-score code's inner loop over k in [0, n_block) ──────────────
        for (int j = core_start; j <= core_end; ++j) {
            int global_j = start + j;

            for (int k = 0; k < n_block; ++k) {
                if (k == j) continue;
                int global_k = start + k;

                // Row-distance window: only count neighbours within
                // max_neighbors rows (mirrors the 1Mb bp check in the
                // LD-score code, but using row/gene index distance since
                // burden files don't carry bp position).
                if (std::abs(global_j - global_k) > max_neighbors) continue;

                double r = B.row(j).dot(B.row(k)) / (double)denom;
                if (r > 1.0) r = 1.0;
                if (r < -1.0) r = -1.0;

                double r2 = r * r;
                sum_r2(global_j) += r2;

                // Only tally each unordered pair once for the diagnostic mean
                if (global_k > global_j) {
                    total_r2_sum += r2;
                    total_pairs++;
                }
            }
        }

        start += max_neighbors;   // slide forward, same as `start_snp += max_block`
    }

    // ── Derive weights ──────────────────────────────────────────────────
    NumericVector weights(n_genes);
    for (int j = 0; j < n_genes; ++j) {
        weights[j] = 1.0 / (1.0 + sum_r2(j));
    }

    double m_eff = 0.0;
    for (int j = 0; j < n_genes; ++j) m_eff += weights[j];

    double mean_r2 = (total_pairs > 0) ? (total_r2_sum / (double)total_pairs) : 0.0;

    int n_zero_var = 0;
    for (bool z : is_zero_var_global) if (z) n_zero_var++;

    Rcout << "Done. m_raw = " << n_genes
          << ", m_eff = " << m_eff
          << ", mean pairwise r^2 (within window) = " << mean_r2
          << ", block size = " << read_block_size << " genes\n";

    return List::create(
        Named("weights")  = weights,
        Named("m_raw")    = n_genes,
        Named("m_eff")    = m_eff,
        Named("mean_r2")  = mean_r2,
        Named("n_zero_var_burdens") = n_zero_var,
        Named("block_size") = read_block_size
    );
}