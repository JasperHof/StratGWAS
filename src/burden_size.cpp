// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;
using Eigen::VectorXd;

// ---------------------------------------------------------------------------
// One window: a contiguous range of SNP indices in the .bim file
// ---------------------------------------------------------------------------
struct Window {
    int chr;
    int start_idx;    // 0-based, inclusive, index into .bim
    int end_idx;      // 0-based, exclusive (half-open)
    long bp_start;
    long bp_end;
};

// ---------------------------------------------------------------------------
// Read .bim into parallel vectors (chr, pos per SNP)
// ---------------------------------------------------------------------------
struct BimData {
    std::vector<int>  chr;
    std::vector<long> pos;
    int n_snps;
};

BimData read_bim_data(const IntegerVector& chr_in, const NumericVector& pos_in) {
    BimData b;
    b.n_snps = chr_in.size();
    b.chr.resize(b.n_snps);
    b.pos.resize(b.n_snps);
    for (int i = 0; i < b.n_snps; ++i) {
        b.chr[i] = chr_in[i];
        b.pos[i] = (long)pos_in[i];
    }
    return b;
}

// ---------------------------------------------------------------------------
// Read a two-column effects file (SNP, effect size) and align it to .bim
// order. SNPs not found in the effects file get effect size 0.
// ---------------------------------------------------------------------------
std::vector<double> read_effects_file(const std::string& effects_file,
                                       const CharacterVector& snp_ids) {
    int n = snp_ids.size();
    std::vector<double> beta(n, 0.0);

    std::unordered_map<std::string, double> eff_map;
    std::ifstream in(effects_file);
    if (!in.is_open()) Rcpp::stop("Could not open effects file: " + effects_file);

    std::string snp;
    double b;
    while (in >> snp >> b) eff_map[snp] = b;
    in.close();

    int n_matched = 0;
    for (int i = 0; i < n; ++i) {
        std::string s = as<std::string>(snp_ids[i]);
        auto it = eff_map.find(s);
        if (it != eff_map.end()) {
            beta[i] = it->second;
            ++n_matched;
        }
    }
    Rcout << "Effects file: matched " << n_matched << "/" << n << " SNPs\n";
    return beta;
}

// ---------------------------------------------------------------------------
// Squared Pearson correlation between two vectors (e.g. burden score vs
// true genetic score), used to convert true regional heritability into
// the heritability actually captured by the burden score.
// ---------------------------------------------------------------------------
double compute_r2(const VectorXd& a, const VectorXd& b) {
    double mean_a = a.mean();
    double mean_b = b.mean();
    VectorXd ac = a.array() - mean_a;
    VectorXd bc = b.array() - mean_b;

    double denom = ac.norm() * bc.norm();
    if (denom < 1e-12) return 0.0;   // one side has no variance (e.g. window has no effects)

    double r = ac.dot(bc) / denom;
    return r * r;
}

// ---------------------------------------------------------------------------
// Window-generation strategies
// Each returns a vector of Window structs covering all SNPs in the bim,
// leaving no SNPs unassigned (last window may be smaller than the target).
// ---------------------------------------------------------------------------

// Strategy 1: fixed basepair window size
std::vector<Window> windows_by_kb(const BimData& bim, double kb_size) {
    long bp_size = (long)(kb_size * 1000.0);
    std::vector<Window> windows;

    int n = bim.n_snps;
    int i = 0;
    while (i < n) {
        int this_chr  = bim.chr[i];
        long tile_start = (bim.pos[i] / bp_size) * bp_size;
        long tile_end   = tile_start + bp_size - 1;

        Window w;
        w.chr       = this_chr;
        w.start_idx = i;
        w.bp_start  = tile_start;
        w.bp_end    = tile_end;

        // Advance while still on same chr and within this bp tile
        while (i < n && bim.chr[i] == this_chr && bim.pos[i] <= tile_end) ++i;

        w.end_idx = i;   // half-open
        windows.push_back(w);
    }
    return windows;
}

// Strategy 2: fixed number of SNPs per window
std::vector<Window> windows_by_nsnps(const BimData& bim, int n_snps_per_window) {
    std::vector<Window> windows;
    int n = bim.n_snps;
    int i = 0;
    while (i < n) {
        Window w;
        w.chr       = bim.chr[i];
        w.start_idx = i;
        w.bp_start  = bim.pos[i];

        // Advance n_snps_per_window steps, but stop at chromosome boundary
        int count = 0;
        while (i < n && bim.chr[i] == w.chr && count < n_snps_per_window) {
            ++i;
            ++count;
        }

        w.end_idx = i;
        w.bp_end  = bim.pos[i - 1];
        windows.push_back(w);
    }
    return windows;
}

// Strategy 3: target variation level, measured as MAC/n_inds per window.
// Accumulates SNPs until total MAC / n_inds reaches the target, then closes
// the window. This requires a first pass over the .bed to compute per-SNP
// MAC, which is done here via the .frq file if available, or estimated from
// the bim MAF otherwise. Since we have a .frq file in your data
// (all_rare.frq), we read MAC = round(2 * n_inds * MAF) from there.
// If no .frq file is found, falls back to computing MAC from genotype data
// directly (slower but always correct).
std::vector<Window> windows_by_variation(
    const BimData& bim,
    double target_mac_per_ind,
    int n_inds,
    const std::string& bed_prefix,
    int chunk_size = 5000
) {
    std::vector<Window> windows;
    int n = bim.n_snps;

    // Try to read MAC from .frq file first (fast path)
    // .frq format: CHR SNP A1 A2 MAF NCHROBS
    std::vector<double> mac(n, -1.0);
    std::string frq_path = bed_prefix + ".frq";
    std::ifstream frq_in(frq_path);
    bool has_frq = frq_in.is_open();

    if (has_frq) {
        Rcout << "Reading MAF from .frq file for variation windows\n";
        std::string line;
        std::getline(frq_in, line);   // skip header
        int idx = 0;
        while (std::getline(frq_in, line) && idx < n) {
            std::istringstream iss(line);
            std::string chr_str, snp, a1, a2;
            double maf;
            int nchrobs;
            if (iss >> chr_str >> snp >> a1 >> a2 >> maf >> nchrobs) {
                mac[idx] = maf * nchrobs;   // nchrobs = 2*n_inds for diploid
            }
            ++idx;
        }
        frq_in.close();
    } else {
        // Slow path: compute MAC from genotype chunks
        Rcout << "No .frq file found - computing MAC from genotype data\n";
        int n_chunks = (n + chunk_size - 1) / chunk_size;
        for (int c = 0; c < n_chunks; ++c) {
            int snp_start = c * chunk_size;
            int snp_end   = std::min(n - 1, snp_start + chunk_size - 1);
            int n_blk     = snp_end - snp_start + 1;

            IntegerMatrix geno_block = readBedBlock(
                bed_prefix + ".bed", n_inds, n,
                0, n_inds - 1, snp_start, snp_end);

            for (int s = 0; s < n_blk; ++s) {
                double snp_mac = 0.0;
                for (int i = 0; i < n_inds; ++i) {
                    int g = geno_block(i, s);
                    if (g != -1) snp_mac += g;
                }
                // Use minor allele: if mac > n_inds, flip to other allele
                if (snp_mac > n_inds) snp_mac = 2.0 * n_inds - snp_mac;
                mac[snp_start + s] = snp_mac;
            }
            if ((c + 1) % 100 == 0)
                Rcout << "MAC scan: chunk " << c+1 << "/" << n_chunks << "\n";
        }
    }

    // Now build windows by accumulating MAC until target is reached
    double target_mac = target_mac_per_ind * n_inds;
    int i = 0;
    while (i < n) {
        Window w;
        w.chr       = bim.chr[i];
        w.start_idx = i;
        w.bp_start  = bim.pos[i];

        double accum_mac = 0.0;
        while (i < n && bim.chr[i] == w.chr) {
            accum_mac += std::max(0.0, mac[i]);
            ++i;
            if (accum_mac >= target_mac) break;
        }

        w.end_idx = i;
        w.bp_end  = bim.pos[i - 1];
        windows.push_back(w);
    }
    return windows;
}

// ---------------------------------------------------------------------------
// Shared output routine: given a set of windows, read genotypes and write
// burden scores. Identical logic for all three strategies. If beta is
// non-empty, also computes true per-window heritability from known effects
// (assumes var(phenotype) = 1) and writes it to out_file + ".hers".
// ---------------------------------------------------------------------------
void write_burdens(
    const std::vector<Window>& windows,
    const std::string& bed_prefix,
    int n_inds,
    int n_snps,
    const std::string& out_file,
    int write_buffer_size,
    const std::vector<double>& beta,   // empty = no heritability calculation
    bool write_burden_scores,          // false = skip writing burden score / kept_windows files
    const std::vector<int>& annot,     // flat n_snps x n_cat, column-major; empty = no categories
    const std::vector<std::string>& cat_names   // one name per annotation column
) {
    bool has_effects = !beta.empty();

    std::ofstream out;
    std::ofstream kept_out;
    std::string kept_file = out_file + ".kept_windows";
    if (write_burden_scores) {
        out.open(out_file);
        if (!out.is_open()) Rcpp::stop("Could not open output file: " + out_file);

        kept_out.open(kept_file);
        if (!kept_out.is_open()) Rcpp::stop("Could not open kept-windows file: " + kept_file);
    }

    std::ofstream hers_out;
    if (has_effects) {
        std::string hers_file = out_file + ".hers";
        hers_out.open(hers_file);
        if (!hers_out.is_open()) Rcpp::stop("Could not open hers file: " + hers_file);
        hers_out << "window\tchr\tbp_start\tbp_end\tcategory\tn_snps\th2\tr2_burden\th2_burden\n";
    }

    // Running genetic score for the whole genome (only tracked if we have
    // effects), used to report the actual total heritability alongside the
    // naive sum of per-window heritabilities
    VectorXd g_total = VectorXd::Zero(n_inds);
    VectorXd burden_total = VectorXd::Zero(n_inds);
    double sum_window_h2 = 0.0;
    double sum_window_h2_burden = 0.0;

    // Same quantities, tracked separately for each annotation category
    const int n_cat = (int) cat_names.size();
    std::vector<VectorXd> g_total_c(n_cat, VectorXd::Zero(n_inds));
    std::vector<VectorXd> burden_total_c(n_cat, VectorXd::Zero(n_inds));
    std::vector<double>   sum_h2_c(n_cat, 0.0), sum_h2b_c(n_cat, 0.0);

    std::vector<std::string> row_buf;
    std::vector<std::string> win_buf;
    row_buf.reserve(write_buffer_size);
    win_buf.reserve(write_buffer_size);

    int n_windows    = (int)windows.size();
    int n_written    = 0;
    int n_skipped    = 0;

    for (int wi = 0; wi < n_windows; ++wi) {
        const Window& w = windows[wi];
        int n_snps_win = w.end_idx - w.start_idx;

        if (n_snps_win == 0) { n_skipped++; continue; }

        IntegerMatrix geno_block = readBedBlock(
            bed_prefix + ".bed", n_inds, n_snps,
            0, n_inds - 1,
            w.start_idx, w.end_idx - 1   // inclusive end for readBedBlock
        );

        VectorXd burden_vec = VectorXd::Zero(n_inds);
        VectorXd g_win = VectorXd::Zero(n_inds);   // true genetic score for this window

        // Per-category versions of the same two vectors, plus the SNP count
        std::vector<VectorXd> burden_c(n_cat, VectorXd::Zero(n_inds));
        std::vector<VectorXd> g_c(n_cat, VectorXd::Zero(n_inds));
        std::vector<int> m_c(n_cat, 0);
        std::vector<int> mycats;                   // categories of the current SNP

        for (int s = 0; s < n_snps_win; ++s) {
            int gidx = w.start_idx + s;
            double b = has_effects ? beta[gidx] : 0.0;

            // Resolve category membership once per SNP, not once per individual.
            // Only needed when we have effects; otherwise mycats stays empty and
            // the per-individual loop below is exactly as cheap as before.
            mycats.clear();
            if (has_effects)
                for (int c = 0; c < n_cat; ++c)
                    if (annot[(size_t) c * n_snps + gidx]) { mycats.push_back(c); m_c[c]++; }

            const bool use_b = has_effects && b != 0.0;
            for (int i = 0; i < n_inds; ++i) {
                int g = geno_block(i, s);
                if (g == -1) continue;
                burden_vec(i) += g;
                double bg = use_b ? b * g : 0.0;
                if (use_b) g_win(i) += bg;
                for (size_t k = 0; k < mycats.size(); ++k) {
                    burden_c[mycats[k]](i) += g;
                    if (use_b) g_c[mycats[k]](i) += bg;
                }
            }
        }

        // Window label: chr_bpstart_bpend_nsnps
        std::ostringstream wlabel;
        wlabel << w.chr << "_" << w.bp_start << "_" << w.bp_end
               << "_" << n_snps_win;

        if (write_burden_scores) {
            // Row: tab-separated burden scores, no labels
            std::ostringstream row;
            for (int i = 0; i < n_inds; ++i) {
                if (i > 0) row << "\t";
                row << burden_vec(i);
            }
            row_buf.push_back(row.str());
            win_buf.push_back(wlabel.str());
        }
        n_written++;

        if (has_effects) {
            double mean_g = g_win.mean();
            double h2_win = (g_win.array() - mean_g).square().sum() / (n_inds - 1);

            // How much of this window's true genetic signal the (unweighted)
            // burden score actually captures: h2_burden = r2 * h2_true
            double r2_win = compute_r2(burden_vec, g_win);
            double h2_win_burden = r2_win * h2_win;

            hers_out << wlabel.str() << "\t" << w.chr << "\t" << w.bp_start << "\t"
                      << w.bp_end << "\tALL\t" << n_snps_win << "\t" << h2_win << "\t"
                      << r2_win << "\t" << h2_win_burden << "\n";

            sum_window_h2 += h2_win;
            sum_window_h2_burden += h2_win_burden;
            g_total += g_win;
            burden_total += burden_vec;

            // One row per category present in this window (empty ones are skipped,
            // matching he_reg_spa, so downstream joins line up)
            for (int c = 0; c < n_cat; ++c) {
                if (m_c[c] == 0) continue;
                double mean_gc = g_c[c].mean();
                double h2_c    = (g_c[c].array() - mean_gc).square().sum() / (n_inds - 1);
                double r2_c    = compute_r2(burden_c[c], g_c[c]);
                double h2b_c   = r2_c * h2_c;

                hers_out << wlabel.str() << "\t" << w.chr << "\t" << w.bp_start << "\t"
                          << w.bp_end << "\t" << cat_names[c] << "\t" << m_c[c] << "\t"
                          << h2_c << "\t" << r2_c << "\t" << h2b_c << "\n";

                sum_h2_c[c]  += h2_c;
                sum_h2b_c[c] += h2b_c;
                g_total_c[c]      += g_c[c];
                burden_total_c[c] += burden_c[c];
            }
        }

        if (write_burden_scores && (int)row_buf.size() >= write_buffer_size) {
            for (const auto& r : row_buf) out << r << "\n";
            for (const auto& wl : win_buf) kept_out << wl << "\n";
            row_buf.clear();
            win_buf.clear();
        }

        if ((wi + 1) % 1000 == 0 || wi + 1 == n_windows)
            Rcout << "Written " << n_written << "/" << n_windows << " windows\n";
    }

    // Final flush
    if (write_burden_scores) {
        for (const auto& r : row_buf) out << r << "\n";
        for (const auto& wl : win_buf) kept_out << wl << "\n";
        out.close();
        kept_out.close();
    }

    if (has_effects) {
        double mean_total = g_total.mean();
        double h2_total_actual = (g_total.array() - mean_total).square().sum() / (n_inds - 1);

        // Genome-wide burden-captured h2, computed directly from the full
        // accumulated burden score vs the full accumulated genetic score
        // (not the same as summing per-window h2_burden, for the same
        // reason TOTAL_SUM and TOTAL_ACTUAL can differ - see TOTAL_ACTUAL)
        double r2_total = compute_r2(burden_total, g_total);
        double h2_total_burden_actual = r2_total * h2_total_actual;

        hers_out << "TOTAL_SUM\tNA\tNA\tNA\tALL\tNA\t" << sum_window_h2
                  << "\tNA\t" << sum_window_h2_burden << "\n";
        hers_out << "TOTAL_ACTUAL\tNA\tNA\tNA\tALL\tNA\t" << h2_total_actual
                  << "\t" << r2_total << "\t" << h2_total_burden_actual << "\n";
        for (int c = 0; c < n_cat; ++c) {
            double mg   = g_total_c[c].mean();
            double h2a  = (g_total_c[c].array() - mg).square().sum() / (n_inds - 1);
            double r2a  = compute_r2(burden_total_c[c], g_total_c[c]);
            hers_out << "TOTAL_SUM\tNA\tNA\tNA\t" << cat_names[c] << "\tNA\t"
                      << sum_h2_c[c] << "\tNA\t" << sum_h2b_c[c] << "\n";
            hers_out << "TOTAL_ACTUAL\tNA\tNA\tNA\t" << cat_names[c] << "\tNA\t"
                      << h2a << "\t" << r2a << "\t" << r2a * h2a << "\n";
        }
        hers_out.close();

        Rcout << "Sum of per-window h2 (naive): " << sum_window_h2 << "\n";
        Rcout << "Actual h2 of combined genetic score: " << h2_total_actual << "\n";
        Rcout << "Sum of per-window burden-captured h2 (naive): " << sum_window_h2_burden << "\n";
        Rcout << "Actual burden-captured h2 (genome-wide r2 x h2): " << h2_total_burden_actual << "\n";
    }

    Rcout << "Done. " << n_written << " windows written";
    if (n_skipped > 0) Rcout << ", " << n_skipped << " skipped (empty)";
    if (write_burden_scores) {
        Rcout << "\nOutput: " << out_file << "\nWindow labels: " << kept_file;
    }
    if (has_effects) {
        Rcout << "\nHeritabilities: " << out_file << ".hers";
    }
    Rcout << "\n";
}

// Public R-callable function. Exactly one of kb_size, n_snps_per_window,
// or target_mac_per_ind should be positive; the others should be left at
// their default of -1 (meaning "not used"). effects_file is optional: if
// given, also writes <out_file>.hers with per-window and total heritability
// computed from the known true effect sizes (assumes var(phenotype) = 1).
// write_burden_scores = false skips writing the burden score matrix and
// .kept_windows file entirely, e.g. when only the .hers file is wanted.
// [[Rcpp::export]]
void compute_burden_windows(
    const std::string& bed_prefix,
    const std::string& out_file,
    double kb_size             = -1,   // strategy 1: window size in kb
    int    n_snps_per_window   = -1,   // strategy 2: fixed SNP count
    double target_mac_per_ind  = -1,   // strategy 3: target MAC/n_inds
    int    write_buffer_size   = 500,
    int    chunk_size          = 5000, // only used for strategy 3 MAC scan
    std::string effects_file   = "",   // optional: SNP, effect size
    bool   write_burden_scores = true,
    Rcpp::Nullable<Rcpp::IntegerMatrix>   annotation  = R_NilValue,  // n_snps x n_cat, 0/1
    Rcpp::Nullable<Rcpp::CharacterVector> annot_names = R_NilValue   // one name per column
) {
    // -- Read .bim --
    List bim_list = read_bim_file(bed_prefix);
    IntegerVector snp_chr = bim_list["chr"];
    NumericVector snp_pos = bim_list["pos"];
    int n_snps = snp_chr.size();
    BimData bim = read_bim_data(snp_chr, snp_pos);

    // -- Read .fam for n_inds --
    List fam = read_fam_file(bed_prefix);
    int n_inds = as<CharacterVector>(fam["iid"]).size();

    Rcout << "Individuals: " << n_inds << "\n";
    Rcout << "SNPs: "        << n_snps << "\n";

    // -- Optional effects file --
    std::vector<double> beta;   // stays empty if no effects_file given
    if (effects_file != "") {
        CharacterVector snp_id = bim_list["snp"];
        beta = read_effects_file(effects_file, snp_id);
    }

    // -- Optional annotation: every column is a category, as in he_reg_spa.
    //    Heritability is reported for the window as a whole ("ALL") and for
    //    each category separately. --
    std::vector<int> annot;
    std::vector<std::string> cat_names;
    if (annotation.isNotNull()) {
        Rcpp::IntegerMatrix A(annotation.get());
        if (A.nrow() != n_snps)
            Rcpp::stop("annotation must have one row per SNP in the .bim (" +
                       std::to_string(n_snps) + " rows), got " + std::to_string(A.nrow()));
        int n_cat = A.ncol();
        if (annot_names.isNotNull()) {
            Rcpp::CharacterVector nm(annot_names.get());
            if (nm.size() != n_cat)
                Rcpp::stop("annot_names length must equal the number of annotation columns");
            for (int c = 0; c < n_cat; ++c) cat_names.push_back(Rcpp::as<std::string>(nm[c]));
        } else {
            for (int c = 0; c < n_cat; ++c) cat_names.push_back("CAT_" + std::to_string(c + 1));
        }
        annot.resize((size_t) n_snps * n_cat);
        std::vector<long> per_cat(n_cat, 0);
        for (int c = 0; c < n_cat; ++c)
            for (int s = 0; s < n_snps; ++s) {
                int v = (A(s, c) != 0) ? 1 : 0;
                annot[(size_t) c * n_snps + s] = v;
                per_cat[c] += v;
            }
        Rcout << "Annotation: " << n_cat << " categories (";
        for (int c = 0; c < n_cat; ++c)
            Rcout << (c ? ", " : "") << cat_names[c] << ": " << per_cat[c];
        Rcout << " SNPs)\n";
    }

    // -- Check exactly one window strategy is specified --
    int n_strategies = (kb_size > 0) + (n_snps_per_window > 0) + (target_mac_per_ind > 0);
    if (n_strategies != 1)
        Rcpp::stop("Specify exactly one of: kb_size, n_snps_per_window, target_mac_per_ind");

    // -- Generate windows --
    std::vector<Window> windows;

    if (kb_size > 0) {
        Rcout << "Strategy: fixed window size (" << kb_size << " kb)\n";
        windows = windows_by_kb(bim, kb_size);
    } else if (n_snps_per_window > 0) {
        Rcout << "Strategy: fixed SNP count (" << n_snps_per_window << " SNPs/window)\n";
        windows = windows_by_nsnps(bim, n_snps_per_window);
    } else {
        Rcout << "Strategy: target variation (MAC/n_inds = " << target_mac_per_ind << ")\n";
        windows = windows_by_variation(bim, target_mac_per_ind, n_inds, bed_prefix, chunk_size);
    }

    Rcout << "Windows generated: " << windows.size() << "\n";

    // -- Compute and write burdens (and heritability, if effects given) --
    write_burdens(windows, bed_prefix, n_inds, n_snps, out_file, write_buffer_size, beta,
                  write_burden_scores, annot, cat_names);
}
