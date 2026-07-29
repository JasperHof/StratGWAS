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
#include <cstdlib>
#include <map>
#include <unordered_map>

using namespace Rcpp;
using Eigen::VectorXd;

// ===========================================================================
// Enrichment-weighted burden association.
//
// Builds burden scores exactly like compute_burden_windows() in burden_size.cpp
// (same three window strategies: kb_size / n_snps_per_window / target_mac_per_ind),
// but with each SNP's allele count scaled by a per-category ENRICHMENT weight
// taken from the output of he_sliding_window_part() / reml_sliding_window_part()
// (he_sliding_window_enrich_alpha.cpp). Immediately after a window's (weighted)
// burden score is built, it is regressed against a user-supplied phenotype
// matrix (individual IDs as rownames, same convention as he_sliding_window_part()
// / linear_gwas_parallel()) -- one row of association output per window per
// phenotype column. Nothing is written to disk except the association result
// table: no burden-score matrix, no .kept_windows file (burdens live only in
// memory, per window).
//
// ENRICHMENT WEIGHTING (no LOCO, per current request): each category's
// genome-wide heritability enrichment is
//     enrichment_c = h2_c / (SNPs_c / SNPs_total)
// read directly off the TOTAL rows of the HE/REML output file (h2_c) and the
// summed per-window SNP counts (m_c) for that category. This is computed ONCE
// from the whole genome (not leave-one-chromosome-out) and applied directly to
// scale allele counts everywhere. A LOCO version would instead recompute
// enrichment_c once per chromosome, excluding that chromosome's windows from
// both the h2 and the expected-proportion sums, and use the matching map when
// building burdens for windows on that chromosome -- left as a future
// extension (see the note above compute_category_enrichment()).
//
// Everything below lives in an anonymous namespace (internal linkage) so its
// helper names (Window, BimData, windows_by_kb, ...) cannot collide with the
// identically-named, externally-linked symbols already defined in
// burden_size.cpp when both files are compiled together into the same package
// (same reasoning as documented at the top of he_sliding_window_enrich_alpha.cpp).
// ===========================================================================
namespace {

// ---------------------------------------------------------------------------
// Window + .bim plumbing (same conventions/logic as burden_size.cpp)
// ---------------------------------------------------------------------------
struct Window {
    int chr;
    int start_idx;    // 0-based, inclusive, index into .bim
    int end_idx;      // 0-based, exclusive (half-open)
    long bp_start;
    long bp_end;
};

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

        while (i < n && bim.chr[i] == this_chr && bim.pos[i] <= tile_end) ++i;

        w.end_idx = i;
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
std::vector<Window> windows_by_variation(
    const BimData& bim,
    double target_mac_per_ind,
    int n_inds,
    const std::string& bed_prefix,
    int chunk_size = 5000
) {
    std::vector<Window> windows;
    int n = bim.n_snps;

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
                mac[idx] = maf * nchrobs;
            }
            ++idx;
        }
        frq_in.close();
    } else {
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
                if (snp_mac > n_inds) snp_mac = 2.0 * n_inds - snp_mac;
                mac[snp_start + s] = snp_mac;
            }
            if ((c + 1) % 100 == 0)
                Rcout << "MAC scan: chunk " << c+1 << "/" << n_chunks << "\n";
        }
    }

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
// Parse a he_sliding_window_part() / reml_sliding_window_part() output file
// (long format: chr, start, end, n_left, n_right, n_common, phenotype,
// category, m_c, vg, se_vg, h2, se_h2, enrichment, fit, [converged, n_iter])
// and compute genome-wide per-category heritability enrichment for ONE
// phenotype, as the ratio of a category's SHARE OF HERITABILITY to its share
// of SNPs:
//     enrichment_c = [ h2_c / sum_c h2_c ] / [ sum_windows m_c / sum_all m_c ]
// Both numerator and denominator are proportions summing to 1, so an
// unenriched category gets enrichment 1.
//
// NOTE: h2_c must be normalized by the total h2 the model explains, NOT used
// raw. Using raw h2_c silently rescales every enrichment by 1/h2_total, so for
// a trait with h2_total ~ 0.5 all enrichments come out roughly halved; combined
// with the min_enrichment floor that collapses several categories onto exactly
// the floor value and leaves the rest in a narrow band, making the weighting
// nearly uniform (and the results nearly identical to an unweighted run).
//
// Falls back to 1.0 (no reweighting) for a category whose TOTAL h2 or expected
// proportion could not be determined. Negative / unstable enrichments (HE and
// REML h2 estimates can go negative) are clamped to `min_enrichment` (default
// 0, i.e. such a category never gets down-weighted below "contributes nothing").
//
// NOTE (LOCO): this reads the file exactly as produced by a single full-genome
// he_sliding_window_part()/reml_sliding_window_part() run, so a category's
// enrichment here is informed by windows on every chromosome, including
// whichever one is later tested for association below -- i.e. it is NOT
// leave-one-chromosome-out. A LOCO version would instead build one enrichment
// map per chromosome, each excluding that chromosome's own rows from both the
// h2 and expected-proportion sums, and select the matching map when weighting
// a window's SNPs. Skipped here per the current request (direct genome-wide
// weighting); this function returns a single genome-wide map.
// ---------------------------------------------------------------------------
std::unordered_map<std::string, double> compute_category_enrichment(
    const std::string& he_file,
    const std::string& trait_name,
    double min_enrichment,
    std::vector<std::string>& cat_order   // out: categories in first-seen order (for logging)
) {
    std::ifstream in(he_file);
    if (!in.is_open()) Rcpp::stop("Could not open HE/REML output file: " + he_file);

    std::string header_line;
    if (!std::getline(in, header_line)) Rcpp::stop("HE/REML output file is empty: " + he_file);

    std::vector<std::string> cols;
    {
        std::istringstream hs(header_line);
        std::string tok;
        while (std::getline(hs, tok, '\t')) cols.push_back(tok);
    }
    auto col_idx = [&](const std::string& name) -> int {
        for (size_t k = 0; k < cols.size(); ++k) if (cols[k] == name) return (int)k;
        Rcpp::stop("HE/REML output file is missing expected column: " + name);
        return -1;
    };
    int i_chr = col_idx("chr"), i_ph = col_idx("phenotype"), i_cat = col_idx("category");
    int i_mc  = col_idx("m_c"), i_h2 = col_idx("h2");
    int max_needed = i_chr;
    max_needed = std::max(max_needed, i_ph);
    max_needed = std::max(max_needed, i_cat);
    max_needed = std::max(max_needed, i_mc);
    max_needed = std::max(max_needed, i_h2);

    std::unordered_map<std::string, double> total_h2;
    std::unordered_map<std::string, double> sum_mc;
    double sum_mc_all = 0.0;
    std::unordered_map<std::string, bool> seen_cat;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        {
            std::istringstream ls(line);
            std::string tok;
            while (std::getline(ls, tok, '\t')) f.push_back(tok);
        }
        if ((int)f.size() <= max_needed) continue;
        if (f[i_ph] != trait_name) continue;

        const std::string& cat = f[i_cat];
        if (!seen_cat.count(cat)) { seen_cat[cat] = true; cat_order.push_back(cat); }

        if (f[i_chr] == "TOTAL") {
            if (f[i_h2] != "NA") {
                try { total_h2[cat] = std::stod(f[i_h2]); } catch (...) {}
            }
        } else {
            if (f[i_mc] != "NA") {
                try {
                    double mc = std::stod(f[i_mc]);
                    sum_mc[cat] += mc;
                    sum_mc_all  += mc;
                } catch (...) {}
            }
        }
    }
    in.close();

    if (sum_mc_all <= 0.0)
        Rcpp::stop("Could not find any per-window rows for phenotype '" + trait_name + "' in " + he_file);

    auto fmt = [](double v) { return ISNAN(v) ? std::string("NA") : std::to_string(v); };

    // Total h2 the model explains, summed over categories. Enrichment compares
    // SHARE of h2 against SHARE of SNPs, so this is the correct denominator.
    double h2_total = 0.0;
    int n_h2 = 0;
    for (const auto& c : cat_order) {
        auto it = total_h2.find(c);
        if (it != total_h2.end() && !ISNAN(it->second)) { h2_total += it->second; ++n_h2; }
    }
    if (n_h2 == 0)
        Rcpp::stop("No TOTAL rows with a usable h2 found for phenotype '" + trait_name +
                   "' in " + he_file);
    if (h2_total <= 0.0)
        Rcpp::stop("Total h2 across categories is not positive (" + std::to_string(h2_total) +
                   ") for phenotype '" + trait_name + "' -- cannot form enrichment ratios");

    std::unordered_map<std::string, double> enrichment;
    Rcout << "Category enrichment (phenotype = " << trait_name
          << ", total h2 = " << h2_total << "):\n";
    for (const auto& c : cat_order) {
        double h2c = total_h2.count(c) ? total_h2[c] : NA_REAL;
        double mc  = sum_mc.count(c) ? sum_mc[c] : 0.0;
        double expected = (mc > 0) ? mc / sum_mc_all : NA_REAL;
        double share    = ISNAN(h2c) ? NA_REAL : h2c / h2_total;   // share of h2
        double enr = 1.0;
        bool clamped = false;
        if (!ISNAN(share) && !ISNAN(expected) && expected > 0) {
            enr = share / expected;
            if (enr < min_enrichment) { enr = min_enrichment; clamped = true; }
        }
        enrichment[c] = enr;
        Rcout << "  " << c << ": h2=" << fmt(h2c)
              << " h2_share=" << fmt(share)
              << " expected_prop=" << fmt(expected)
              << " enrichment=" << enr;
        if (clamped) Rcout << "  (clamped at min_enrichment)";
        Rcout << "\n";
    }
    return enrichment;
}

// ---------------------------------------------------------------------------
// Read a two-column SNP -> category-name file (whitespace separated, no
// header) and translate it, via the enrichment map, into a per-SNP weight
// aligned to .bim order. SNPs absent from the file (or whose category has no
// enrichment entry) fall back to the "uncategorized" category's enrichment if
// available, else weight 1.0 (no reweighting).
// ---------------------------------------------------------------------------
std::vector<double> read_snp_weights(
    const std::string& category_file,
    const CharacterVector& snp_ids,
    const std::unordered_map<std::string, double>& enrichment
) {
    int n = snp_ids.size();
    std::unordered_map<std::string, std::string> snp_cat;
    std::ifstream in(category_file);
    if (!in.is_open()) Rcpp::stop("Could not open category file: " + category_file);
    std::string snp, cat;
    long n_lines = 0;
    std::map<std::string, long> cnt_file;      // labels present in the category file
    while (in >> snp >> cat) { snp_cat[snp] = cat; ++n_lines; cnt_file[cat]++; }
    in.close();

    double default_w = enrichment.count("uncategorized") ? enrichment.at("uncategorized") : 1.0;

    std::vector<double> w(n, default_w);
    long n_listed = 0, n_unlisted = 0;
    // Per-label tallies of the SNPs actually used, so that a silent label
    // mismatch between the category file and the HE/REML output shows up in the
    // log rather than quietly defaulting every SNP to the fallback weight.
    std::map<std::string, long> cnt_known, cnt_unknown;
    for (int i = 0; i < n; ++i) {
        std::string s = as<std::string>(snp_ids[i]);
        auto it = snp_cat.find(s);
        if (it == snp_cat.end()) { ++n_unlisted; continue; }   // keeps default_w
        ++n_listed;
        auto eit = enrichment.find(it->second);
        if (eit != enrichment.end()) { w[i] = eit->second; cnt_known[it->second]++; }
        else                         { w[i] = default_w;   cnt_unknown[it->second]++; }
    }

    Rcout << "Category file: " << category_file << "\n";
    Rcout << "  " << n_lines << " lines read, " << cnt_file.size() << " distinct labels: ";
    { bool first = true;
      for (const auto& kv : cnt_file) { if (!first) Rcout << ", "; Rcout << kv.first; first = false; } }
    Rcout << "\n";
    Rcout << "  " << n_listed << "/" << n << " .bim SNPs found in the category file";
    if (n_unlisted > 0)
        Rcout << ", " << n_unlisted << " not listed -> fallback weight " << default_w;
    Rcout << "\n";

    Rcout << "SNPs per category (weight applied to allele counts):\n";
    for (const auto& kv : cnt_known)
        Rcout << "  " << kv.first << ": " << kv.second << " SNPs, weight = "
              << enrichment.at(kv.first) << "\n";
    if (n_unlisted > 0)
        Rcout << "  <not in category file>: " << n_unlisted << " SNPs, weight = " << default_w << "\n";

    if (!cnt_unknown.empty()) {
        Rcout << "WARNING: these labels appear in the category file but have NO matching\n"
              << "         category in the HE/REML output -- their SNPs silently fell back\n"
              << "         to weight " << default_w << ". Check that the category names match exactly\n"
              << "         (whitespace, capitalisation) the cat_names used in the HE/REML run:\n";
        for (const auto& kv : cnt_unknown)
            Rcout << "  '" << kv.first << "': " << kv.second << " SNPs\n";
    }

    // Flag the degenerate case where every SNP ends up with the same weight --
    // the burden is then just an unweighted allele count and the run carries no
    // enrichment information at all.
    double wmin = w.empty() ? 0.0 : w[0], wmax = wmin;
    for (int i = 1; i < n; ++i) { if (w[i] < wmin) wmin = w[i]; if (w[i] > wmax) wmax = w[i]; }
    Rcout << "Applied SNP weights: min = " << wmin << ", max = " << wmax << "\n";
    if (wmax - wmin < 1e-12)
        Rcout << "WARNING: every SNP received the SAME weight (" << wmin << ") -- the weighted\n"
              << "         burden is proportional to an unweighted allele count, so results will\n"
              << "         be identical to an unweighted/uncategorized run.\n";
    return w;
}

// ---------------------------------------------------------------------------
// Simple univariate linear regression: y ~ burden. Same statistics as the SNP
// regression in linear_parallel.cpp (burden score standardized to unit
// variance over the NON-MISSING individuals for this phenotype column, Wald
// chisq, two-sided p-value). `has_y` masks out individuals missing THIS
// phenotype column (columns of pheno_mat can have different missingness).
// ---------------------------------------------------------------------------
struct RegResult { double beta, se, chisq, pval; int n; bool ok; };

RegResult regress_burden(const VectorXd& burden, const std::vector<double>& y,
                         const std::vector<char>& has_y) {
    RegResult r; r.ok = false; r.beta = r.se = r.chisq = r.pval = NA_REAL; r.n = 0;
    int n_inds = (int) burden.size();

    double y_sum = 0.0; int n_y = 0;
    for (int i = 0; i < n_inds; ++i) if (has_y[i]) { y_sum += y[i]; ++n_y; }
    if (n_y < 4) return r;
    double y_mean = y_sum / n_y;

    double b_sum = 0.0;
    for (int i = 0; i < n_inds; ++i) if (has_y[i]) b_sum += burden(i);
    double b_mean = b_sum / n_y;

    double b_var = 0.0;
    for (int i = 0; i < n_inds; ++i) if (has_y[i]) { double d = burden(i) - b_mean; b_var += d * d; }
    b_var /= (n_y - 1);
    if (b_var <= 1e-12) return r;   // no variance in burden (e.g. all-zero window)
    double b_sd = std::sqrt(b_var);

    double XtY = 0.0, XtX = 0.0, yTy = 0.0;
    for (int i = 0; i < n_inds; ++i) {
        if (!has_y[i]) continue;
        double xs = (burden(i) - b_mean) / b_sd;
        double yc = y[i] - y_mean;
        XtY += xs * yc;
        XtX += xs * xs;
        yTy += yc * yc;
    }
    if (XtX <= 0) return r;

    double beta = XtY / XtX;
    double rss  = yTy - beta * XtY;
    if (rss <= 0) return r;

    double se = std::sqrt((rss / (n_y - 2)) / XtX);
    double chisq = (beta / se) * (beta / se);
    double pval = std::erfc(std::sqrt(chisq / 2.0));

    r.beta = beta; r.se = se; r.chisq = chisq; r.pval = pval; r.n = n_y; r.ok = true;
    return r;
}

}  // end anonymous namespace (internal linkage)

// ===========================================================================
// Public R-callable function.
//
// Builds enrichment-weighted burden scores window by window (same three
// strategies as compute_burden_windows() in burden_size.cpp -- exactly one of
// kb_size / n_snps_per_window / target_mac_per_ind should be positive) and,
// for each window, immediately regresses the burden score on `pheno_mat`.
// Only the association result table is written to disk; no burden-score
// matrix or window-membership file is produced. Since the (enrichment-
// weighted) burden itself does not depend on phenotype, it is built ONCE per
// window and then regressed against every column of pheno_mat in turn.
//
//   bed_prefix     - genotype fileset the burdens are built from
//   he_file        - out_file produced by he_sliding_window_part() /
//                     reml_sliding_window_part(), used to derive per-category
//                     genome-wide heritability enrichment
//   category_file  - two whitespace-separated columns, SNP then category
//                     name, aligned to bed_prefix's .bim (using the same
//                     category names as cat_names passed to the HE/REML run);
//                     SNPs not listed are treated as "uncategorized"
//   trait_name     - which phenotype column of he_file to use for enrichment
//                     (independent of which pheno_mat column(s) are tested)
//   pheno_mat      - numeric matrix of phenotype(s) with individual IDs as
//                     rownames (same convention as he_sliding_window_part() /
//                     linear_gwas_parallel()); one or more columns, matched to
//                     the .fam by IID; NA entries are excluded per column
//   out_file       - association output, long format: one row per window per
//                     phenotype column
//   min_enrichment - floor applied to (possibly negative/unstable) estimated
//                     enrichment weights; default 0 (never down-weight below
//                     "contributes nothing" to the burden)
// [[Rcpp::export]]
void burden_enrich_association(
    const std::string& bed_prefix,
    const std::string& he_file,
    const std::string& category_file,
    const std::string& trait_name,
    const SEXP pheno_mat,
    const std::string& out_file,
    double kb_size             = -1,   // strategy 1: window size in kb
    int    n_snps_per_window   = -1,   // strategy 2: fixed SNP count
    double target_mac_per_ind  = -1,   // strategy 3: target MAC/n_inds
    double min_enrichment      = 0.0,
    int    write_buffer_size   = 500,
    int    chunk_size          = 5000  // only used for strategy 3 MAC scan
) {
    // -- .bim / .fam for the burden fileset --
    List bim_list = read_bim_file(bed_prefix);
    IntegerVector snp_chr = bim_list["chr"];
    NumericVector snp_pos = bim_list["pos"];
    CharacterVector snp_id = bim_list["snp"];
    int n_snps = snp_chr.size();
    BimData bim = read_bim_data(snp_chr, snp_pos);

    List fam = read_fam_file(bed_prefix);
    CharacterVector geno_iid = fam["iid"];
    int n_total_inds = geno_iid.size();

    Rcout << "Individuals (.fam): " << n_total_inds << "\n";
    Rcout << "SNPs (.bim): "        << n_snps << "\n";

    // -- Per-category enrichment (genome-wide, non-LOCO) and per-SNP weight --
    std::vector<std::string> cat_order;
    std::unordered_map<std::string, double> enrichment =
        compute_category_enrichment(he_file, trait_name, min_enrichment, cat_order);
    std::vector<double> weight = read_snp_weights(category_file, snp_id, enrichment);

    // -- Phenotype matrix, matched to .fam by IID (rownames of pheno_mat) --
    if (!Rf_isMatrix(pheno_mat) || Rf_isNull(rownames(pheno_mat)))
        Rcpp::stop("Phenotype must be a numeric matrix with individual IDs as rownames");
    NumericMatrix pheno = as<NumericMatrix>(pheno_mat);
    CharacterVector pheno_ids = rownames(pheno_mat);
    int n_pheno = pheno.ncol();

    CharacterVector pheno_names;
    SEXP cn = colnames(pheno_mat);
    if (!Rf_isNull(cn)) pheno_names = cn;
    if (pheno_names.size() != n_pheno) {
        pheno_names = CharacterVector(n_pheno);
        for (int j = 0; j < n_pheno; ++j) pheno_names[j] = "trait" + std::to_string(j + 1);
    }

    IntegerVector match_idx = match(geno_iid, pheno_ids);
    std::vector<int> geno_keep;
    std::vector<int> pheno_row;
    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            geno_keep.push_back(i);
            pheno_row.push_back(match_idx[i] - 1);
        }
    }
    if (geno_keep.empty()) Rcpp::stop("No overlapping individuals between genotype and phenotype");
    int n_inds = (int) geno_keep.size();
    Rcout << "Individuals with phenotype data: " << n_inds << "\n";
    Rcout << "Phenotype columns: " << n_pheno << "\n";

    // Per-phenotype value + missingness mask, reordered to match geno_keep
    std::vector< std::vector<double> > Y(n_pheno, std::vector<double>(n_inds));
    std::vector< std::vector<char> >   Ymask(n_pheno, std::vector<char>(n_inds));
    for (int t = 0; t < n_pheno; ++t) {
        for (int i = 0; i < n_inds; ++i) {
            double v = pheno(pheno_row[i], t);
            Y[t][i] = ISNAN(v) ? 0.0 : v;
            Ymask[t][i] = ISNAN(v) ? 0 : 1;
        }
    }

    // -- Check exactly one window strategy is specified --
    int n_strategies = (kb_size > 0) + (n_snps_per_window > 0) + (target_mac_per_ind > 0);
    if (n_strategies != 1)
        Rcpp::stop("Specify exactly one of: kb_size, n_snps_per_window, target_mac_per_ind");

    std::vector<Window> windows;
    if (kb_size > 0) {
        Rcout << "Strategy: fixed window size (" << kb_size << " kb)\n";
        windows = windows_by_kb(bim, kb_size);
    } else if (n_snps_per_window > 0) {
        Rcout << "Strategy: fixed SNP count (" << n_snps_per_window << " SNPs/window)\n";
        windows = windows_by_nsnps(bim, n_snps_per_window);
    } else {
        Rcout << "Strategy: target variation (MAC/n_inds = " << target_mac_per_ind << ")\n";
        windows = windows_by_variation(bim, target_mac_per_ind, n_total_inds, bed_prefix, chunk_size);
    }
    Rcout << "Windows generated: " << windows.size() << "\n";

    // -- Association output (burdens themselves are never written to disk) --
    std::ofstream out(out_file);
    if (!out.is_open()) Rcpp::stop("Could not open output file: " + out_file);
    out << "Window\tChr\tBP_start\tBP_end\tN_SNPs\tSum_Weight\tPhenotype\tBETA\tSE\tChisq\tP\tN\n";

    std::vector<std::string> row_buf;
    row_buf.reserve(write_buffer_size * n_pheno);

    int n_windows = (int) windows.size();
    int n_written = 0, n_skipped = 0;

    for (int wi = 0; wi < n_windows; ++wi) {
        const Window& w = windows[wi];
        int n_snps_win = w.end_idx - w.start_idx;
        if (n_snps_win == 0) { n_skipped++; continue; }

        IntegerMatrix geno_block = readBedBlock(
            bed_prefix + ".bed", n_total_inds, n_snps,
            0, n_total_inds - 1,
            w.start_idx, w.end_idx - 1   // inclusive end for readBedBlock
        );

        // Weighted burden score for THIS window, held only in memory. Built
        // ONCE (it does not depend on phenotype) and reused for every trait.
        VectorXd burden_vec = VectorXd::Zero(n_inds);
        double sum_weight = 0.0;
        for (int s = 0; s < n_snps_win; ++s) {
            double wt = weight[w.start_idx + s];
            sum_weight += wt;
            if (wt == 0.0) continue;
            for (int i = 0; i < n_inds; ++i) {
                int g = geno_block(geno_keep[i], s);
                if (g != -1) burden_vec(i) += wt * g;
            }
        }

        std::ostringstream wlabel;
        wlabel << w.chr << "_" << w.bp_start << "_" << w.bp_end << "_" << n_snps_win;

        for (int t = 0; t < n_pheno; ++t) {
            RegResult reg = regress_burden(burden_vec, Y[t], Ymask[t]);

            std::ostringstream row;
            row << wlabel.str() << "\t" << w.chr << "\t" << w.bp_start << "\t" << w.bp_end
                << "\t" << n_snps_win << "\t" << sum_weight << "\t" << as<std::string>(pheno_names[t]) << "\t";
            if (reg.ok) {
                row << reg.beta << "\t" << reg.se << "\t" << reg.chisq << "\t" << reg.pval << "\t" << reg.n;
            } else {
                row << "NA\tNA\tNA\tNA\t" << reg.n;
            }
            row_buf.push_back(row.str());
        }
        n_written++;

        if ((int) row_buf.size() >= write_buffer_size * n_pheno) {
            for (const auto& r : row_buf) out << r << "\n";
            row_buf.clear();
        }

        if ((wi + 1) % 1000 == 0 || wi + 1 == n_windows)
            Rcout << "Processed " << n_written << "/" << n_windows << " windows\n";
    }

    for (const auto& r : row_buf) out << r << "\n";
    out.close();

    Rcout << "Done. " << n_written << " windows tested";
    if (n_skipped > 0) Rcout << ", " << n_skipped << " skipped (empty)";
    Rcout << "\nAssociation output: " << out_file << "\n";
}
