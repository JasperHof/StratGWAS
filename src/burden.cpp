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

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────
// Gene region as read from the position file
// ─────────────────────────────────────────────────────────────────────────
struct GeneRegion {
    std::string name;
    int chr;
    long start;
    long end;
    char strand;
};

// ─────────────────────────────────────────────────────────────────────────
// Read the gene position file: <name> <chr> <start> <end> <strand>
// ─────────────────────────────────────────────────────────────────────────
std::vector<GeneRegion> read_gene_file(const std::string& path) {
    std::vector<GeneRegion> genes;
    std::ifstream in(path);
    if (!in.is_open()) Rcpp::stop("Could not open gene position file: " + path);

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        GeneRegion g;
        std::string chr_str, strand_str;
        if (!(iss >> g.name >> chr_str >> g.start >> g.end >> strand_str))
            continue;

        bool is_numeric = !chr_str.empty() &&
            std::all_of(chr_str.begin(), chr_str.end(), ::isdigit);
        g.chr = is_numeric ? std::stoi(chr_str) : -1;
        g.strand = strand_str.empty() ? '+' : strand_str[0];

        genes.push_back(g);
    }
    in.close();
    return genes;  // no sorting needed now — we look up indices per gene directly
}

// ─────────────────────────────────────────────────────────────────────────
// SNP position lookup table built once from the .bim file. Assumes the
// .bim file is sorted by (chr, pos) within each chromosome, as is standard
// for PLINK files — this is what makes binary search valid.
// ─────────────────────────────────────────────────────────────────────────
struct BimIndex {
    std::vector<int>  chr;
    std::vector<long> pos;
    int n_snps;

    BimIndex(const IntegerVector& chr_in, const NumericVector& pos_in) {
        n_snps = chr_in.size();
        chr.resize(n_snps);
        pos.resize(n_snps);
        for (int i = 0; i < n_snps; ++i) {
            chr[i] = chr_in[i];
            pos[i] = (long)pos_in[i];
        }
    }

    // Finds the first SNP index whose (chr, pos) is >= (target_chr, target_pos).
    // Standard lower_bound over a comparator on (chr, pos) pairs.
    int lower_bound_idx(int target_chr, long target_pos) const {
        int lo = 0, hi = n_snps;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            bool mid_before_target =
                (chr[mid] < target_chr) ||
                (chr[mid] == target_chr && pos[mid] < target_pos);
            if (mid_before_target) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    // Finds the first SNP index whose (chr, pos) is > (target_chr, target_pos).
    int upper_bound_idx(int target_chr, long target_pos) const {
        int lo = 0, hi = n_snps;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            bool mid_after_or_eq_target =
                (chr[mid] > target_chr) ||
                (chr[mid] == target_chr && pos[mid] > target_pos);
            if (mid_after_or_eq_target) hi = mid;
            else lo = mid + 1;
        }
        return lo;
    }

    // Returns [start_idx, end_idx) — the half-open range of SNP indices
    // whose position falls within [gene.start, gene.end] on gene.chr.
    std::pair<int,int> snp_range_for_gene(const GeneRegion& gene) const {
        int start_idx = lower_bound_idx(gene.chr, gene.start);
        int end_idx   = upper_bound_idx(gene.chr, gene.end);
        return {start_idx, end_idx};
    }
};

//
// Computes per-gene allele-count burden scores using direct index lookup:
// for each gene, binary-search the (sorted) .bim positions to find the
// contiguous block of SNP indices falling inside the gene's coordinates,
// then read only that block from the .bed file and sum allele counts
// across individuals. Genes with zero overlapping SNPs are written as a
// row of zeros (not skipped), so the output always has one row per gene.
//
// Output: tab-separated text file, one row per gene (GENE name first
// column), one column per individual (IID), written incrementally so the
// full n_genes x n_inds matrix is never held in memory at once.
// [[Rcpp::export]]
void compute_gene_burden(
    const std::string& bed_prefix,
    const std::string& gene_file,
    const std::string& out_file,
    int write_buffer_genes = 500   // flush to disk every N gene rows
) {
    // ── Read .fam ─────────────────────────────────────────────────────────
    List fam = read_fam_file(bed_prefix);
    CharacterVector geno_iid = fam["iid"];
    int n_inds = geno_iid.size();

    // ── Read .bim ─────────────────────────────────────────────────────────
    List bim = read_bim_file(bed_prefix);
    IntegerVector snp_chr = bim["chr"];
    NumericVector snp_pos = bim["pos"];   // bp position; adjust field name if different
    int n_snps = snp_chr.size();

    Rcout << "Individuals: " << n_inds << "\n";
    Rcout << "SNPs: "        << n_snps << "\n";

    BimIndex bim_index(snp_chr, snp_pos);

    // ── Read gene file ───────────────────────────────────────────────────
    std::vector<GeneRegion> genes = read_gene_file(gene_file);
    int n_genes = (int)genes.size();
    Rcout << "Genes: " << n_genes << "\n";
    if (n_genes == 0) Rcpp::stop("No genes read from gene position file");

    // ── Open output, write header ────────────────────────────────────────
    std::ofstream out(out_file);
    if (!out.is_open()) Rcpp::stop("Could not open output file for writing: " + out_file);

    out << "GENE";
    for (int i = 0; i < n_inds; ++i) out << "\t" << std::string(geno_iid[i]);
    out << "\n";

    // ── Buffer for batching disk writes (avoids one fstream flush per gene) ─
    std::vector<std::string> write_buffer;
    write_buffer.reserve(write_buffer_genes);

    int n_zero_snp_genes = 0;

    for (int gi = 0; gi < n_genes; ++gi) {

        const GeneRegion& gene = genes[gi];

        // Skip genes with unrecognized (non-numeric) chromosome codes —
        // report once at the end rather than erroring out mid-run
        if (gene.chr < 0) {
            std::ostringstream row;
            row << gene.name;
            for (int i = 0; i < n_inds; ++i) row << "\t0";
            write_buffer.push_back(row.str());
            n_zero_snp_genes++;
            continue;
        }

        auto range = bim_index.snp_range_for_gene(gene);
        int start_idx = range.first;
        int end_idx   = range.second;     // half-open: [start_idx, end_idx)
        int n_snps_gene = end_idx - start_idx;

        VectorXd burden_vec = VectorXd::Zero(n_inds);

        if (n_snps_gene > 0) {
            // Direct read of just this gene's SNP block — no scanning of
            // unrelated chunks, no full-genome pass per gene.
            IntegerMatrix geno_block = readBedBlock(
                bed_prefix + ".bed", n_inds, n_snps,
                0, n_inds - 1,
                start_idx, end_idx - 1   // inclusive end index for readBedBlock
            );

            for (int s = 0; s < n_snps_gene; ++s) {
                for (int i = 0; i < n_inds; ++i) {
                    int allele_count = geno_block(i, s);
                    if (allele_count != -1) {     // skip missing, no imputation
                        burden_vec(i) += allele_count;
                    }
                }
            }
        } else {
            n_zero_snp_genes++;
        }

        // ── Build this gene's output row ─────────────────────────────────
        std::ostringstream row;
        row << gene.name;
        for (int i = 0; i < n_inds; ++i) row << "\t" << burden_vec(i);
        write_buffer.push_back(row.str());

        // ── Flush buffer periodically ────────────────────────────────────
        if ((int)write_buffer.size() >= write_buffer_genes) {
            for (const auto& r : write_buffer) out << r << "\n";
            write_buffer.clear();
        }

        if ((gi + 1) % 1000 == 0 || gi + 1 == n_genes) {
            Rcout << "Processed " << gi + 1 << "/" << n_genes << " genes\n";
        }
    }

    // Flush any remaining buffered rows
    for (const auto& r : write_buffer) out << r << "\n";
    write_buffer.clear();

    out.close();

    Rcout << "Done. " << n_zero_snp_genes << " genes had zero overlapping SNPs.\n";
    Rcout << "Output written to " << out_file << "\n";
}