// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>
#include <random>

using Eigen::Matrix2d;
using Eigen::Vector2d;

// [[Rcpp::export]]
Rcpp::NumericVector he(const Eigen::MatrixXd& genotypes, const Rcpp::NumericVector& pheno) {
  
    int n_inds = genotypes.rows();
    int n_snps = genotypes.cols();

    double tr_K = n_inds; // because of standardization
    double tr_KK = 0.0;

    std::mt19937 rng(12345);
    std::normal_distribution<double> norm(0,1);
    Eigen::VectorXd z(n_inds);                      // dynamic length double vectors
    Eigen::VectorXd Kz(n_inds);
    int nmcmc = 20;

    for (int r = 0; r < nmcmc; ++r) {
        // 1. generate random vector z ~ N(0,1)
        for (int i = 0; i < n_inds; ++i) z[i] = norm(rng);

        // 2. compute K z = X (X^T z)
        Eigen::VectorXd t = genotypes.transpose() * z;
        Kz = genotypes * t;

        // 3. accumulate z^T K^2 z = ||K z||^2
        tr_KK += Kz.squaredNorm() / (nmcmc * n_snps * n_snps);
    }

    Rcpp::Rcout << "solution for tr(KK): " << tr_KK << "\n";

    // compute y'Ky and y'y
    Eigen::VectorXd Ky = genotypes.transpose() * Eigen::VectorXd::Map(pheno.begin(), pheno.size());
    double y_K_yt = Ky.squaredNorm() / n_snps;
    double y_yt = n_inds; // assuming standardized

    Rcpp::Rcout << "solution for y_K_yt: " << y_K_yt << "\n";

    // solve the HE equations
    Matrix2d A;
    Vector2d B, x;

    A << tr_KK, tr_K,
         tr_K, n_inds;
    B << y_K_yt, y_yt;

    x = A.colPivHouseholderQr().solve(B);

    Rcpp::Rcout << "Solution (h2, e2): " << x.transpose() << "\n";

    return Rcpp::NumericVector::create(x[0], x[1]);
}