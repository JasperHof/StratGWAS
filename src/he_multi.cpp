// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>
#include <random>

using Eigen::Matrix2d;
using Eigen::Vector2d;

// [[Rcpp::export]]
Rcpp::NumericMatrix he_multi(const Eigen::MatrixXd& genotypes, const Rcpp::NumericMatrix& pheno) {
  
    int n_inds = genotypes.rows();
    int n_snps = genotypes.cols();
    int n_pheno = pheno.cols();

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

    // compute y'Ky and y'y - for ALL pairs of y!
    Rcpp::NumericMatrix gencors(n_pheno, n_pheno);

    for(int i = 0; i < n_pheno; ++i){
        for(int j = i; j < n_pheno; ++j){
            
            // map the i-th and j-th phenotype to Eigen vectors
            Eigen::VectorXd yi = Eigen::VectorXd::Map(pheno.column(i).begin(), pheno.nrow());
            Eigen::VectorXd yj = Eigen::VectorXd::Map(pheno.column(j).begin(), pheno.nrow());

             // compute Ky = X^T y_i
            Eigen::VectorXd Xy_i = genotypes.transpose() * yi;
            Eigen::VectorXd Xy_j = genotypes.transpose() * yj;
            double yKyt = Xy_i.dot(Xy_j)  / n_snps;
            double yyt = yi.dot(yj);

            // solve 2x2 HE equations for this pair
            Eigen::Matrix2d A;
            Eigen::Vector2d B, x;

            A << tr_KK, tr_K,
                 tr_K, n_inds;
       
            B << yKyt, yyt;

            x = A.colPivHouseholderQr().solve(B);
            gencors(i, j) = x[0];
            gencors(j, i) = x[0]; // symmetric
        }
    }

    return gencors;
}