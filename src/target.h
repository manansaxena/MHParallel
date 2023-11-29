#include <boost/math/distributions/beta.hpp>
#include <limits>

double target(const Eigen::VectorXd& x){
    double alpha = 4; 
    double beta = 2;
    using boost::math::beta_distribution;

    double result = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double xi = x(i);
        if (xi < 0 || xi > 1) {
            return -std::numeric_limits<double>::infinity();  // Return negative infinity if out of bounds
        }
        beta_distribution<double> dist(alpha, beta);
        result += log(boost::math::pdf(dist, xi));  // Accumulate log probability
    }
    return result;
}




// double target(const Eigen::VectorXd& x, double alpha = 4, double beta = 2) {
//     using boost::math::beta_distribution;

//     double result = 0.0;
//     for (int i = 0; i < x.size(); ++i) {
//         double xi = x(i);
//         if (xi < 0 || xi > 1) {
//             return -std::numeric_limits<double>::infinity();  // Return negative infinity if out of bounds
//         }
//         beta_distribution<double> dist(alpha, beta);
//         result += log(boost::math::pdf(dist, xi));  // Accumulate log probability
//     }
//     return result;
// }

// #include <Eigen/Dense>
// #include <cmath>

// double target(const Eigen::VectorXd& x) {
//     using namespace Eigen;

//     Eigen::VectorXd mean(3);  // Mean vector of size 3
//     Eigen::MatrixXd covar(3, 3);  // Covariance matrix 3x3

//     // Initialize the mean vector
//     mean << 1.0, 2.0, 3.0; // Example values

//     // Initialize the covariance matrix
//     covar << 4.0, 2.0, 1.0,  // Row 1
//              2.0, 5.0, 3.0,  // Row 2
//              1.0, 3.0, 6.0;  // Row 3

//     const double log2pi = std::log(2 * M_PI);
//     double log_likelihood = 0.0;

//     LLT<MatrixXd> llt(covar);  // compute the Cholesky decomposition of covar
//     MatrixXd L = llt.matrixL();  // lower triangular matrix

//     VectorXd diff = x - mean;
//     VectorXd solve = L.triangularView<Lower>().solve(diff);
//     double quadform = solve.squaredNorm();

//     log_likelihood -= 0.5 * (covar.cols() * log2pi + 2 * L.diagonal().array().log().sum() + quadform);

//     return log_likelihood;
// }

