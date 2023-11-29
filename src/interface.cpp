#include "RcppEigen.h"
#include "Rcpp.h"
#include "mh_general.h"
#include "target.h"

//[[Rcpp::export]]
Rcpp::List metro_interface(int num_chains,
                           Eigen::MatrixXd initial_states, 
                           int num_steps, 
                           int seed, 
                           int n_cores,
                           double scale_factor_of_proposal,
                           Rcpp::Nullable<Eigen::MatrixXd> covar_object = R_NilValue) {

    std::tuple<std::vector<Eigen::MatrixXd>,Eigen::VectorXd,int> results;
    if(covar_object.isNotNull()) {
        if(Rcpp::as<Eigen::MatrixXd>(covar_object).rows() != 1 || Rcpp::as<Eigen::MatrixXd>(covar_object).cols() != 1){
            if(Rcpp::as<Eigen::MatrixXd>(covar_object).rows() != initial_states.rows() || Rcpp::as<Eigen::MatrixXd>(covar_object).cols() != initial_states.rows()){
                    Rcpp::stop("The covariance matrix should be a square matrix with size equal to the dimension of each state");
            }
        }
        results = metropolis_hastings(num_chains, initial_states, num_steps, seed, n_cores, target, scale_factor_of_proposal, Rcpp::as<Eigen::MatrixXd>(covar_object));
    }
    else {
        results = metropolis_hastings(num_chains, initial_states, num_steps, seed, n_cores, target, scale_factor_of_proposal);
    }
    
    
    return Rcpp::List::create(
        Rcpp::Named("particles") = std::get<0>(results),
        Rcpp::Named("accepted_counts") = std::get<1>(results),
        Rcpp::Named("timer") = std::get<2>(results)
    );
}