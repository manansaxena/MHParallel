#include "RcppEigen.h"
#include "Rcpp.h"
#include "mh_general.h"
#include <string>
#include <vector>

using namespace Rcpp;
 
//[[Rcpp::export]]
List metro_interface(int num_chains,
                     std::vector<Eigen::VectorXd> initial_states, 
                     int num_steps, 
                     int seed, 
                     int n_cores,
                     double scale_factor_of_proposal = 1.0)
{   
    return metropolis_hastings(num_chains,initial_states, num_steps,seed,n_cores, scale_factor_of_proposal);
}
