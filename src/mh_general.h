#include <cmath>
#include <iostream>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <vector>
#include <omp.h>
#include "RcppEigen.h"
#include "Rcpp.h"
#include "Rcpp/Benchmark/Timer.h"

#include "target.h"

/*
TODO    
let the user input a covariance matrix and a scale factor. If they dont input a covariacne matrix then take the identity as default
dont use std::vector isntead use a matrix
so remove any rcpp dependencies from the core c++ code. the use should only have to adda wrapper in interface.h
which would take your function as a parameter and return everything the way its supposed to be
in the end the user should only need to change target.cpp and interface.cpp.


*/
// [[Rcpp::plugins(openmp)]]


/*
Input : 
    1) current_state: vector of current state of the sample in a metropolis chain
    2) random_vector: random vector of eq. distributed standard normal
Ouput :
    vector of new proposals
*/
Eigen::VectorXd generate_proposal(Eigen::VectorXd current_state, Eigen::VectorXd random_vector){
    return current_state + random_vector;
}


template<class RNG>
Eigen::VectorXd generate_random_vector(RNG& rng, int dim, Eigen::MatrixXd covar){
    Eigen::LLT<Eigen::MatrixXd> llt(covar);
    Eigen::MatrixXd L = llt.matrixL();
    assert(llt.info() != Eigen::NumericalIssue);
                
    boost::normal_distribution<> norm(0, 1);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > randn(rng, norm);
                
    Eigen::VectorXd z(dim);
    for (int i = 0; i < dim; ++i) {
        z(i) = randn();
    }
    return L * z;
}


template<class RNG>
Eigen::VectorXd mh_step(RNG& rng, Eigen::VectorXd current_state, double* proposal_accepted_cnt, double scale_factor_of_proposal){
    Eigen::MatrixXd covar = scale_factor_of_proposal * Eigen::MatrixXd::Identity(current_state.size(),current_state.size());
    Eigen::VectorXd random_vector = generate_random_vector(rng, current_state.size(), covar);
    Eigen::VectorXd proposed_state = generate_proposal(current_state, random_vector);

    double acceptance_ratio = std::min(double(1.0), std::exp(target(proposed_state) - target(current_state)));
    boost::random::uniform_real_distribution<double> uni_gen(0.0, 1.0);

    if(uni_gen(rng)<acceptance_ratio){
        *proposal_accepted_cnt += 1;
        return proposed_state;
    }
    else{
        return current_state; 
    }
}


/*
Input : 
    1) num_chains : number of chains to be run in the metropolis algorithm
    2) initial_states: a vector of vectors of intial states. Shape: num_chains X size of each initial state
    3) num_steps : number of steps each chain should take during the metropolis algorithm
    4) seed: seed used to initialize
    5) n_cores : number of cores
    6) target_ll : a Rcpp functional type to calculate the log likelihood of our target distribution
*/
Rcpp::List metropolis_hastings (int num_chains,
                                std::vector<Eigen::VectorXd> initial_states,
                                int num_steps,
                                int seed,
                                int n_cores,
                                double scale_factor_of_proposal
                                )

{
    Rcpp::Timer timer;
    timer.step("Overall Start");

    Rcpp::Rcout << "scale_factor_of_proposal: " << scale_factor_of_proposal << std::endl;
    Rcpp::List out(3);
    out.names() = Rcpp::CharacterVector::create("particles","accepted_counts","timer");

    if (n_cores > 0) {
      omp_set_num_threads(n_cores);
    } else {
      omp_set_num_threads(omp_get_max_threads());
    }

    // num_chains * (length of each state vector * num_steps) 
    std::vector<Eigen::MatrixXd> final_states(num_chains);
    std::vector<double> accepted_count(num_chains,0.0);
    
    #pragma omp parallel for 
    for(int chain_id = 0; chain_id < num_chains; chain_id++){
        boost::random::mt19937 rng(seed);
        rng.discard(omp_get_thread_num()*num_steps);

        Eigen::MatrixXd interim_states(initial_states[chain_id].size(),num_steps+1);
        interim_states.col(0) = initial_states[chain_id];
        double proposal_accepted_cnt = 0.0;
        for(int i=0;i<num_steps;i++){
            interim_states.col(i+1) = mh_step(rng, interim_states.col(i),&proposal_accepted_cnt, scale_factor_of_proposal);
        }
        final_states[chain_id] = interim_states;
        accepted_count[chain_id] = proposal_accepted_cnt;
    }

    timer.step("Overall Stop");
    Rcpp::NumericVector t(timer);
    out[0] = final_states;
    out[1] = accepted_count;
    out[2] = timer;
    return out;
}