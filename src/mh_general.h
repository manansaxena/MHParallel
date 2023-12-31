#include <cmath>
#include <iostream>
#include <chrono>
#include <tuple>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <vector>
#include <omp.h>
#include <functional>

using targetFunction = std::function<double(const Eigen::VectorXd&)>;

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
Eigen::VectorXd mh_step(RNG& rng, Eigen::VectorXd current_state, double* proposal_accepted_cnt, Eigen::MatrixXd covar_matrix, targetFunction target){
    Eigen::VectorXd random_vector = generate_random_vector(rng, current_state.size(), covar_matrix);
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


std::tuple<std::vector<Eigen::MatrixXd>,Eigen::VectorXd,int> metropolis_hastings (int num_chains,
                                Eigen::MatrixXd initial_states,
                                int num_steps,
                                int seed,
                                int n_cores,
                                targetFunction target,
                                double scale_factor_of_proposal = 1,
                                Eigen::MatrixXd covar = Eigen::MatrixXd::Identity(1,1)
                                )

{   
    Eigen::MatrixXd covar_matrix;
    if(covar.rows() == 1 && covar.cols() == 1){
        covar_matrix = scale_factor_of_proposal * Eigen::MatrixXd::Identity(initial_states.rows(), initial_states.rows());
    }
    else {
        covar_matrix = scale_factor_of_proposal * covar;
    }
    
    auto start = std::chrono::high_resolution_clock::now();

    if (n_cores > 0) {
      omp_set_num_threads(n_cores);
    } else {
      omp_set_num_threads(omp_get_max_threads());
    }

    // num_chains * (length of each state vector * num_steps) 
    std::vector<Eigen::MatrixXd> final_states(num_chains);
    Eigen::VectorXd accepted_count = Eigen::VectorXd::Zero(num_chains);
    
    #pragma omp parallel for 
    for(int chain_id = 0; chain_id < num_chains; chain_id++){
        boost::random::mt19937 rng(seed);
        rng.discard(omp_get_thread_num()*num_steps);

        Eigen::MatrixXd interim_states(initial_states.rows(),num_steps+1);
        interim_states.col(0) = initial_states.col(chain_id);
        double proposal_accepted_cnt = 0.0;
        for(int i=0;i<num_steps;i++){
            interim_states.col(i+1) = mh_step(rng, interim_states.col(i),&proposal_accepted_cnt, covar_matrix, target);
        }
        final_states[chain_id] = interim_states;
        accepted_count[chain_id] = proposal_accepted_cnt;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    return std::make_tuple(final_states, accepted_count, duration);
}