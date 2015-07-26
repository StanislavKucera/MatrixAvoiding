#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "Matrix.hpp"
#include "AvoidanceTests.hpp"

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t, general_vector_pattern&, matrix<size_t>& big_matrix);

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t, general_set_pattern&, matrix<size_t>& big_matrix);

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t, walking_pattern&, matrix<size_t>& big_matrix);

#endif