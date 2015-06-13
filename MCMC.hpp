#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "Matrix.hpp"
#include "AvoidanceTests.hpp"

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t, const size_t, general_pattern&);
// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t, const size_t, walking_pattern&);

#endif