#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "Matrix.hpp"
#include "AvoidanceTests.hpp"

// Generates random-ish matrix of given size, which is avoiding given pattern. Uses iter iterations on markov chain.
matrix<int> MCMCgenerator(const size_t n, const size_t iter, grandfather_pattern* test);

#endif