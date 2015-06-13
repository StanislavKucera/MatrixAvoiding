#ifndef MCMC_cpp_
#define MCMC_cpp_

#include "MCMC.hpp"
#include <random>

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t n, const size_t iter, general_pattern& test)
{
	matrix<size_t> N(n, n, 0);		// generatated matrix

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> uni(0, n - 1);

	size_t r, c;	// coords of changed element
	for (size_t i = 0; i < iter; i++)
	{
		r = uni(rng);
		c = uni(rng);
		N.at(r, c) = N.at(r, c) ? 0 : 1;	// switch 0 and 1 entry of the element
		if (!test.avoid(N))					// avoid returns true, if pattern avoids the matrix
			N.at(r, c) = N.at(r, c) ? 0 : 1;
	}
	return N;
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t n, const size_t iter, walking_pattern& test)
{
	matrix<size_t> N(n, n, 0);		// generatated matrix

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;    
	std::mt19937 rng(rd());   
	std::uniform_int_distribution<size_t> uni(0, n - 1); 

	size_t r, c;	// coords of changed element
	for (size_t i = 0; i < iter; i++)
	{
		r = uni(rng);
		c = uni(rng);
		N.at(r, c) = N.at(r, c) ? 0 : 1;	// switch 0 and 1 entry of the element
		if (!test.avoid(r, c, N))			// avoid returns true, if pattern avoids the matrix
		{
			N.at(r, c) = N.at(r, c) ? 0 : 1;
			test.revert(r, c, N);
		}
	}
	return N;
}

#endif