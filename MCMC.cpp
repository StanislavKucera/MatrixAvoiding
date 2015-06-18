#ifndef MCMC_cpp_
#define MCMC_cpp_

#include "MCMC.hpp"
#include <random>

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t n, const size_t iter, general_pattern& test)
{
	// generatated matrix
	matrix<size_t> N(n, n, 0);		

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> uni(0, n - 1);

	// coordinates of changed element
	size_t r, c;

	// go through iterations
	for (size_t i = 0; i < iter; i++)
	{
		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		N.at(r, c) = N.at(r, c) ? 0 : 1;	

		// test if the changed matrix still avoids the pattern
		if (!test.avoid(N))					
			// if not return to the previous step
			N.at(r, c) = N.at(r, c) ? 0 : 1;
	}

	// return the resulting matrix
	return std::move(N);
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
matrix<size_t> MCMCgenerator(const size_t n, const size_t iter, walking_pattern& test)
{
	// generatated matrix
	matrix<size_t> N(n, n, 0);		

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;    
	std::mt19937 rng(rd());   
	std::uniform_int_distribution<size_t> uni(0, n - 1); 

	// coordinates of changed element
	size_t r, c;

	// go through iterations
	for (size_t i = 0; i < iter; i++)
	{
		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		N.at(r, c) = N.at(r, c) ? 0 : 1;	

		// test if the changed matrix still avoids the pattern
		if (!test.avoid(r, c, N))			
		{
			// if not return to the previous step
			N.at(r, c) = N.at(r, c) ? 0 : 1;
			// and recalculate used structures
			test.revert(r, c, N);
		}
	}

	// return the resulting matrix
	return std::move(N);
}

#endif