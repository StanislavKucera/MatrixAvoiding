#ifndef MCMC_cpp_
#define MCMC_cpp_

#include "MCMC.hpp"
#include <random>

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t iter, general_vector_pattern& pattern, matrix<size_t>& big_matrix)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// go through iterations
	for (size_t i = 0; i < iter; ++i)
	{
		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;

		// test if the changed matrix still avoids the pattern, it obviously holds, when one-entry changes to zero-entry
		if (big_matrix.at(r, c) && !pattern.avoid(r, c, big_matrix))
			// if not return to the previous step
			big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;
	}
}

// Generates random-ish matrix of given size, which is avoiding given general pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t iter, general_set_pattern& pattern, matrix<size_t>& big_matrix)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// go through iterations
	for (size_t i = 0; i < iter; ++i)
	{
		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;

		// test if the changed matrix still avoids the pattern, it obviously holds, when one-entry changes to zero-entry
		if (big_matrix.at(r, c) && !pattern.avoid(big_matrix))
			// if not return to the previous step
			big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;
	}
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
void MCMCgenerator(size_t iter, walking_pattern& pattern, matrix<size_t>& big_matrix)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;    
	std::mt19937 rng(rd());   
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1); 

	// coordinates of changed element
	size_t r, c;

	// go through iterations
	for (size_t i = 0; i < iter; ++i)
	{
		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;	

		// test if the changed matrix still avoids the pattern
		if (!pattern.avoid(r, c, big_matrix))			
		{
			// if not return to the previous step
			big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;
			// and recalculate used structures
			pattern.revert(r, c, big_matrix);
		}
	}
}

#endif