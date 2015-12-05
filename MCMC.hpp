#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "PatternHeaders.hpp"

#include <random>
#include <assert.h>

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void MCMCgenerator(size_t iter, Pattern& pattern, Matrix<size_t>& big_matrix)
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
		if (!pattern.avoid(big_matrix, r, c))
		{
			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;
			// and recalculate used structures if needed
			bool ok = pattern.revert(big_matrix, r, c);

			if (!ok) {
				assert(!"Matrix after reverting contains the pattern!");
				throw my_exception("Matrix after reverting contains the pattern!");
			}
		}
	}
}

#endif