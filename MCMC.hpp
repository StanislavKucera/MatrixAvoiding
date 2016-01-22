#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "PatternHeaders.hpp"
#include "Statistics.hpp"

#include <random>
#include <assert.h>

#include <time.h>
#include <iostream>

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void MCMCgenerator(size_t iter, Pattern& pattern, Matrix<size_t>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(1993);
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// this is gonna be here for now
	std::vector<Counter> sizes;
	clock_t t;
	bool success;
	size_t last_perc = 0;

	// to show used order in the statistics
	perf_stats.set_order(pattern.get_order());

	// matrix statistics purposes
	size_t ones = big_matrix.getOnes();

	// go through iterations
	for (size_t i = 0; i < iter; ++i)
	{
		success = true;
		sizes.clear();

		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);

		t = clock();

		// test if the changed matrix still avoids the pattern
		if (!pattern.avoid(big_matrix, sizes, r, c))
		{
			success = false;

			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);
			// and recalculate used structures if needed
			bool ok = pattern.revert(big_matrix, sizes, r, c);

			if (!ok)
			{
				assert(!"Matrix after reverting contains the pattern!");
				throw my_exception("Matrix after reverting contains the pattern!");
			}
		}

		t = clock() - t;
		
		matrix_stats.add_data(i, ones, big_matrix);
		perf_stats.add_data(i, success, t, sizes);

		const size_t current_it = (i + 1) * 10 / iter;
		
		if (current_it > last_perc)
		{
			last_perc = current_it;
			std::cout << last_perc * 10 << " %\n";
		}
	}
}

#endif