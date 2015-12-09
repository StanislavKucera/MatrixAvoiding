#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "PatternHeaders.hpp"

#include <random>
#include <assert.h>

#include <time.h>
#include <iostream>

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void MCMCgenerator(size_t iter, Pattern& pattern, Matrix<size_t>& big_matrix)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// this is gonna be here for now
	std::vector<std::pair<std::pair<size_t, size_t>, size_t> > sizes;
	clock_t t;
	bool success;
	unsigned long long success_counter = 0, fail_counter = 0;
	unsigned long long success_levels = 0;
	std::vector<std::pair<std::pair<unsigned long long, unsigned long long>, unsigned long long> > success_sizes(big_matrix.getCol() + big_matrix.getRow()),
		fail_sizes(big_matrix.getCol() + big_matrix.getRow());
	unsigned long long success_time = 0, fail_time = 0;

	// go through iterations
	for (size_t i = 0; i < iter; ++i)
	{
		success = true;
		sizes.clear();

		r = uni(rng);
		c = uni(rng);
		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;

		t = clock();

		// test if the changed matrix still avoids the pattern
		if (!pattern.avoid(big_matrix, sizes, r, c))
		{
			success = false;

			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? 0 : 1;
			// and recalculate used structures if needed
			bool ok = pattern.revert(big_matrix, sizes, r, c);

			if (!ok) {
				assert(!"Matrix after reverting contains the pattern!");
				throw my_exception("Matrix after reverting contains the pattern!");
			}
		}

		t = clock() - t;
		
		if (success)
		{
			++success_counter;
			success_time += t;
			success_levels += sizes.size();

			for (size_t i = 0; i < sizes.size(); ++i) {
				success_sizes[i].first.first += sizes[i].first.first;
				success_sizes[i].first.second += sizes[i].first.second;
				success_sizes[i].second += sizes[i].second;
			}
		}
		else
		{
			++fail_counter;
			fail_time += t;

			for (size_t i = 0; i < sizes.size(); ++i) {
				fail_sizes[i].first.first += sizes[i].first.first;
				fail_sizes[i].first.second += sizes[i].first.second;
				fail_sizes[i].second += sizes[i].second;
			}
		}
	}

	// performance statistics:
	std::cout << "Successful calls of avoid (matrix avoids the pattern):\n";
	std::cout << "Count: " << success_counter << " - " << (double)success_counter / iter * 100 << "%\n";
	std::cout << "Average time per call: " << (double)success_time / success_counter / CLOCKS_PER_SEC << " sec\n";
	std::cout << "Average number of lines mapped: " << (double)success_levels / success_counter << "\n";
	// statistics for each level ...
	std::cout << "\nUnsuccessful calls of avoid (matrix doesn't avoid the pattern):\n";
	std::cout << "Count: " << fail_counter << " - " << (double)fail_counter / iter * 100 << "%\n";
	std::cout << "Average time per call: " << (double)fail_time / success_counter / CLOCKS_PER_SEC << " sec\n\n";
}

#endif