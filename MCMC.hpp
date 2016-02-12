#ifndef MCMC_hpp_
#define MCMC_hpp_

#include "PatternHeaders.hpp"
#include "Statistics.hpp"

#include <random>
#include <assert.h>

#include <time.h>
#include <iostream>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <tuple>
#include <deque>
#include <queue>

#include <fstream>

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void MCMCgenerator(const size_t iter, Patterns& patterns, Matrix<size_t>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(1993);
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// this is gonna be here for now
	std::vector<std::vector<Counter> > sizes;
	clock_t t;
	bool success;
	size_t last_perc = 0;

	// to show used order in the statistics
	perf_stats.set_order(patterns.get_order());

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
		if (!patterns.avoid(big_matrix, sizes, r, c))
		{
			success = false;
			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);
			// and recalculate used structures if needed
			patterns.revert(big_matrix, r, c);
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

void parallel_avoid(Patterns& patterns, std::vector<std::vector<Counter> >& sizes, const size_t r, const size_t c, const std::vector<bool>::reference forced_end,
	std::vector<bool>::reference ret, std::vector<bool>::reference done, std::vector<bool>::reference ret_read)
{
	ret = patterns.avoid(sizes, r, c, forced_end);
	ret_read = false;
	done = true;
}

void parallel_revert(Patterns& patterns, const size_t r, const size_t c, std::vector<bool>::reference done, const Matrix<size_t>& mat)
{
	patterns.revert(r, c, mat);
	done = true;
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void parallelMCMCgenerator(const size_t iter, Patterns& patterns, Matrix<size_t>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const size_t threads_count)
{
	if (threads_count == 0)
	{
		MCMCgenerator(iter, patterns, big_matrix, perf_stats, matrix_stats);
		return;
	}

	// calculation ids - current means currently being checked (waited for) by the main thread, last is the lastly assigned id
	size_t current = 1, last = 0;
	// calculations ordered by id - use this to find the order in which I deal with the threads
	std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t> >, std::greater<std::pair<size_t, size_t> > > priority;

	std::vector<std::thread> threads(threads_count);
	// when parallel, each patterns has its own copy of a big_matrix, the copies may differ during the calculations by a little bit
	patterns.set_matrix(big_matrix);
	// the same patterns for each thread
	std::vector<Patterns> patterns_v(threads_count, patterns);
	// this forces the end of calculation if the main thread finds out there was a successful call of avoid in a calculation with lower id
	std::vector<bool> force_end(threads_count);
	// return value of avoid function for each thread and indicator whether the calculation is done. Those are the only variables the threads can modify (besides the patterns itself)
	std::vector<bool> ret(threads_count), state(threads_count);
	// indicator whether the returned value was already read. Wouldn't be nice if I occasionally changed the returned value of a random call.
	std::vector<bool> ret_read(threads_count);
	// for each thread there is an pair saying which calculations need to be reverted
	// reverting[i] = a; means that in the i-th thread all calculations with id greater than or equal to a need to be reverted 
	std::vector<size_t> reverting(threads_count);
	// everytime there is a successful avoid taken, all the threads need to get its own matrix into a valid state - this is the queue of all changes
	std::vector<std::queue<std::pair<size_t, size_t> > > sync(threads_count);
	// for each thread there is a list of tasks done. Deque is used because I want to pop_front and pop_back
	// queue[thread][i] = <id, r, c, return_val, reverted>; means that the i-th not checked (by the main thread) calculation has id, changed the entry at position [r,c],
	// the avoid ended up with return_val and it has been reverted already (that is needed either when avoid returns false or when calculation with smaller id succeeds)
	std::vector<std::deque<std::tuple<size_t, size_t, size_t, bool, bool, bool> > > queue(threads_count);

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(1993);
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// this is gonna be here for now
	std::vector<std::vector<Counter> > sizes;
	clock_t t;
	bool success;
	size_t last_perc = 0;

	std::ofstream oFile("changes.txt"), ar("avoidrevert.txt");

	// to show used order in the statistics
	//perf_stats.set_order(patterns.get_order());

	// matrix statistics purposes
	size_t ones = big_matrix.getOnes();
	
	for (int i = 0; i != threads_count; ++i)
	{
		r = uni(rng);
		c = uni(rng);

		force_end[i] = false;
		state[i] = false;
		queue[i].push_back(std::make_tuple(++last, r, c, false, false, true));
		ar << i << ": avoid [" << r << "," << c << "]" << std::endl;
		threads[i] = std::thread(parallel_avoid, std::ref(patterns_v[i]), std::ref(sizes), r, c, force_end[i], ret[i], state[i], ret_read[i]);
	}

	// go through iterations
	size_t iterations = 0;
	while (iterations < iter)
	{
		//success = true;
		//sizes.clear();

		for (size_t i = 0; i < threads_count; ++i)
		{
			priority.push(std::make_pair(std::get<0>(queue[i].front()), i));
		}

		while (!priority.empty())
		{
			const size_t index = priority.top().second;
			priority.pop();

			// special value indicating that there is a calculation with a smaller id using the same thread
			if (index == threads_count)
				continue;

			// required calculation isn't done yet
			if (queue[index].size() == 1 && !state[index])
				continue;

			// special id value indicating that there was a synchronization running
			if (std::get<0>(queue[index].front()) == 0)
			{
				if (queue[index].size() > 1)
					assert(!"After sync there are more than one calculations assigned to a thread.");

				threads[index].join();
				goto sync_block;
			}
			
			// just in case front = back
			if (queue[index].size() == 1 && !ret_read[index])
			{
				std::get<3>(queue[index].back()) = ret[index];
				ret_read[index] = true;
			}

			// if other threads are waiting for this one or there is nothing to revert
			if (std::get<0>(queue[index].front()) == current || !force_end[index])
			{
				// call of the avoid succeeded
				if (std::get<3>(queue[index].front()))
				{
					const size_t id = std::get<0>(queue[index].front());

					/*// all the calculations with higher id expected this one to fail, so I need to revert them
					if (std::get<5>(queue[index].front()))
					{
						std::get<5>(queue[index].front()) = false;

						for (size_t j = 0; j < threads_count; ++j)
						{
							if (std::get<0>(queue[j].back()) > id)
							{
								// this only affects calls of avoid, not revert
								force_end[j] = true;

								if (reverting[j] > id || reverting[j] == 0)
									reverting[j] = id;
							}
						}
					}*/

					// the calculation I was waiting for is done
					if (std::get<0>(queue[index].front()) == current)
					{
						// all the calculations including the last one are useless
						current = last + 1;
						++iterations;
						big_matrix.flip(std::get<1>(queue[index].front()), std::get<2>(queue[index].front()));
						oFile << std::get<1>(queue[index].front()) << " " << std::get<2>(queue[index].front()) << std::endl;
						

						if (iterations == iter)
							goto while_end;

						for (size_t j = 0; j != threads_count; ++j)
							if (j != index)
								sync[j].push(std::make_pair(std::get<1>(queue[index].front()), std::get<2>(queue[index].front())));

						queue[index].pop_front();

						for (size_t j = 0; j < threads_count; ++j)
						{
							if (j == index && queue[index].empty())
								continue;

							if (std::get<0>(queue[j].back()) > id)
							{
								// this only affects calls of avoid, not revert
								force_end[j] = true;

								if (reverting[j] > id || reverting[j] == 0)
									reverting[j] = id;
							}
						}
					}
				}
				else
				{
					if (std::get<4>(queue[index].front()))
					{
						// the calculation I was waiting for is done
						if (std::get<0>(queue[index].front()) == current)
						{
							++current;
							++iterations;
							if (iterations == iter)
								goto while_end;

							queue[index].pop_front();
						}
					}
					else if (!state[index])
						continue;
					else
					{
						threads[index].join();

						std::get<4>(queue[index].front()) = true;
						state[index] = false;
						ar << index << ": revert [" << std::get<1>(queue[index].front()) << "," << std::get<2>(queue[index].front()) << "] - " << std::get<0>(queue[index].front()) << std::endl;
						threads[index] = std::thread(parallel_revert, std::ref(patterns_v[index]), std::get<1>(queue[index].front()), std::get<2>(queue[index].front()), state[index], std::ref(big_matrix));
						continue;
					}
				}
			}

			if (!state[index])
				continue;

			// it should be done
			threads[index].join();
			
			if (state[index] && !ret_read[index])
			{
				std::get<3>(queue[index].back()) = ret[index];
				ret_read[index] = true;
			}

			// reverting
			if (force_end[index])
			{
				// this should not happen, but just in case
				if (std::get<0>(queue[index].back()) <= reverting[index] || reverting[index] == 0)
					assert(!"The thread is in a forced_end state while having no incorrect calculations.");

				// it was already reverted
				if (/*!std::get<3>(queue[index].back()) ||*/ std::get<4>(queue[index].back()))
				{
					queue[index].pop_back();

					while (!queue[index].empty() && std::get<0>(queue[index].back()) > reverting[index] && (/*!std::get<3>(queue[index].back()) ||*/ std::get<4>(queue[index].back())))
						queue[index].pop_back();

					// I've reverted everything I had to
					if (queue[index].empty() || std::get<0>(queue[index].back()) <= reverting[index])
					{
						force_end[index] = false;
						reverting[index] = 0;
					}
					else
					{
						/*if (!std::get<3>(queue[index].back()))
							std::cout << "Something is bad ";*/
						std::get<4>(queue[index].back()) = true;
						state[index] = false;
						ar << index << ": revert [" << std::get<1>(queue[index].back()) << "," << std::get<2>(queue[index].back()) << "] - " << std::get<0>(queue[index].back()) << std::endl;
						threads[index] = std::thread(parallel_revert, std::ref(patterns_v[index]), std::get<1>(queue[index].back()), std::get<2>(queue[index].back()), state[index], std::ref(big_matrix));
						break;
					}
				}
				else
				{
					/*if (!std::get<3>(queue[index].back()))
						std::cout << "Something is bad ";*/
					std::get<4>(queue[index].back()) = true;
					state[index] = false;
					ar << index << ": revert [" << std::get<1>(queue[index].back()) << "," << std::get<2>(queue[index].back()) << "] - " << std::get<0>(queue[index].back()) << std::endl;
					threads[index] = std::thread(parallel_revert, std::ref(patterns_v[index]), std::get<1>(queue[index].back()), std::get<2>(queue[index].back()), state[index], std::ref(big_matrix));
					break;
				}
			}

		sync_block:

			if (!queue[index].empty() && std::get<0>(queue[index].front()) == 0)
				queue[index].pop_front();

			if (!sync[index].empty())
			{
				if (!queue[index].empty())
					assert(!"Trying to sync while there is still something else to do.");

				std::pair<size_t, size_t> pair = sync[index].front();
				sync[index].pop();
				
				state[index] = false;
				queue[index].push_front(std::make_tuple(0, pair.first, pair.second, false, false, true));
				ar << index << ": sync [" << pair.first << "," << pair.second << "] - " << 0 << std::endl;
				threads[index] = std::thread(parallel_revert, std::ref(patterns_v[index]), pair.first, pair.second, state[index], std::ref(big_matrix));
				continue;
			}
			//if (queue[index].empty())
			//	if (patterns_v[index].check_matrix(big_matrix))

			r = uni(rng);
			c = uni(rng);

			//t = clock();

			if (force_end[index] || reverting[index] != 0)
				assert(!"The thread is in a forced_end state while calling a new avoid.");
			state[index] = false;
			queue[index].push_back(std::make_tuple(++last, r, c, false, false, true));
			ar << index << ": avoid [" << r << "," << c << "] - " << last << std::endl;
			threads[index] = std::thread(parallel_avoid, std::ref(patterns_v[index]), std::ref(sizes), r, c, force_end[index], ret[index], state[index], ret_read[index]);

			//t = clock() - t;

		}

	while_end:

		//matrix_stats.add_data(i, ones, big_matrix);
		//perf_stats.add_data(i, success, t, sizes);

		const size_t current_it = (iterations + 1) * 10 / iter;

		if (current_it > last_perc)
		{
			last_perc = current_it;
			std::cout << last_perc * 10 << " %\n";
		}
	}

	// need to clean up
	for (size_t i = 0; i != threads_count; ++i)
		force_end[i] = true;
	for (size_t i = 0; i != threads_count; ++i)
		threads[i].join();
	// everything besides patterns are std containers or basic variables so they will destruct themselves
	// since the patterns is send by l-value reference (not sure yet whether it is a good thing), I'm not the only owner
}

#endif