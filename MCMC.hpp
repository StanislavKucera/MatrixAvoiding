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

#include <array>
#include <deque>
#include <queue>

#include <fstream>

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void MCMCgenerator(const size_t iter, Patterns& patterns, Matrix<size_t>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const size_t threads_count)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(1993);
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	size_t r, c;

	// initialize and run threads if parallel computation is requested
	if (threads_count > 1)
		patterns.construct_threads(big_matrix);

	// this is gonna be here for now
	std::vector<std::vector<Counter> > sizes;
	std::chrono::system_clock::time_point start, end;
	bool success;
	int last_perc = -1;

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

		start = std::chrono::system_clock::now();

		// test if the changed matrix still avoids the pattern
		if ((threads_count <= 1 && !patterns.avoid(big_matrix, sizes, r, c)) ||
			(threads_count > 1 && !patterns.parallel_avoid(threads_count, big_matrix, sizes, r, c)))
		{
			success = false;
			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);
			// and recalculate used structures if needed
			if (threads_count <= 1)
				patterns.revert(big_matrix, r, c);
			else
				patterns.parallel_revert(2, big_matrix, r, c);
		}
		else
			end = std::chrono::system_clock::now();

		end = std::chrono::system_clock::now();
		
		matrix_stats.add_data(i, ones, big_matrix);
		perf_stats.add_data(i, success, std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0, sizes);

		const int current_it = (i + 1) * 10 / iter;
		
		if (current_it > last_perc)
		{
			last_perc = current_it;

			if (last_perc == 0)
				std::cout << "Generating started" << std::endl;
			else
				std::cout << last_perc * 10 << " %" << std::endl;
		}
	}

	// destruct the threads if parallel computation was requested
	if (threads_count > 1)
		patterns.destruct_threads();
}

void parallel_avoid(Patterns& patterns, std::vector<std::vector<Counter> >& sizes, const Job& job, const std::atomic_bool& forced_end, std::atomic_bool& ret, std::atomic_bool& done,
	std::atomic_bool& ret_read, const std::atomic_bool& end, std::atomic_bool& main_job, std::condition_variable& my_cv, std::mutex& my_mtx, std::condition_variable& cv, std::mutex& mtx)
{
	while (!end)
	{
		{
			std::unique_lock<std::mutex> lck(my_mtx);

			while (done)
				my_cv.wait(lck);
		}

		if (job.avoid)
		{
			ret = patterns.avoid(sizes, job.r, job.c, forced_end);
			ret_read = false;
		}
		else
			patterns.revert(job.r, job.c);

		done = true;

		{
			std::unique_lock<std::mutex> lck2(mtx);
			main_job = true;
			cv.notify_one();
		}
	}
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void parallelMCMCgenerator(const size_t iter, Patterns& patterns, Matrix<size_t>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const size_t threads_count)
{
	// I wouldn't accomplish anothing using 0 workers
	if (threads_count == 0)
	{
		MCMCgenerator(iter, patterns, big_matrix, perf_stats, matrix_stats, 1);
		return;
	}
	
	/////////////////////////////////
	// workers and their variables //
	/////////////////////////////////
	std::vector<std::thread> threads(threads_count);
	// when parallel, each patterns has its own copy of a big_matrix, the copies may differ during the calculations by a little bit
	// accessed by its worker only - no need for a mutex
	patterns.set_matrix(big_matrix);
	// the same patterns for each thread
	// accessed by its worker only - no need for a mutex
	std::vector<Patterns> patterns_v(threads_count, patterns);
	// this forces the end of calculation if the main thread finds out there was a successful call of avoid in a calculation with lower id
	// accesed by its worker and the main thread - read by the worker; read and write by the main thread - thread safety by atomicity
	std::vector<std::atomic_bool> force_end(threads_count);
	// returned value of avoid function for each thread
	// accesed by its worker and the main thread - read by the worker; read by the main thread - thread safety by atomicity
	std::vector<std::atomic_bool> ret(threads_count);
	// indicator whether the calculation is done
	// accesed by its worker and the main thread - read and write by the worker; read and write by the main thread - thread safety by atomicity
	std::vector<std::atomic_bool> done(threads_count);
	// indicator whether the returned value was already read
	// accesed by its worker and the main thread - write by the worker; read and write by the main thread - thread safety by atomicity
	std::vector<std::atomic_bool> ret_read(threads_count);
	// job for a worker - (size_t r, size_t c, bool avoid): flip the bit at [r,c] and test avoid (avoid = true) or revert (avoid = false)
	// accesed by its worker and the main thread - read by the worker; write by the main thread - thread safety - changed only when worker sleeps
	std::vector<Job> jobs(threads_count);
	std::vector<std::condition_variable> cvs(threads_count);
	std::vector<std::mutex> mtxs(threads_count);
	std::condition_variable cv;
	std::mutex mtx;
	std::atomic_bool end;
	std::atomic_bool main_job(false);

	///////////////////////////
	// main thread variables //
	///////////////////////////
	// for each thread there is an pair saying which calculations need to be reverted
	// reverting[i] = a; means that in the i-th thread all calculations with id greater than or equal to a need to be reverted 
	std::vector<int> reverting(threads_count);
	// everytime there is a successful avoid taken, all the threads need to get its own matrix into a valid state - this is the queue of all changes
	std::vector<std::set<std::pair<int, Job> > > sync(threads_count);
	// for each thread there is a list of tasks done. Deque is used because I want to pop_front and pop_back
	// queue[thread][i] = <job, id, returned_val, reverted, synced>; means that the i-th not checked (by the main thread) calculation has id, changed the entry at position [job.r,job.c],
	// the avoid (if it was not revert) ended up with returned_val and it has been reverted already (that is needed either when avoid returns false or when calculation with smaller id succeeds),
	// if synced is true, the calculation was already checked by the main thread and all running avoid call know the result of this one, working with it
	std::vector<std::deque<Task> > queue(threads_count);
	// calculation ids - current means currently being checked (waited for) by the main thread, last is the lastly assigned id
	int current_id = 1, last_id = 0;
	// calculations ordered by id - use this to find the order in which I deal with the threads
	//std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t> >, std::greater<std::pair<size_t, size_t> > > priority;
	//std::set<std::pair<size_t, size_t> > priority;

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(1993);
	std::uniform_int_distribution<size_t> uni(0, big_matrix.getRow() - 1);

	// this is gonna be here for now
	std::vector<std::vector<std::vector<Counter> > > sizes(threads_count);
	//clock_t t;
	//bool success;
	int last_perc = -1;

	std::ofstream oFile("changes.txt"), ar("avoidrevert.txt");

	// to show used order in the statistics
	//perf_stats.set_order(patterns.get_order());

	// matrix statistics purposes
	size_t ones = big_matrix.getOnes();
	
	end = false;

	for (size_t i = 0; i != threads_count; ++i)
	{
		jobs[i].r = uni(rng);
		jobs[i].c = uni(rng);
		jobs[i].avoid = true;

		force_end[i] = false;
		done[i] = false;
		ret_read[i] = true;
		++last_id;
		queue[i].push_back(Task(jobs[i], last_id, last_id + 1, false, false, false));
		ar << i << ": avoid [" << jobs[i].r << "," << jobs[i].c << "]" << std::endl;
		threads[i] = std::thread(parallel_avoid, std::ref(patterns_v[i]), std::ref(sizes[i]), std::ref(jobs[i]), std::ref(force_end[i]), std::ref(ret[i]), std::ref(done[i]), std::ref(ret_read[i]),
			std::ref(end), std::ref(main_job), std::ref(cvs[i]), std::ref(mtxs[i]), std::ref(cv), std::ref(mtx));
	}

	size_t iterations = 0;

	// do as many iterations as requested, only counting those that propagate to the generated matrix
	while (iterations < iter)
	{
		//priority.clear();

		//for (size_t i = 0; i < threads_count; ++i)
		//{
		//	priority.push(std::make_pair(std::get<0>(queue[i].front()), i));
		//	priority.insert(std::make_pair(std::get<0>(queue[i].front()), i));
		//}

		{
			std::unique_lock<std::mutex> lck(mtx);

			// If there is no job for the main thread (every worker calculates its task), wait for some calculation to end
			while (!main_job)
				cv.wait(lck);

			main_job = false;
		}

		//while (!priority.empty())
		for (size_t index = 0; index != threads_count; ++index)
		//for (auto& prior : priority)
		{
			// there is always atleast one task in the worker queue - if not there is a bug in the implementation
			// exception can happen right before new avoid is called, but that is not at the beginning of the main loop
			if (queue[index].empty())
				assert(!"Queue is empty at the beginning of a cycle");

			//const size_t index = prior.second;
			//const size_t index = priority.top().second;
			//priority.pop();

			// there is only one job and it is still running
			if (queue[index].size() == 1 && !done[index])
				continue;

			// if there is a result of avoid call that hasn't been read yet
			if (!ret_read[index])
			{
				// ret_read is only set false after an avoid call and after the returned value is set - it always contains the right value
				queue[index].back().returned = ret[index];
				ret_read[index] = true;
			}

			// if this is the calculation the process is waiting for or there is nothing to revert and there wasn't synchronization running
			if ((queue[index].front().id == current_id || !force_end[index]) && queue[index].front().id != 0)
			{
				// pop from the queue all tasks that have been synced and already confirmed - the change was taken the to generated matrix
				while (!queue[index].empty() && queue[index].front().id < 0 && -queue[index].front().id <= current_id)
					queue[index].pop_front();

				// if after poping there are no tasks left, create another one
				if (queue[index].empty())
					goto sync_block;

				// there is only one job and it is still running
				if (queue[index].size() == 1 && !done[index])
					continue;

				// call of the avoid succeeded
				if (queue[index].front().returned)
				{
					const int id = queue[index].front().id;

					// if other workers haven't been notified that the avoid succeded (ids with negative ids are not avoid calls)
					if (!queue[index].front().synced && queue[index].front().id > 0)
					{
						queue[index].front().synced = true;
						// since other workers didn't know this calculation will successed they have been computing wrong things
						queue[index].front().next_id = last_id + 1;

						for (size_t j = 0; j < threads_count; ++j)
						{
							// the worker calculated (or is calculating) something it wouldn't have in serial case
							if (queue[j].back().id > id || -queue[j].back().id > id)
							{
								// stop the calculation
								force_end[j] = true;

								// and tell the worker to revert every computation with id greater than id
								if (reverting[j] > id || reverting[j] == 0)
									reverting[j] = id;
							}

							std::set<std::pair<int, Job> >::iterator it = sync[j].begin();

							// go through the queue of things the worker needs to synchronize
							for (; it != sync[j].end(); ++it)
								if (id < it->first)
									break;

							// if there is anything with higher id, delete it
							if (it != sync[j].end())
								sync[j].erase(it, sync[j].end());

							// don't have to synchronize myself
							if (j != index)
								// and add the most recently changed position to the list
								sync[j].emplace_hint(sync[j].end(), std::make_pair(id, queue[index].front().job));
						}
					}

					// the calculation I was waiting for is done
					if (queue[index].front().id == current_id)
					{
						// the next calculation I will wait for has its id equal to next_id
						current_id = queue[index].front().next_id;
						++iterations;

						// flip the bit in the generated matrix
						if (big_matrix.flip(queue[index].front().job.r, queue[index].front().job.c))
							++ones;
						else
							--ones;
						
						matrix_stats.add_data(iterations, ones, big_matrix);
						//oFile << queue[index].front().job.r << " " << queue[index].front().job.c << std::endl;
						
						// this was the last iteration of the generator
						if (iterations == iter)
							goto while_end;

						// the task was dealt with
						queue[index].pop_front();
					}
				}
				else
				{
					// the calculation I was waiting for is done
					if (queue[index].front().id == current_id)
					{
						// since the avoid call wasn't successful and all the calculations counted on it, the next on is the right one in the chain
						++current_id;
						++iterations;
						matrix_stats.add_data(iterations, ones, big_matrix);

						// this was the last iteration of the generator
						if (iterations == iter)
							goto while_end;

						// the task was dealt with
						queue[index].pop_front();
					}
				}
			}

			// cannot assign another task to a worker when the last one isn't done yet
			if (!done[index])
				continue;

			// if there is a result of avoid call that hasn't been read yet
			if (!ret_read[index])
			{
				if (queue[index].empty())
					assert("!Queue is empty while the last result wasn't read.");

				// ret_read is only set false after an avoid call and after the returned value is set - it always contains the right value
				queue[index].back().returned = ret[index];
				ret_read[index] = true;
			}

			// reverting changes of the matrix that wouldn't happen in the serial case
			if (force_end[index])
			{
				// no job needs to be reverted (I shouldn't be in a force_end state)
				if ((queue[index].back().id > 0 && queue[index].back().id <= reverting[index]) || (queue[index].back().id < 0 && -queue[index].back().id <= reverting[index]) || reverting[index] == 0)
					assert(!"The thread is in a forced_end state while having no incorrect calculations.");

				// even if this happens I won't influence the calculation but it still shouldn't happen
				if (!ret_read[index])
					assert(!"Somehow didn't notice the returned value of a previous avoid call.");

				// delete everything that has been reverted already
				while (!queue[index].empty() && (queue[index].back().id > reverting[index] || -queue[index].back().id > reverting[index]) && (!queue[index].back().returned || queue[index].back().reverted))
					queue[index].pop_back();

				// I've reverted everything I had to -> the reverting is over
				if (queue[index].empty() || (queue[index].back().id <= reverting[index] && queue[index].back().id > 0) || (-queue[index].back().id <= reverting[index] && queue[index].back().id < 0))
				{
					force_end[index] = false;
					reverting[index] = 0;
				}
				// there is still something to revert -> revert the last thing changed
				else
				{
					queue[index].back().reverted = true;
					jobs[index].r = queue[index].back().job.r;
					jobs[index].c = queue[index].back().job.c;
					jobs[index].avoid = false;
					//ar << index << ": revert [" << jobs[index].r << "," << jobs[index].c << "] - " << queue[index].back().id << std::endl;
					
					{
						std::unique_lock<std::mutex> lck(mtxs[index]);
						done[index] = false;
						// notify the worker
						cvs[index].notify_one();
					}

					continue;
				}
			}

		sync_block:

			// the last job was to synchronize and it is done
			if (!queue[index].empty() && queue[index].front().id == 0)
				queue[index].pop_front();

			// there are synchrozations to make - other worker computed something that propageted to the generated matrix so the worker needs to make the change in its internal matrix as well
			if (!sync[index].empty())
			{
				std::set<std::pair<int, Job> >::iterator it = sync[index].begin();

				// even if this happens I won't influence the calculation but it still shouldn't happen
				if (!ret_read[index])
					assert(!"Somehow didn't notice the returned value of a previous avoid call.");

				jobs[index].r = it->second.r;
				jobs[index].c = it->second.c;
				// I call revert since this change was already proved to be successful by another worker
				jobs[index].avoid = false;
				// put the task at the beginning of the list - POSSIBLY USELESS
				queue[index].push_front(Task(jobs[index], 0, it->first, false, false, false));

				// if the result isn't certain to be in the generated matrix I put the calculation (with negative id) to the end of the queue so it can still be reverted if needed
				if (queue[index].front().next_id > current_id)
					queue[index].push_back(Task(jobs[index], -it->first, 0, true, false, false));

				//ar << index << ": sync [" << it->second.r << "," << it->second.c << "] - " << 0 << std::endl;

				{
					std::unique_lock<std::mutex> lck(mtxs[index]);
					done[index] = false;
					cvs[index].notify_one();
				}

				sync[index].erase(it);
				continue;
			}
			
			// I somehow skipped reverting and I am going to create a completely new avoid task
			if (force_end[index] || reverting[index] != 0)
				assert(!"The thread is in a forced_end state while calling a new avoid.");

			// the result of the last avoid call wasn't read and I'm going to call another one
			if (!ret_read[index])
				assert(!"Somehow didn't notice the returned value of a previous avoid call.");

			jobs[index].r = uni(rng);
			jobs[index].c = uni(rng);
			jobs[index].avoid = true;
			++last_id;
			queue[index].push_back(Task(jobs[index], last_id, last_id + 1, false, false, false));
			//ar << index << ": avoid [" << jobs[index].r << "," << jobs[index].c << "] - " << last_id << std::endl;

			{
				std::unique_lock<std::mutex> lck(mtxs[index]);
				done[index] = false;
				cvs[index].notify_one();
			}
		}

	while_end:

		//perf_stats.add_data(i, success, t, sizes);

		// showing progress to the user - maybe let user choose if he want to see the progress?
		const int current_it = (iterations + 1) * 10 / iter;

		if (current_it > last_perc)
		{
			last_perc = current_it;

			if (last_perc == 0)
				std::cout << "Generating started" << std::endl;
			else
				std::cout << last_perc * 10 << " %" << std::endl;
		}
	}

	end = true;

	// threads clean up
	for (size_t i = 0; i != threads_count; ++i)
		force_end[i] = true;

	for (size_t i = 0; i != threads_count; ++i)
	{
		{
			std::unique_lock<std::mutex> lck(mtxs[i]);
			done[i] = false;
			cvs[i].notify_one();
		}

		threads[i].join();
	}
	// everything besides patterns are std containers or basic variables so they will destruct themselves
	// since the patterns is sent by l-value reference (not sure yet whether it is a good thing), I'm not the only owner
}

#endif