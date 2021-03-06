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
inline void MCMCgenerator(const int iter, Patterns& patterns, Matrix<bool>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const int /*threads_count*/, const int random_seed)
{
	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(random_seed == -1 ? rd() : random_seed);
	std::uniform_int_distribution<int> uni(0, big_matrix.getRow() - 1);

	// coordinates of changed element
	int r, c;

	// initialize and run threads if parallel computation is requested
	//if (threads_count > 1)
	//	patterns.construct_threads(big_matrix);

	// this is gonna be here for now
	std::vector<std::vector<Counter> > sizes;
	std::chrono::system_clock::time_point start, end;
	bool success;
	int last_perc = -1;

	// to show used order in the statistics
	perf_stats.set_order(patterns.get_order());

	// matrix statistics purposes
	int ones = big_matrix.getOnes();

	// go through iterations
	for (int i = 0; i < iter; ++i)
	{
		success = true;
		sizes.clear();

		r = uni(rng);
		c = uni(rng);

		// switch 0 and 1 entry of the element
		big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);

		start = std::chrono::system_clock::now();

		// test if the changed matrix still avoids the pattern
		if (/*(threads_count <= 1 &&*/ !patterns.avoid(big_matrix, r, c, sizes))// ||
			//(threads_count > 1 && !patterns.parallel_avoid(big_matrix, r, c, sizes, threads_count)))
		{
			success = false;
			// if not return to the previous matrix
			big_matrix.at(r, c) = big_matrix.at(r, c) ? (--ones, 0) : (++ones, 1);
			// and recalculate used structures if needed
		//	if (threads_count <= 1)
				patterns.revert(big_matrix, r, c);
		//	else
		//		patterns.parallel_revert(big_matrix, r, c, threads_count);
		}

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
	//if (threads_count > 1)
	//	patterns.destruct_threads();
}

struct worker_state
{
	worker_state(const Patterns& p) : patterns(p) {}
	worker_state(const worker_state& ws) : patterns(ws.patterns) {}

	// the same patterns for each thread
	// accessed by its worker only - no need for a mutex
	Patterns patterns;
	std::vector<std::vector<Counter> > sizes;
	// job for a worker - (int r, int c, bool avoid): flip the bit at [r,c] and test avoid (avoid = true) or revert (avoid = false)
	// accesed by its worker and the main thread - read by the worker; write by the main thread - thread safety - changed only when worker sleeps
	Job jobs;
	// this forces the end of calculation if the main thread finds out there was a successful call of avoid in a calculation with lower id
	// accesed by its worker and the main thread - read by the worker; read and write by the main thread - thread safety by atomicity
	std::atomic<bool> force_end;
	// returned value of avoid function for each thread
	// accesed by its worker and the main thread - read by the worker; read by the main thread - thread safety by atomicity
	std::atomic<bool> ret;
	// indicator whether the calculation is done
	// accesed by its worker and the main thread - read and write by the worker; read and write by the main thread - thread safety by atomicity
	std::atomic<bool> done;
	// condition variable allowing the worker to wait passively when there is nothing to do
	std::condition_variable cvs;
	// mutex for the condition variable to prevent the situation when the worker is being woken up while it goes to sleep
	std::mutex mtxs;
};

void parallel_avoid(worker_state& worker_state, const std::atomic<bool>& end, std::atomic<bool>& main_job, std::condition_variable& cv, std::mutex& mtx)
{
	while (!end)
	{
		{
			std::unique_lock<std::mutex> lck(worker_state.mtxs);

			while (worker_state.done)
				worker_state.cvs.wait(lck);
		}

		worker_state.ret = worker_state.patterns.lazy_avoid(worker_state.jobs.r, worker_state.jobs.c, worker_state.sizes, worker_state.force_end);
		worker_state.done = true;

		{
			std::unique_lock<std::mutex> lck2(mtx);
			main_job = true;
			cv.notify_one();
		}
	}
}

// Generates random-ish matrix of given size, which is avoiding given walking pattern. Uses iter iterations on markov chain.
inline void parallelMCMCgenerator(const int iter, Patterns& patterns, Matrix<bool>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const int threads_count, const int random_seed)
{
	// I wouldn't accomplish anothing using 0 workers
	if (threads_count == 0)
	{
		MCMCgenerator(iter, patterns, big_matrix, perf_stats, matrix_stats, 1, random_seed);
		return;
	}

	/////////////////////////////////
	// workers and their variables //
	/////////////////////////////////
	std::vector<std::thread> threads(threads_count);
	// when parallel, each patterns has its own copy of a big_matrix, the copies may differ during the calculations by a little bit
	// accessed by its worker only - no need for a mutex
	patterns.set_matrix(big_matrix);

	std::vector<worker_state> worker_states(threads_count, worker_state(patterns));
	std::condition_variable cv;
	std::mutex mtx;
	std::atomic<bool> end;
	std::atomic<bool> main_job(false);

	///////////////////////////
	// main thread variables //
	///////////////////////////
	// for each thread there is an pair saying which calculations need to be reverted
	// reverting[i] = a; means that in the i-th thread all calculations with id greater than or equal to a need to be reverted 
	std::vector<int> reverting(threads_count);
	// everytime there is a successful avoid taken, all the threads need to get its own matrix into a valid state - this is the queue of all changes
	std::vector<std::vector<std::pair<int, Job> > > sync(threads_count);
	// for each thread there is a list of tasks done. Deque is used because I want to pop_front and pop_back
	// queue[thread][i] = <job, id, returned_val, reverted, synced>; means that the i-th not checked (by the main thread) calculation has id, changed the entry at position [job.r,job.c],
	// the avoid (if it was not revert) ended up with returned_val and it has been reverted already (that is needed either when avoid returns false or when calculation with smaller id succeeds),
	// if synced is true, the calculation was already checked by the main thread and all running avoid call know the result of this one, working with it
	std::vector<std::deque<Task> > queue(threads_count);
	// calculation ids - current means currently being checked (waited for) by the main thread, last is the lastly assigned id
	int current_id = 1, last_id = 0;
	int last_perc = -1;
	int iterations = 0;
	// calculations ordered by id - use this to find the order in which I deal with the threads
	//std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t> >, std::greater<std::pair<size_t, size_t> > > priority;
	//std::set<std::pair<size_t, size_t> > priority;

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(random_seed == -1 ? rd() : random_seed);
	std::uniform_int_distribution<int> uni(0, big_matrix.getRow() - 1);

	//std::ofstream oFile("changes.txt"), ar("avoidrevert.txt");

	// to show used order in the statistics
	perf_stats.set_order(patterns.get_order(), true);

	// matrix statistics purposes
	int ones = big_matrix.getOnes();

	end = false;

	for (int i = 0; i != threads_count; ++i)
	{
		worker_states[i].jobs.r = uni(rng);
		worker_states[i].jobs.c = uni(rng);

		worker_states[i].force_end = false;
		worker_states[i].done = false;
		++last_id;
		queue[i].emplace_back(worker_states[i].jobs, last_id, last_id + 1, false, false);
		//ar << i << ": avoid [" << worker_states[i].jobs.r << "," << worker_states[i].jobs.c << "]" << std::endl;
		threads[i] = std::thread(parallel_avoid, std::ref(worker_states[i]), std::ref(end), std::ref(main_job), std::ref(cv), std::ref(mtx));
	}

	// do as many iterations as requested, only counting those that propagate to the generated matrix
	while (iterations < iter)
	{
		//priority.clear();

		//for (size_t i = 0; i < threads_count; ++i)
		//	priority.emplace(queue[i].front().id, i);

		{
			std::unique_lock<std::mutex> lck(mtx);

			// If there is no job for the main thread (every worker calculates its task), wait for some calculation to end
			while (!main_job)
				cv.wait(lck);

			main_job = false;
		}

		//while (!priority.empty())
		for (int index = 0; index != threads_count; ++index)
			//for (auto& prior : priority)
		{
			// there is always atleast one task in the worker queue - if not there is a bug in the implementation
			// exception can happen right before new avoid is called, but that is not at the beginning of the main loop
			assert(!queue[index].empty() || !"Queue is empty at the beginning of a cycle");

			//const size_t index = prior.second;
			//const size_t index = priority.top().second;
			//priority.pop();

			// if this is the calculation the process is waiting for or there is nothing to revert and there wasn't synchronization running
			if ((queue[index].front().id == current_id || !worker_states[index].force_end) && queue[index].front().id != 0)
			{
				// pop from the queue all tasks that have been synced and already confirmed - the change was taken the to generated matrix
				while (!queue[index].empty() && queue[index].front().id < 0 && -queue[index].front().id <= current_id)
					queue[index].pop_front();

				// if after poping there are no tasks left, create another one
				if (queue[index].empty())
					goto sync_block;

				// there is only one job and it is still running
				if (queue[index].size() == 1 && !worker_states[index].done)
					continue;

				// if there is a result of avoid call that hasn't been read yet
				if (worker_states[index].done)
					queue[index].back().returned = worker_states[index].ret;

				// call of the avoid succeeded
				if (queue[index].front().returned)
				{
					const int id = queue[index].front().id;

					// if other workers haven't been notified that the avoid succeded (ids with negative ids are not avoid calls)
					if (!queue[index].front().synced && queue[index].front().id > 0)
					{
						queue[index].front().synced = true;
						// since other workers didn't know this calculation will succeed they have been computing wrong things
						queue[index].front().next_id = last_id + 1;

						for (int j = 0; j < threads_count; ++j)
						{
							// the worker calculated (or is calculating) something it wouldn't have in serial case
							if (queue[j].back().id > id || -queue[j].back().id > id)
							{
								// stop the calculation
								worker_states[j].force_end = true;

								// and tell the worker to revert every computation with id greater than id
								if (reverting[j] > id || reverting[j] == 0)
									reverting[j] = id;
							}

							std::vector<std::pair<int, Job> >::iterator it = sync[j].begin();

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
								sync[j].emplace_back(id, queue[index].front().job);
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
			if (!worker_states[index].done)
				continue;

			// if there is a result of avoid call that hasn't been read yet
			if (!queue[index].empty())
				queue[index].back().returned = worker_states[index].ret;

			// reverting changes of the matrix that wouldn't happen in the serial case
			if (worker_states[index].force_end)
			{
				// no job needs to be reverted (I shouldn't be in a force_end state)
				assert((((queue[index].back().id > 0 && queue[index].back().id > reverting[index]) || (queue[index].back().id < 0 && -queue[index].back().id > reverting[index])) && reverting[index] != 0) ||
					!"The thread is in a forced_end state while having no incorrect calculations.");

				while (!queue[index].empty() && (queue[index].back().id > reverting[index] || -queue[index].back().id > reverting[index]))
				{
					if (queue[index].back().returned)
						worker_states[index].patterns.lazy_revert(queue[index].back().job.r, queue[index].back().job.c);

					queue[index].pop_back();
				}

				worker_states[index].force_end = false;
				reverting[index] = 0;
			}

		sync_block:

			if (!sync[index].empty())
			{
				for (const auto& s : sync[index])
				{
					worker_states[index].patterns.lazy_revert(s.second.r, s.second.c);

					// if the result isn't certain to be in the generated matrix I put the calculation (with negative id) to the end of the queue so it can still be reverted if needed
					if (s.first > current_id)
						queue[index].emplace_back(s.second, -s.first, 0, true, true);
				}

				sync[index].clear();
			}

			// I somehow skipped reverting and I am going to create a completely new avoid task
			assert((!worker_states[index].force_end && reverting[index] == 0) || !"The thread is in a forced_end state while calling a new avoid.");

			worker_states[index].jobs.r = uni(rng);
			worker_states[index].jobs.c = uni(rng);
			++last_id;
			queue[index].emplace_back(worker_states[index].jobs, last_id, last_id + 1, false, false);
			//ar << index << ": avoid [" << worker_states[index].jobs.r << "," << worker_states[index].jobs.c << "] - " << last_id << " (" << current_id << ")" << std::endl;

			{
				std::unique_lock<std::mutex> lck(worker_states[index].mtxs);
				worker_states[index].done = false;
				worker_states[index].cvs.notify_one();
			}
		}

	while_end:

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
	for (int i = 0; i != threads_count; ++i)
		worker_states[i].force_end = true;

	for (int i = 0; i != threads_count; ++i)
	{
		{
			std::unique_lock<std::mutex> lck(worker_states[i].mtxs);
			worker_states[i].done = false;
			worker_states[i].cvs.notify_one();
		}

		threads[i].join();
	}
	// everything besides patterns are std containers or basic variables so they will destruct themselves
	// since the patterns is sent by l-value reference (not sure yet whether it is a good thing), I'm not the only owner
}

/*
inline void get_sync(std::deque<std::pair<int, Job> >& sync, std::queue<std::pair<int, Job> >& new_syncs, std::mutex& new_syncs_mutex, std::atomic<bool>& synchronize, int revert = 0)
{
	// first throw away synchronizations that has been reverted
	while (!sync.empty() && sync.back().first > revert && revert != 0)
		sync.pop_back();

	std::pair<int, Job> job;

	while (true)
	{
		{
			std::unique_lock<std::mutex> lck(new_syncs_mutex);

			if (new_syncs.empty())
			{
				synchronize = false;
				break;
			}

			job = new_syncs.front();
			new_syncs.pop();
		}

		// the change I was supposed to synchronize was reverted
		if (revert != 0 && revert < job.first)
			continue;

		if (!sync.empty() && sync.back().first > job.first)
			sync.pop_back();

		sync.emplace_back(job);
	}
}

inline void make_iteration(const int my_index, std::deque<Task>& queue, Matrix<bool>& big_matrix, std::vector<std::atomic<bool>>& forced_ends, std::atomic<bool>& end,
	std::atomic<int>& current_id, std::atomic<int>& iterations, const int iter, std::atomic<int>& ones, Matrix_Statistics& matrix_stats, std::ostream& oFile, int& last_perc)
{
	while (!queue.empty() && queue.front().id == current_id)
	{
		// there is exactly one worker with a task having current_id - until current_id is changed the worker is the only one to access big_matrix, iterations and matrix_stats
		++iterations;

		// the change was successful
		if (queue.front().returned)
		{
			// flip the bit in the generated matrix
			if (big_matrix.flip(queue.front().job.r, queue.front().job.c))
				++ones;
			else
				--ones;
		}
		oFile << my_index << ": " << queue.front().job.r << " " << queue.front().job.c << " : " << queue.front().returned << " - " << queue.front().id << std::endl;

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

		matrix_stats.add_data(iterations, ones, big_matrix);

		// this was the last iteration of the generator
		if (iterations == iter)
		{
			end = true;

			for (size_t i = 0; i != forced_ends.size(); ++i)
				forced_ends[i] = true;

			// breaking here means current_id won't get changed and this is the last worker to access matrix_stats as well as big_matrix
			break;
		}

		// the next calculation I will wait for has its id equal to next_id
		current_id = queue.front().next_id;
		//std::cout << current_id << std::endl;
		// the task was dealt with
		queue.pop_front();
	}
}

void parallel_avoid2(const int my_index, Patterns patterns, Matrix<bool>& big_matrix, std::vector<std::vector<Counter> >& sizes, std::vector<std::atomic<bool>>& forced_ends,
	std::vector<std::atomic<bool>>& synchronize, std::atomic<bool>& end, std::atomic<int>& current_id, std::atomic<int>& last_id, std::atomic<int>& iterations, std::vector<int>& revertings,
	std::vector<std::queue<std::pair<int, Job> > >& syncs, const int N, const int iter, std::atomic<int>& ones, Matrix_Statistics& matrix_stats,
	std::vector<std::mutex>& syncs_mutexes, std::vector<std::mutex>& revertings_mutexes, std::mutex& last_id_mutex, int& last_perc, std::vector<std::atomic<int>>& last_job_id,
	std::ostream& oFile, std::ostream& ar, std::vector<std::atomic<bool>>& last_change_noted, std::vector<std::atomic<bool>>& reverts, const int random_seed)
{
	// everytime there is a successful avoid taken, all the threads need to get its own matrix into a valid state - this is the queue of all changes
	std::deque<std::pair<int, Job> > sync;
	// list of tasks done. Deque is used because I want to pop_front and pop_back
	// queue[thread][i] = <job, id, returned_val, reverted, synced>; means that the i-th not checked (by the main thread) calculation has id, changed the entry at position [job.r,job.c],
	// the avoid (if it was not revert) ended up with returned_val and it has been reverted already (that is needed either when avoid returns false or when calculation with smaller id succeeds),
	// if synced is true, the calculation was already checked by the main thread and all running avoid call know the result of this one, working with it
	std::deque<Task> queue; 
	int revert = 0;

	// random generator from uniform distribution [0, n-1]
	std::random_device rd;
	std::mt19937 rng(random_seed);//1993);//random_seed == -1 ? rd() : random_seed + my_index);
	std::uniform_int_distribution<int> uni(0, N - 1);

	while (!end)
	{
		// synchronized changes that propagated to the big matrix
		while (!queue.empty() && ((queue.front().id < 0 && queue.front().id >= -current_id)))// || (queue.front().id < current_id && queue.front().id > 0 && !forced_ends[my_index])))
			queue.pop_front();

		// this should not happen but it does all the time. Current_id is bigger than my first job while Im not supposed to revert it
		if (!queue.empty() && queue.front().id < current_id && queue.front().id > 0 && !reverts[my_index])
			revert = revert;

		// there is something to revert
		while (reverts[my_index])
		{
			// check the first job result if its the one the process is waiting for
			if (!queue.empty() && queue.front().id == current_id)
				make_iteration(my_index, queue, big_matrix, forced_ends, end, current_id, iterations, iter, ones, matrix_stats, oFile, last_perc);

			if (end)
				break;

			{
				std::unique_lock<std::mutex> lck(revertings_mutexes[my_index]);
				revert = revertings[my_index];
			}

			// revert all tasks with higher id than revert
			while (!end && !queue.empty() && (queue.back().id > revert || -queue.back().id > revert))
			{
				// if the change was successful
				if (queue.back().returned)
					patterns.revert(queue.back().job.r, queue.back().job.c);

				ar << my_index << ": revert [" << queue.back().job.r << "," << queue.back().job.c << "]" << std::endl;
				queue.pop_back();
			}

			// there is something to put into the synchronization list
			if (synchronize[my_index])
				get_sync(sync, syncs[my_index], syncs_mutexes[my_index], synchronize[my_index], revert);

			{
				std::unique_lock<std::mutex> lck(revertings_mutexes[my_index]);

				// it wasn't changed what to revert since the last time I checked
				if (revert == revertings[my_index])
				{
					// reverting is over
					revertings[my_index] = 0;
					reverts[my_index] = false;
					forced_ends[my_index] = false;
				}
			}
		}

		// there is something to put into the synchronization list
		if (synchronize[my_index])
			get_sync(sync, syncs[my_index], syncs_mutexes[my_index], synchronize[my_index]);

		// there is something in the synchronization list
		while (!reverts[my_index] && !sync.empty())
		{
			// check the first job result if its the one the process is waiting for
			if (!queue.empty() && queue.front().id == current_id)
				make_iteration(my_index, queue, big_matrix, forced_ends, end, current_id, iterations, iter, ones, matrix_stats, oFile, last_perc);

			if (end)
				break;

			auto job = sync.front();

			if (!queue.empty() && job.first < queue.back().id)
				revert = revert;

			sync.pop_front();
			ar << my_index << ": sync [" << job.second.r << "," << job.second.c << "]" << std::endl;
			patterns.revert(job.second.r, job.second.c);
			last_job_id[my_index] = job.first;

			// the synchronization is not certain to be taken to the big matrix - I add it to tasks done so it can be reverted later if needed
			if (job.first > current_id)
				queue.emplace_back(job.second, -job.first, 0, true, true);
		}
		
		// check the first job result if its the one the process is waiting for
		if (!queue.empty() && queue.front().id == current_id)
			make_iteration(my_index, queue, big_matrix, forced_ends, end, current_id, iterations, iter, ones, matrix_stats, oFile, last_perc);

		// generator ends or there is something to revert or there is something to add to the synchronization list
		if (end || reverts[my_index] || synchronize[my_index])
			continue;

		// create a new job
		int r = uni(rng);
		int c = uni(rng);
		Job new_job(r, c);
		// id of the newly created job
		int my_id = 0;

		{
			std::unique_lock<std::mutex> lck(last_id_mutex);

			// generator ends or there is something to revert or there is something to add to the synchronization list
			if (end || reverts[my_index] || synchronize[my_index] || !last_change_noted[my_index])
			{
				last_change_noted[my_index] = true;
				continue;
			}
			else
				my_id = ++last_id;
		}

		// add the new task to the list
		queue.emplace_back(new_job, my_id, my_id + 1, false, false);
		last_job_id[my_index] = my_id;
		queue.back().returned = patterns.avoid(r, c, sizes, forced_ends[my_index]);
		ar << my_index << ": avoid [" << r << "," << c << "] : " << queue.back().returned << " - " << my_id << " - " << current_id << std::endl;
		
		//if (queue.size() > 1000)
		//{
		//	end = true;

		//	for (size_t i = 0; i != forced_ends.size(); ++i)
		//		forced_ends[i] = true;

			// breaking here means current_id won't get changed and this is the last worker to access matrix_stats as well as big_matrix
		//	break;
		//}

		// call of the avoid succeeded
		if (queue.back().returned)
		{
			// if I was asked to revert something in the meantime
			if (end || reverts[my_index])
			{ 
				std::unique_lock<std::mutex> lck(revertings_mutexes[my_index]);

				// and the reverting result in the task having no impact
				if (end || (revertings[my_index] < my_id && revertings[my_index] != 0))
					continue;
			}

			for (size_t j = 0; j < forced_ends.size(); ++j)
			{
				// no need to synchronize/revert myself
				if ((int)j == my_index)
					continue;

				{
					std::unique_lock<std::mutex> lck(syncs_mutexes[j]);
					// synchronize the worker
					synchronize[j] = true;
					syncs[j].emplace(my_id, new_job);
				}

				{
					std::unique_lock<std::mutex> lck(revertings_mutexes[j]);

					// if the worker doesn't revert anything or it reverts something with higher id and it computed/computes a task with higher id (don't want to force end a computation with smaller id)
					if (revertings[j] > my_id || revertings[j] == 0)// && my_id <= last_job_id[j])
					{
						if (my_id <= last_job_id[j])
							// stop the calculation
							forced_ends[j] = true;

						reverts[j] = true;
						// and tell the worker to revert every computation with id greater than id
						revertings[j] = my_id;
					}
				}
			}

			{
				std::unique_lock<std::mutex> lck(last_id_mutex);
				// since other workers didn't know this calculation will succeed they have been computing wrong things and I need to skip those results
				for (size_t i = 0; i != last_change_noted.size(); ++i)
					if (i != (size_t)my_index)
						last_change_noted[i] = false;

				queue.back().next_id = last_id + 1;
			}
		}
	}
}

inline void parallelMCMCgenerator2(const int iter, Patterns& patterns, Matrix<bool>& big_matrix, Performance_Statistics& perf_stats, Matrix_Statistics& matrix_stats, const int threads_count, const int random_seed)
{
	// I wouldn't accomplish anothing using 0 workers
	if (threads_count == 0)
	{
		MCMCgenerator(iter, patterns, big_matrix, perf_stats, matrix_stats, 1, random_seed);
		return;
	}

	/////////////////////////////////
	// workers and their variables //
	/////////////////////////////////
	std::vector<std::thread> threads(threads_count - 1);
	// when parallel, each patterns has its own copy of a big_matrix, the copies may differ during the calculations by a little bit
	// accessed by its worker only - no need for a mutex
	patterns.set_matrix(big_matrix);
	// this forces the end of calculation if the main thread finds out there was a successful call of avoid in a calculation with lower id
	// accesed by its worker and the main thread - read by the worker; read and write by the main thread - thread safety by atomicity
	std::vector<std::atomic<bool> > force_end(threads_count);
	std::vector<std::atomic<bool> > reverts(threads_count);
	std::vector<std::queue<std::pair<int, Job> > > syncs(threads_count);
	std::vector<std::atomic<int> > last_job_id(threads_count);
	std::vector<std::atomic<bool> > synchronize(threads_count);
	std::vector<std::mutex> syncs_mutexes(threads_count);
	std::vector<std::atomic<bool> > last_change_noted(threads_count);

	///////////////////////////
	// main thread variables //
	///////////////////////////
	// for each thread there is a pair saying which calculations need to be reverted
	// reverting[i] = a; means that in the i-th thread all calculations with id greater than or equal to a need to be reverted 
	std::vector<int> reverting(threads_count);
	std::vector<std::mutex> reverting_mutexes(threads_count);
	// calculation ids - current means currently being checked (waited for) by the main thread, last is the lastly assigned id
	std::atomic<int> current_id(1), last_id(0);
	std::mutex last_id_mutex;
	std::atomic<int> iterations(0);
	// calculations ordered by id - use this to find the order in which I deal with the threads
	//std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t> >, std::greater<std::pair<size_t, size_t> > > priority;
	//std::set<std::pair<size_t, size_t> > priority;
	
	// this is gonna be here for now
	std::vector<std::vector<std::vector<Counter> > > sizes(threads_count);
	int last_perc = -1;

	std::ofstream oFile("changes.txt"), ar("avoidrevert.txt");

	// to show used order in the statistics
	//perf_stats.set_order(patterns.get_order());

	// matrix statistics purposes
	std::atomic<int> ones(big_matrix.getOnes());

	std::atomic<bool> end(false);

	const int size = big_matrix.getCol();

	for (int i = 1; i != threads_count; ++i)
	{
		force_end[i] = false;
		last_change_noted[i] = true;
		threads[i - 1] = std::thread(parallel_avoid2, i, patterns, std::ref(big_matrix), std::ref(sizes[i]), std::ref(force_end), std::ref(synchronize), std::ref(end),
			std::ref(current_id), std::ref(last_id), std::ref(iterations), std::ref(reverting), std::ref(syncs), size, iter, std::ref(ones), std::ref(matrix_stats),
			std::ref(syncs_mutexes), std::ref(reverting_mutexes), std::ref(last_id_mutex), std::ref(last_perc), std::ref(last_job_id), std::ref(oFile), std::ref(ar),
			std::ref(last_change_noted), std::ref(reverts), random_seed);
	}

	force_end[0] = false;
	last_change_noted[0] = true;
	parallel_avoid2(0, patterns, big_matrix, sizes[0], force_end, synchronize, end, current_id, last_id, iterations, reverting, syncs,
		size, iter, ones, matrix_stats, syncs_mutexes, reverting_mutexes, last_id_mutex, last_perc, last_job_id, oFile, ar, last_change_noted, reverts, random_seed);

	for (int i = 1; i != threads_count; ++i)
		threads[i - 1].join();

	oFile.close();
	ar.close();
}
*/

#endif