#ifndef GeneralPatternFunctions_hpp_
#define GeneralPatternFunctions_hpp_

#include "PatternHeaders.hpp"

#include <queue>
#include <algorithm>

template<typename T>
General_pattern<T>::General_pattern(const Matrix<size_t>& pattern, const size_t threads_count, const Order order, const Map map_approach, std::vector<size_t>&& custom_order)
	: row_(pattern.getRow()),
	col_(pattern.getCol()), 
	lines_(row_ + col_),
	order_(row_ + col_),
	what_to_remember_(row_ + col_ + 1),
	parallel_bound_indices_(row_ + col_),
	extending_order_(row_ + col_),
	map_index_(row_ + col_ + 1),
	building_tree_(2),
	steps_(row_ + col_),
	empty_lines_(0),
	map_approach_(map_approach), 
	threads_(threads_count),
	qs_(threads_count),
	mutexes_(threads_count),
	cvs_(threads_count),
	mtxs_(threads_count),
	sleeps_(threads_count),
	done_(true),
	end_(false)
{
	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			// reading lines of pattern and storing them memory efficiently
			if (pattern.at(i, j))
			{
				lines_[i] |= 1 << j;
				lines_[j + row_] |= 1 << i;
			}
		}

		// if the row has no one-entries, I don't have to map it and can make less steps
		if (lines_[i] == 0)
		{
			empty_lines_ |= 1 << i;
			--steps_;
		}
	}

	// if the column has no one-entries, I don't have to map it and can make less steps
	for (size_t j = row_; j < row_ + col_; j++)
	{
		if (lines_[j] == 0)
		{
			empty_lines_ |= 1 << j;
			--steps_;
		}
	}

	// finding the best order of line mapping
	switch (order)
	{
	case DESC:
		find_DESC_order();
		break;
	case SUM:
		find_SUM_order();
		break;
	case MAX:
		find_MAX_order();
		break;
	case CUSTOM:
		steps_ = custom_order.size();
		order_ = custom_order;
		break;
	case AUTO:
		assert(!"Order AUTO has nothing to do in the pattern constructor.");
		throw my_exception("Order AUTO has nothing to do in the pattern constructor.");
		break;
	default:
		assert(!"Unsupported order was given in the pattern constructor.");
		throw my_exception("Unsupported order was given in the pattern constructor.");
		break;
	}

	// finding mapped lines I need to remember after each line mapping
	find_what_to_remember();
	find_parralel_bound_indices();
	find_extending_order();
}

template<typename T>
bool General_pattern<T>::avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& sizes, const size_t r, const size_t c, const size_t& force_end)
{
	// There was a change from one-entry to zero-entry
	if (r != (size_t)-1 && c != (size_t)-1 && big_matrix.at(r, c) == 0)
		return true;

	// I start with empty mapping - no lines are mapped
	building_tree_[0].init();

	Counter counter;
	size_t from, to;
	const size_t big_matrix_rows = big_matrix.getRow(),
		big_matrix_cols = big_matrix.getCol();

	// main loop through added lines (loops 2*k times)
	for (size_t level = 0; level < steps_; ++level)
	{
		// the function is forced to end from outside
		if (force_end == 1)
			return false;

		counter.maps = 0;
		counter.tries = 0;
		// I cannot even map the first (ordered) i lines of the pattern, I definitely cannot map all lines of the pattern
		if (building_tree_[level % 2].size() == 0)
			return true;

		building_tree_[(level + 1) % 2].clear();

		// loop through the mappings found in the last iteration
		for (const std::vector<size_t>& mapping : building_tree_[level % 2])
		{
			// the function is forced to end from outside
			if (force_end == 1)
				return false;

			// find boundaries of added line:
			find_parallel_bounds(order_[level], level, mapping, big_matrix_rows, big_matrix_cols, from, to, r, c);

			// map orders_[level] to j-th line of N if possible
			for (size_t big_line = from; big_line < to; ++big_line)
			{
				++counter.tries;

				// if currenly added line can be mapped to big_line in mapping map
				if (map(map_approach_ > LAZY, order_[level], level, big_line, mapping, big_matrix))
				{
					++counter.maps;

					// I have mapped the last line so I have found the mapping of the pattern into the big matrix - it doesn't avoid the pattern
					if (level == steps_ - 1) {
						counter.uniques = building_tree_[(level + 1) % 2].size();
						sizes.push_back(counter);
						return false;
					}

					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line, and add it to the building_tree_
					building_tree_[(level + 1) % 2].insert_without_duplicates(extend(level, big_line, mapping));
				}
			}
		}

		counter.uniques = building_tree_[(level + 1) % 2].size();
		sizes.push_back(counter);
	}

	// after the last important line is mapped, I find out that there is no mapping of the whole pattern - matrix avoids the pattern
	return true;
}

template<typename T>
void General_pattern<T>::parallel_map(const size_t index, const Matrix<size_t>& big_matrix)
{
	// the loop end (and the thread dies) at the end of MCMCgenerator
	while (!end_)
	{
		{
			// the mutex to not allow waiting when the condition variable is being notified. This situation won't happen but this syntax is "best practice"
			std::unique_lock<std::mutex> lck(mtxs_[index]);
			// indicator that the thread sleeps
			sleeps_[index] = true;
			
			// before waiting I notify the main thread - either to let it know something was mapped or that the last mapping is being computed
			{
				// to prevent the situation main threads fall asleep right after I notify it
				std::unique_lock<std::mutex> lck2(mtx_);
				something_is_mapped_or_done_ = true;
				cv_.notify_one();
			}

			// the last mapping of current level is being computed or was computed and MCMCgenerator doesn't end 
			while (done_ && !end_)
				cvs_[index].wait(lck);

			// the worker doesn't sleep, it computes another mapping
			sleeps_[index] = false;
		}

		// a line of the big_matrix to which I want to try to map currently mapped line of the pattern
		size_t big_line = ++big_line_;

		// the line I have is beyong the last one which is worth trying or MCMCgenerator ended
		if (big_line >= big_line_to_ || end_)
		{
			done_ = true;
			continue;
		}

		// compute whether the currenly mapped line can be mapped to big_line and extend the mapping (mapping_ptr)
		if (map(map_approach_ > LAZY, order_[level_], level_, big_line, *mapping_ptr_, big_matrix))
		{
			// I have mapped last line so I have found the mapping of the pattern into the big matrix - it doesn't avoid the pattern
			if (level_ == steps_ - 1)
				done_ = true;

			// compute the extension of the mapping (mapping_ptr_) with the big_line
			auto extended = extend(level_, big_line, *mapping_ptr_);

			// and add it to the queue
			mutexes_[index].lock();
			qs_[index].push(std::move(extended));
			mutexes_[index].unlock();
		}
	}
}

template<typename T>
bool General_pattern<T>::parallel_avoid(const size_t threads_count, const Matrix<size_t>& big_matrix, std::vector<Counter>& sizes, const size_t r, const size_t c, const size_t& force_end)
{
	// There was a change from one-entry to zero-entry
	if (r != (size_t)-1 && c != (size_t)-1 && big_matrix.at(r, c) == 0)
		return true;

	// I start with empty mapping - no lines are mapped
	building_tree_[0].init();

	Counter counter; 
	size_t level = 0, from, to = 0;
	const size_t big_matrix_rows = big_matrix.getRow(),
		big_matrix_cols = big_matrix.getCol();

	// main loop through added lines (loops 2*k times)
	for (level_ = 0; level_ < steps_; ++level_)
	{
		// the function is forced to end from outside
		if (force_end == 1)
			return false;

		counter.maps = 0;
		counter.tries = 0;
		// I cannot even map the first (ordered) i lines of the pattern, I definitely cannot map all lines of the pattern
		if (building_tree_[level_ % 2].size() == 0)
			return true;

		building_tree_[(level_ + 1) % 2].clear();

		// loop through the mappings found in the last iteration
		for (const std::vector<size_t>& mapping : building_tree_[level_ % 2])
		{
			// this is how I give access to a mapping to threads
			mapping_ptr_ = &mapping;

			// the function is forced to end from outside
			if (force_end == 1)
				return false;

			// find boundaries of added line:
			find_parallel_bounds(order_[level_], level_, mapping, big_matrix_rows, big_matrix_cols, from, to, r, c);

			// this value will be wrong if the last is mapped not using the last possible line
			counter.tries += to - from;

			// the first line I want to try to map the currently mapped line to - 1
			big_line_ = from - 1;
			// the last line I want to try to map the currently mapped line to  + 1
			big_line_to_ = to;
			// all the attempts haven't been computed
			done_ = false;
			something_is_mapped_or_done_ = false;

			// wake up the threads
			for (size_t index = 0; index < threads_count - 1; ++index)
			{
				// again using a mutex here is useless since I know the thread sleeps
				std::unique_lock<std::mutex> lck(mtxs_[index]);
				cvs_[index].notify_one();
			}
			
			// indicator whether I managed to map the last mapped line
			bool extended = false;

			// waits until there is a mapping found and there is one, it adds the mapping to the building_tree_
			// if workers are done (done_) and they all sleep - meaning no more mappings are being computed - the loop breaks
			while (true)
			{
				{
					// this mutex is here to prevent the case when the process starts waiting without noticing it was just notified
					std::unique_lock<std::mutex> lck(mtx_);

					// nothing new mapped and something is still running
					while (!something_is_mapped_or_done_)
						cv_.wait(lck);

					something_is_mapped_or_done_ = false;
				}

				// the last mapping is computed or being computed
				if (done_)
				{
					// all maps are computed
					bool level_done = true;

					for (size_t i = 0; i < threads_count - 1; ++i)
						if (!sleeps_[i]) {
							level_done = false;
							break;
						}

					// the working threads all wait, meaning all map attempts of current level are computed - this is the only way to get out of the loop
					if (level_done)
						break;
				}

				std::vector<size_t> current;

				// iterate through queues and if one is not empty, take its first element and add it to building_tree_
				for (size_t index = 0; index < threads_count - 1; ++index)
				{
					{
						// need a mutex since accessing a queue which is being filled by a different thread
						std::unique_lock<std::mutex> lck(mutexes_[index]);

						if (qs_[index].empty())
							continue;

						// take the first element, pop it and unlock the mutex
						current = qs_[index].front();
						qs_[index].pop();
					}

					// I found a mapping of the pattern since the queue was not empty and I mapped the last line
					if (level_ == steps_ - 1)
					{
						done_ = true;
						extended = true;
						break;
					}

					++counter.maps;
					building_tree_[(level_ + 1) % 2].insert_without_duplicates(std::move(current));
				}
			}

			// I found atleast one mapping of the last line -> the matrix contains the pattern
			if (level_ == steps_ - 1 && extended) {
				counter.uniques = 1;
				sizes.push_back(counter);

				return false;
			}
			
			// now I know all the workers sleep but it is possible the queues are not empty - but no need to lock them up because the main thread is the only one to access them now
			for (size_t index = 0; index < threads_count - 1; ++index)
			{
				// no need to lock mutex since all workers sleep
				while (!qs_[index].empty())
				{
					// I found a mapping of the last line
					if (level_ == steps_ - 1) {
						counter.uniques = 1;
						sizes.push_back(counter);

						return false;
					}

					auto current = qs_[index].front();
					qs_[index].pop();

					++counter.maps;
					building_tree_[(level_ + 1) % 2].insert_without_duplicates(std::move(current));
				}
			}
		}

		counter.uniques = building_tree_[(level_ + 1) % 2].size();
		sizes.push_back(counter);
	}

	// after the last important line is mapped, I find out that there is no mapping of the whole pattern - matrix avoids the pattern
	return true;
}

template<typename T>
void General_pattern<T>::find_DESC_order()
{
	// vector of pairs (number of one-entries, line_ID) of nonempty lines
	std::vector<std::pair<size_t, size_t> > pairs(row_ + col_ - bit_count(empty_lines_));

	// index to the pairs vector to next element
	size_t index = 0;

	// loop through lines of the pattern, for nonempty lines it computes number of one-entries and stores it to pairs vector
	for (size_t i = 0; i < row_ + col_; ++i)
	{
		// if the line has no one-entries I won't be mapping it
		if ((empty_lines_ >> i) & 1)
			continue;

		pairs[index] = std::make_pair(row_ + col_ - bit_count(lines_[index]), index);
		++index;
	}

	// sorts the vector according to the first elements of each pair
	std::sort(pairs.begin(), pairs.end());

	// I just take the ordered vector as a order of mapping
	for (size_t i = 0; i < index; ++i)
	{
		order_[i] = pairs[i].second;
	}
}

template<typename T>
void General_pattern<T>::find_SUM_order()
{
	// queue for subsets of lines
	std::queue<size_t> q;
	q.push(0);

	// vector of distances for each subset of all lines indices - distance is a number of mapped lines I need to remember throughout the whole algorithm
	std::vector<size_t> distances(1 << (row_ + col_), (size_t)-1);
	distances[0] = 0;

	// vector for retrieving the order of line mapped after I find the shortest path from 0 lines mapped to all lines mapped
	std::vector<size_t> back_trace(1 << (row_ + col_), (size_t)-1);

	// the set of all lines I want to map
	size_t position = ((1 << (row_ + col_)) - 1) ^ empty_lines_;

	while (!q.empty())
	{
		size_t current = q.front();
		q.pop();

		// extend current subset by one if possible
		for (size_t i = 0; i < row_ + col_; ++i)
		{
			// this line is already an element of the subset or the line is empty
			if ((current >> i) & 1 || (empty_lines_ >> i) & 1)
				continue;

			// extend current subset by i-th line
			size_t supset = current | (1 << i);
			// distance from 0 is equal to distance to current + number of elements I need to remember in this step
			size_t count = (supset == position) ? distances[current] : distances[current] + count_what_to_remember(supset);

			if (count < distances[supset])
			{
				// I haven't seen this subset yet, if I did there is not reason to add it to queue, since algorithm goes through layers of the same number of lines
				if (back_trace[supset] == (size_t)-1 && supset != position)
					q.push(supset);

				// update improved distance
				distances[supset] = count;
				back_trace[supset] = i;
			}
		}
	}

	// get the best order backtracing found shortest path from 0 lines to all lines
	for (size_t i = 0; i < steps_; ++i)
	{
		order_[steps_ - 1 - i] = back_trace[position];
		position = position - (1 << back_trace[position]);
	}
}

template<typename T>
void General_pattern<T>::find_MAX_order()
{
	// queue for subsets of lines
	std::queue<size_t> q;
	q.push(0);

	// vector of distances for each subset of all lines indices - distance is a number of mapped lines I need to remember throughout the whole algorithm
	std::vector<std::pair<size_t, size_t> > distances(1 << (row_ + col_), std::make_pair((size_t)-1, 0));
	distances[0] = std::make_pair(0, 1);

	// vector for retrieving the order of line mapped after I find the shortest path from 0 lines mapped to all lines mapped
	std::vector<size_t> back_trace(1 << (row_ + col_), (size_t)-1);

	// the set of all lines I want to map
	size_t position = ((1 << (row_ + col_)) - 1) ^ empty_lines_;

	while (!q.empty())
	{
		size_t current = q.front();
		q.pop();

		// extend current subset by one if possible
		for (size_t i = 0; i < row_ + col_; ++i)
		{
			// this line is already an element of the subset or the line is empty
			if ((current >> i) & 1 || (empty_lines_ >> i) & 1)
				continue;

			// extend current subset by i-th line
			size_t supset = current | (1 << i);
			size_t number = count_what_to_remember(supset);
			// distance from 0 is equal to distance to current + number of elements I need to remember in this step
			std::pair<size_t, size_t> count = (supset == position || number < distances[current].first) ? distances[current] :
				(number == distances[current].first) ? std::make_pair(distances[current].first, distances[current].second + 1) :
				std::make_pair(number, (size_t)1);

			if (count < distances[supset])
			{
				// I haven't seen this subset yet, if I did there is not reason to add it to queue, since algorithm goes through layers of the same number of lines
				if (back_trace[supset] == (size_t)-1 && supset != position)
					q.push(supset);

				// update improved distance
				distances[supset] = count;
				back_trace[supset] = i;
			}
		}
	}

	// get the best order backtracing found shortest path from 0 lines to all lines
	for (size_t i = 0; i < steps_; ++i)
	{
		order_[steps_ - 1 - i] = back_trace[position];
		position = position - (1 << back_trace[position]);
	}
}

template<typename T>
size_t General_pattern<T>::count_what_to_remember(const size_t current) const
{
	bool needed;
	size_t count = bit_count(current);

	// go through all lines
	for (size_t j = 0; j < row_ + col_; ++j)
	{
		// if the line is not an element of current, then do not care about it
		if ((current >> j) & 1)
		{
			// I'm not adding the only one row or column - if so, I don't need to remember the line for parellel purposes 
			if (!((j == 0 && row_ == 1) || (j == row_ && col_ == 1)))
			{
				// I'm adding the first row or column and I don't remember the second one
				if (((j == 0 || j == row_) && !((current >> (j + 1)) & 1))
					||
					// I'm adding the last row or column and I don't remember the previous one
					((j == row_ - 1 || j == row_ + col_ - 1) && !((current >> (j - 1)) & 1))
					||
					// I don't remember either the previous line or the next one
					(((j > 0 && j < row_) || (j > row_ && j < row_ + col_)) &&
					(!((current >> (j - 1)) & 1) || !((current >> (j + 1)) & 1))))
					// I cannot forget this line
					continue;
			}

			needed = false;

			// go through the line's elements
			for (size_t l = 0; l < ((j < row_) ? col_ : row_); ++l)
			{
				// if I find a one-entry and I don't remember the line, which itersects with computed line at the one-entry
				if (((lines_[j] >> l) & 1) &&
					((j < row_ && !((current >> (l + row_)) & 1)) ||
					(j >= row_ && !((current >> l) & 1))))
				{
					// I cannot forget the line
					needed = true;
					break;
				}
			}

			// the line is needed since it intersects with a line, which I don't know yet, at a one-entry
			if (needed)
				continue;

			// else I don't need to remember the line and can forget it
			count--;
		}
	}

	return count;
}

template<typename T>
void General_pattern<T>::find_what_to_remember()
{
	bool needed;
	// at the beginning I don't know any line
	size_t what_do_I_know = 0;
	what_to_remember_[0] = 0;

	// go through each step of line adding and find out, what to remember in that step
	for (size_t i = 1; i <= steps_; ++i)
	{
		what_do_I_know |= 1 << order_[i - 1];
		what_to_remember_[i] = what_to_remember_[i - 1] | (1 << order_[i - 1]);

		// index into created mapping vector
		size_t index = (size_t)-1;
		map_index_[i].resize(row_ + col_);

		// go through each line and if it is something I believe I need to remember, try if I can forget it
		for (size_t j = 0; j < row_ + col_; ++j)
		{
			// if the line is not an element of current, then do not care about it
			if ((what_to_remember_[i] >> j) & 1)
			{
				// I'm not adding the only one row or column - if so, I don't need to remember the line for parellel purposes 
				if (!((j == 0 && row_ == 1) || (j == row_ && col_ == 1)))
				{
					// I'm adding the first row or column and I don't remember the second one
					if (((j == 0 || j == row_) && !((what_do_I_know >> (j + 1)) & 1))
						||
						// I'm adding the last row or column and I don't remember the previous one
						((j == row_ - 1 || j == row_ + col_ - 1) && !((what_do_I_know >> (j - 1)) & 1))
						||
						// I don't remember either the previous line or the next one
						(((j > 0 && j < row_) || (j > row_ && j < row_ + col_)) &&
						(!((what_do_I_know >> (j - 1)) & 1) || !((what_do_I_know >> (j + 1)) & 1))))
					{
						// index of j-th line in i-th level is index
						map_index_[i][j] = ++index;
						// I cannot forget this line
						continue;
					}
				}

				needed = false;

				// go through the line's elements
				for (size_t l = 0; l < ((j < row_) ? col_ : row_); ++l)
				{
					// if I find a one-entry and I don't remember the line, which itersects with computed line in the one-entry
					if (((lines_[j] >> l) & 1) &&
						(((j < row_) && !((what_do_I_know >> (l + row_)) & 1)) ||
						((j >= row_) && !((what_do_I_know >> l) & 1))))
					{
						// I cannot forget the line
						needed = true;
						break;
					}
				}

				// the line is needed since it intersects with a line, which I don't know yet, in a one-entry
				if (needed)
				{
					// index of j-th line in i-th level is index
					map_index_[i][j] = ++index;
					continue;
				}

				// else I don't need to remember the line and can forget it
				what_to_remember_[i] ^= 1 << j;
			}
		}
	}
}

template<typename T>
void General_pattern<T>::find_extending_order()
{
	for (size_t i = 0; i < steps_; ++i)
	{
		// extended mapping - elements are those indices of the big matrix I need to remember in the (i+1)-th step
		std::vector<size_t> extended;
		size_t index = 0;

		// go through all lines
		for (size_t l = 0; l < row_ + col_; ++l)
		{
			// if I remembered that line in the previous step
			if ((what_to_remember_[i] >> l) & 1)
			{
				// and I want to remember it now as well
				if ((what_to_remember_[i + 1] >> l) & 1)
					// then add it to the extended vector
					extended.push_back(index);
				// if I don't need to remember it now, just increase the index to the previous mapping
				++index;
			}
			// this only happens when I need to remember currenly adding line
			else if ((what_to_remember_[i + 1] >> l) & 1)
				extended.push_back((size_t)-1);
		}

		extending_order_[i] = extended;
	}
}

template<typename T>
void General_pattern<T>::find_parralel_bound_indices()
{
	for (size_t i = 0; i < steps_; ++i)
	{
		parallel_bound_indices_[i].resize(row_ + col_);
		find_bound_indices(order_[i], i);

		for (size_t j = 0; j < ((order_[i] < row_) ? col_ : row_); ++j)
		{
			if (((lines_[order_[i]] >> j) & 1) &&
				(!((order_[i] < row_ && ((what_to_remember_[i] >> (j + row_)) & 1)) ||
				(order_[i] >= row_ && ((what_to_remember_[i] >> j) & 1)))))
				find_bound_indices((order_[i] < row_) ? j + row_ : j, i);
		}
	}
}

template<typename T>
void General_pattern<T>::find_bound_indices(const size_t line, const size_t level)
{
	size_t	index = (size_t)-1,	// index into map m
			bot = (size_t)-1,	// index to the lower bound for currently added line in map vector
			top = (size_t)-1,	// index to the upper bound for currently added line in map vector
			i_top = 0,			// index of the line of the pattern which is a upper bound of currently added line 
			i_bot = (size_t)-1;	// index of the line of the pattern which is a lower bound of currently added line

	// go through all lines and find the lower and upper bound for "line", which I remember in i-th step of the algorithm
	for (i_top = 0; i_top < row_ + col_; ++i_top)
	{
		// iterating through rows and adding a row
		if (((line < row_ && i_top < row_) ||
			// iterating through columns and adding a column
			(line >= row_ && i_top >= row_)) &&
			// I remember this line
			((what_to_remember_[level] >> i_top) & 1))
		{
			++index;

			// I change the lower bound, I might do it a few times (up to k times)
			if (i_top < line)
			{
				i_bot = i_top;
				bot = index;
			}
			// I only set top once
			else if (i_top > line)
			{
				top = index;
				break;
			}
		}
		// I need to remeber this line, but I'm adding column and going through rows or vice versa
		else if ((what_to_remember_[level] >> i_top) & 1)
			++index;
		// if there is no upper bound
		else if ((line < row_ && i_top + 1 == row_) ||
			(line >= row_ && i_top + 1 == row_ + col_))
		{
			i_top = (size_t)-1;
			break;
		}
	}

	parallel_bound_indices_[level][line].first = std::make_pair(bot, top);
	parallel_bound_indices_[level][line].second = std::make_pair(i_bot, i_top);
}

template<typename T>
void General_pattern<T>::find_parallel_bounds(const size_t line, const size_t level, const std::vector<size_t>& mapping, const size_t rows, const size_t columns,
	size_t& from, size_t& to, const size_t r, const size_t c) const
{
	const size_t	bot = parallel_bound_indices_[level][line].first.first,		// index to the lower bound for currently added line in map vector
					top = parallel_bound_indices_[level][line].first.second,	// index to the upper bound for currently added line in map vector
					i_bot = parallel_bound_indices_[level][line].second.first,	// index of the line of the pattern which is a upper bound of currently added line 
					i_top = parallel_bound_indices_[level][line].second.second;	// index of the line of the pattern which is a lower bound of currently added line

	// i have parallel boundaries of (size_t)-1 if there are not any, return correct bounds in both cases:
	// if there is no lower bound
	if (bot == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_) {
			from = 0 + line;

			// I know r-th row is involved in every occurence of the pattern; therefore, I cannot map the last row before r-th one
			if (line == row_ - 1 && r != (size_t)-1 && r > from)
				from = r;
		}
		// "line" is a column
		else {
			from = rows - row_ + line;

			// I know c-th column is involved in every occurence of the pattern; therefore, I cannot map the last column before c-th one
			if (line == row_ + col_ - 1 && c != (size_t)-1 && c + rows > from)
				from = rows + c;
		}
	}
	// else return found lower bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
		from = mapping[bot] - i_bot + line;

	// if there is no upper bound
	if (top == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_) {
			to = rows - row_ + line + 1;

			// I know r-th row is involved in every occurence of the pattern; therefore, I cannot map the first row after r-th one
			if (line == 0 && r != (size_t)-1 && r + 1 < to)
				to = r + 1;
		}
		// "line" is a column
		else {
			to = rows + columns - row_ - col_ + line + 1;

			// I know c-th column is involved in every occurence of the pattern; therefore, I cannot map the first column after c-th one
			if (line == row_ && c != (size_t)-1 && c + rows + 1 < to)
				to = c + rows + 1;
		}
	}
	// else return found upper bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
		to = mapping[top] - i_top + line + 1;

	if (top != (size_t)-1) {
		// If I remember the next row, it has to be mapped to "top" and that line is mapped after "r" then I cannot map previous row before r-th one
		if (line < row_ && line != row_ - 1 && r != (size_t)-1 && ((what_to_remember_[level] >> (line + 1)) & 1) && mapping[top] > r && r + 1 > from)
			// if r > to the for cycle in map() will iterate through an empty set, which is ok
			from = r;
		// If I remember the next column, it has to be mapped to "top" and that line is mapped after "c" then I cannot map then previous column before c-th one
		else if (line < row_ && line != row_ + col_ - 1 && c != (size_t)-1 && ((what_to_remember_[level] >> (line + 1)) & 1) && mapping[top] > c + rows && c + rows > from)
			// if rows + c > to the for cycle in map() will iterate through an empty set, which is ok
			from = rows + c;
	}

	if (bot != (size_t)-1) {
		// If I remember the previous row, it has to be mapped to "bot" and that line is mapped before "r" then I cannot map then next row after r-th one
		if (line < row_ && line != 0 && r != (size_t)-1 && ((what_to_remember_[level] >> (line - 1)) & 1) && mapping[bot] < r && r < to)
			// if r < from the for cycle in map() will iterate through an empty set, which is ok
			to = r + 1;
		// If I remember the previous column, it has to be mapped to "bot" and that line is mapped before "c" then I cannot map then next column after c-th one
		else if (line < row_ && line != row_ && c != (size_t)-1 && ((what_to_remember_[level] >> (line - 1)) & 1) && mapping[bot] < c + rows && c + rows + 1 < to)
			// if r < from the for cycle in map() will iterate through an empty set, which is ok
			to = c + rows + 1;
	}
}

template<typename T>
bool General_pattern<T>::check_orthogonal_bounds(const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping,
	const size_t orthogonal_line, const size_t big_orthogonal_line, const Matrix<size_t>& big_matrix) const
{
	const size_t bot = parallel_bound_indices_[level][line].first.first,	// index to the lower bound for currently added line in map vector
		top = parallel_bound_indices_[level][line].first.second,			// index to the upper bound for currently added line in map vector
		i_bot = parallel_bound_indices_[level][line].second.first,			// index of the line of the pattern which is a upper bound of currently added line 
		i_top = parallel_bound_indices_[level][line].second.second;			// index of the line of the pattern which is a lower bound of currently added line

	size_t from, to, current;

	// i have parallel boundaries of (size_t)-1 if there are not any, return correct bounds in both cases:
	// if there is no lower bound
	if (bot == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_)
		{
			from = 0;
			current = 0;
		}
		// "line" is a column
		else
		{
			from = big_matrix.getRow();
			current = row_;
		}
	}
	// else return found lower bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
	{
		from = mapping[bot];
		current = i_bot;
	}

	while (current != line)
	{
		// I didn't manage to find enough one-entries
		if (big_line == from)
			return false;

		if ((orthogonal_line >> current) & 1)
		{
			if ((line < row_ && big_matrix.at(from, big_orthogonal_line)) ||
				(line >= row_ && big_matrix.at(big_orthogonal_line, from - big_matrix.getRow())))
				++current;
		}
		else
			++current;

		++from;
	}

	// if there is no upper bound
	if (top == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_) {
			to = big_line + 1;
			current = line + 1;

			while (current != row_ + 1)
			{
				// I didn't manage to find enough one-entries
				if (big_matrix.getRow() + 1 == to)
					return false;

				if ((orthogonal_line >> current) & 1)
				{
					if (big_matrix.at(to, big_orthogonal_line))
						++current;
				}
				else
					++current;

				++to;
			}
		}
		// "line" is a column
		else {
			to = big_line + 1;
			current = line + 1;

			while (current != row_ + col_ + 1)
			{
				// I didn't manage to find enough one-entries
				if (big_matrix.getRow() + big_matrix.getCol() + 1 == to)
					return false;

				if ((orthogonal_line >> current) & 1)
				{
					if (big_matrix.at(big_orthogonal_line, to - big_matrix.getRow()))
						++current;
				}
				else
					++current;

				++to;
			}
		}
	}
	// else return found upper bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
	{
		to = big_line + 1;
		current = line + 1;

		while (current != i_top + 1)
		{
			// I didn't manage to find enough one-entries
			if (to == mapping[top] + 1)
				return false;

			if ((orthogonal_line >> current) & 1)
			{
				if ((line < row_ && big_matrix.at(to, big_orthogonal_line)) ||
					(line >= row_ && big_matrix.at(big_orthogonal_line, to - big_matrix.getRow())))
					++current;
			}
			else
				++current;

			++to;
		}
	}

	return from < to;
}

template<typename T>
bool General_pattern<T>::map(const bool backtrack, const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping, const Matrix<size_t>& big_matrix)
{
	size_t from,
		to,
		last_one = 0;
	bool know_bounds = false,
		atleast1 = false;
	const bool line_is_row = line < row_,
		line_is_col = line >= row_;

	// go through all elements of "line"
	for (size_t l = 0; l < (line_is_row ? col_ : row_); ++l)
	{
		++last_one;

		// real index of l as a line of the big_matrix
		const size_t l_index = (line_is_row ? l + row_ : l);

		// if there is a one-entry at the l-th position of added line in the pattern
		if ((lines_[line] >> l) & 1)
		{
			// and I remember the l-th line
			if ((what_to_remember_[level] >> l_index) & 1)
			{
				const size_t index = map_index_[level][l_index];
				know_bounds = false;

				// and there is no one-entry at their intersection
				if ((line_is_row && (!big_matrix.at(big_line, mapping[index] - big_matrix.getRow()))) ||
					(line_is_col && (!big_matrix.at(mapping[index], big_line - big_matrix.getRow()))))
					// I can't map the line here
					return false;

				// or there is not enough one-entries on the intersected line
				if ((map_approach_ == LAZY || map_approach_ > SEMILAZY) && !check_orthogonal_bounds(line, level, big_line, mapping, lines_[l],
					(line_is_row ? mapping[index] - big_matrix.getRow() : mapping[index]), big_matrix))
					return false;
			}
			// I don't want to call myself recursively anymore - I only backtrack the function calls from avoid function
			else if (!backtrack)
				continue;
			// I have the bounds from previously found one-entry, need to find another one from the "last one" to "to"
			else if (know_bounds)
			{
				find_parallel_bounds(l_index, level, mapping, big_matrix.getRow(), big_matrix.getCol(), from, to);
				last_one = (last_one < from) ? from : last_one;
				atleast1 = false;

				// go through all possible lines of the big matrix and find a one-entry if there is any
				for (; last_one < to; ++last_one)
				{
					// is there a one-entry where I need it?
					if ((line_is_row && big_matrix.at(big_line, last_one - big_matrix.getRow())) ||
						(line_is_col && big_matrix.at(last_one, big_line - big_matrix.getRow())))
					{
						// can l-th line be mapped to last_one?
						if (map_approach_ == SUPERACTIVE || map(false, l_index, level, last_one, mapping, big_matrix))
						{
							// yes, it can
							atleast1 = true;
							break;
						}
					}
				}

				// there is not enough one-entries in the interval (last_one, to) - I cannot map "line" to j-th line
				if (!atleast1)
					return false;
			}
			// I don't have the bounds from previously found one-entry, need to find the bounds and then a one-entry in that interval
			else
			{
				// find the bounds
				find_parallel_bounds(l_index, level, mapping, big_matrix.getRow(), big_matrix.getCol(), from, to);
				know_bounds = true;
				atleast1 = false;

				// go through the lines of a big matrix, where I can map l to and check if there is a one-entry 
				for (last_one = from; last_one < to; ++last_one)
				{
					// if there is a one-entry
					if ((line_is_row && big_matrix.at(big_line, last_one - big_matrix.getRow())) ||
						(line_is_col && big_matrix.at(last_one, big_line - big_matrix.getRow())))
					{
						// can l-th line be mapped to last_one?
						if (map_approach_ == SUPERACTIVE || map(false, l_index, level, last_one, mapping, big_matrix))
						{
							// yes, it can
							atleast1 = true;
							break;
						}
					}
				}

				// there is not enough one-entries in the interval [from, to) - I cannot map "line" to the j-th line
				if (!atleast1)
					return false;
			}
		}
	}

	// "line" can be mapped to j-th line of the big matrix
	return true;
}

template<typename T>
std::vector<size_t> General_pattern<T>::extend(const size_t level, const size_t value, const std::vector<size_t>& mapping) const
{
	const size_t max = extending_order_[level].size();

	// extended mapping - elements are those indices of the big matrix I need to remember in the (i+1)-th step
	std::vector<size_t> extended(max);

	// go through elements of extended mapping and write those elements from previous mapping that should be there
	for (size_t l = 0; l < max; ++l)
	{
		if (extending_order_[level][l] == (size_t)-1)
			extended[l] = value;
		else
			extended[l] = mapping[extending_order_[level][l]];
	}

	return extended;
}

#endif