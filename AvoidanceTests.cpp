#ifndef AvoidanceTest_cpp_
#define AvoidanceTest_cpp_

#include "AvoidanceTests.hpp"
#include <queue>
#include <algorithm>
#include <assert.h>

/* General pattern */
general_pattern::general_pattern(const matrix<size_t>& pattern, Order order, Map map_approach, std::vector<size_t>&& custom_order)
	: map_approach_(map_approach),
	  row_(pattern.getRow()),
	  col_(pattern.getCol()),
	  steps_(row_ + col_),
	  empty_lines_(0),
	  lines_(row_ + col_),
	  order_(row_ + col_),
	  what_to_remember_(row_ + col_ + 1),
	  parallel_bound_indices_(row_ + col_),
	  extending_order_(row_ + col_),
	  map_index_(row_ + col_ + 1)
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
		throw new std::exception("Order AUTO has nothing to do in the pattern constructor.");
		break;
	default:
		assert(!"Unsupported order was given in the pattern constructor.");
		throw new std::exception("Unsupported order was given in the pattern constructor.");
		break;
	}

	// finding mapped lines I need to remember after each line mapping
	find_what_to_remember();
	find_parralel_bound_indices();
	find_extending_order();
}

bool general_vector_pattern::avoid(const matrix<size_t>& N)
{
	// I start with empty mapping - no lines are mapped
	building_tree_[0].clear();
	building_tree_[0].push_back(std::vector<size_t>(0));	

	size_t from, to,
		big_matrix_rows = N.getRow(),
		big_matrix_cols = N.getCol();

	// main loop through added lines (loops 2*k times)
	for (size_t level = 0; level < steps_; ++level)
	{
		// I cannot even map the first (ordered) i lines of the pattern, I definitely cannot map all lines of the pattern
		if (building_tree_[level % 2].size() == 0)
			return true;

		building_tree_[(level + 1) % 2].clear();

		// loop through the mappings found in the last iteration
		for (auto& mapping : building_tree_[level % 2])
		{
			// find boundaries of added line:
			find_parallel_bounds(order_[level], level, mapping, big_matrix_rows, big_matrix_cols, from, to);

			// map orders_[level] to j-th line of N if possible
			for (size_t big_line = from; big_line < to; ++big_line) 
			{	
				// if currenly added line can be mapped to big_line in mapping map
				if (map((map_approach_ == NORECURSION) ? false : true, order_[level], level, big_line, mapping, N))
				{
					// I have mapped last line so I have found the mapping of the pattern into the big matrix - it doesn't avoid the pattern
					if (level == steps_ - 1)
						return false;

					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line
					std::vector<size_t> extended = extend(level, big_line, mapping);

					bool mapped = false;

					// go through all already found mappings in (i+1)-th step and check if extended is not already in there
					for (auto& mapping2 : building_tree_[(level + 1) % 2])
					{
						mapped = true;

						// go through m2 mapping
						for (size_t l2 = 0; l2 < building_tree_[(level + 1) % 2][0].size(); ++l2)
						{
							// and check whether the index in extended and m2 are the same
							if (extended[l2] != mapping2[l2])
							{
								// if not, extended is a different mapping then m2
								mapped = false;
								break;
							}
						}

						// extended has already been added (atleast its different class) - I won't add it for the second time
						if (mapped)
							break;
					}

					// if extended is not yet an element, add it to the tree
					if (!mapped)
						building_tree_[(level + 1) % 2].push_back(extended);
				}
			}
		}
	}

	// after the last important line is mapped, I find out that there is no mapping of the whole pattern - matrix avoids the pattern
	return true;
}

bool general_set_pattern::avoid(const matrix<size_t>& N)
{
	// I start with empty mapping - no lines are mapped
	building_tree_[0].clear();
	building_tree_[0].insert(std::vector<size_t>(0));

	size_t from, to,
		big_matrix_rows = N.getRow(),
		big_matrix_cols = N.getCol();

	// main loop through added lines (loops 2*k times)
	for (size_t level = 0; level < steps_; ++level)
	{
		// I cannot even map the first (ordered) i lines of the pattern, I definitely cannot map all lines of the pattern
		if (building_tree_[level % 2].size() == 0)
			return true;

		building_tree_[(level + 1) % 2].clear();

		// loop through the mappings found in the last iteration
		for (auto& mapping : building_tree_[level % 2])
		{
			// find boundaries of added line:
			find_parallel_bounds(order_[level], level, mapping, big_matrix_rows, big_matrix_cols, from, to);

			// map orders_[level] to j-th line of N if possible
			for (size_t big_line = from; big_line < to; ++big_line)
			{
				// if currenly added line can be mapped to big_line in mapping map
				if (map((map_approach_ == NORECURSION) ? false : true, order_[level], level, big_line, mapping, N))
				{
					// I have mapped last line so I have found the mapping of the pattern into the big matrix - it doesn't avoid the pattern
					if (level == steps_ - 1)
						return false;

					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line
					std::vector<size_t> extended = extend(level, big_line, mapping);

					building_tree_[(level + 1) % 2].insert(extended);
				}
			}
		}
	}

	// after the last important line is mapped, I find out that there is no mapping of the whole pattern - matrix avoids the pattern
	return true;
}

void general_pattern::find_DESC_order()
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
		index++;
	}

	// sorts the vector according to the first elements of each pair
	std::sort(pairs.begin(), pairs.end());	

	// I just take the ordered vector as a order of mapping
	for (size_t i = 0; i < index; ++i)
	{
		order_[i] = pairs[i].second;
	}
}

void general_pattern::find_SUM_order()
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
		order_[row_ + col_ - 1 - i] = back_trace[position];
		position = position - (1 << back_trace[position]);
	}
}

void general_pattern::find_MAX_order()
{
	// queue for subsets of lines
	std::queue<size_t> q;
	q.push(0);

	// vector of distances for each subset of all lines indices - distance is a number of mapped lines I need to remember throughout the whole algorithm
	std::vector<std::pair<size_t, size_t> > distances(1 << (row_ + col_), std::make_pair((size_t)-1, 0));
	distances[0] = std::make_pair(0,1);

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
			std::pair<size_t, size_t> count =	(supset == position || number < distances[current].first) ? distances[current] :
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
		order_[row_ + col_ - 1 - i] = back_trace[position];
		position = position - (1 << back_trace[position]);
	}
}

size_t general_pattern::count_what_to_remember(const size_t current)
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

void general_pattern::find_what_to_remember()
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

void general_pattern::find_extending_order()
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

void general_pattern::find_parralel_bound_indices()
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

void general_pattern::find_bound_indices(const size_t line, const size_t level)
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

void general_pattern::find_parallel_bounds(const size_t line, const size_t level, const std::vector<size_t>& mapping, const size_t rows, const size_t columns,
		size_t& from, size_t& to)
{
	size_t	bot = parallel_bound_indices_[level][line].first.first,		// index to the lower bound for currently added line in map vector
			top = parallel_bound_indices_[level][line].first.second,	// index to the upper bound for currently added line in map vector
			i_bot = parallel_bound_indices_[level][line].second.first,	// index of the line of the pattern which is a upper bound of currently added line 
			i_top = parallel_bound_indices_[level][line].second.second;	// index of the line of the pattern which is a lower bound of currently added line

	// i have parallel boundaries of (size_t)-1 if there are not any, return correct bounds in both cases:
	// if there is no lower bound
	if (bot == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_)
			from = 0 + line;
		// "line" is a column
		else
			from = rows - row_ + line;
	}
	// else return found lower bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
		from = mapping[bot] - i_bot + line;

	// if there is no upper bound
	if (top == (size_t)-1)
	{
		// and "line" is a row
		if (line < row_)
			to = rows - row_ + line + 1;
		// "line" is a column
		else
			to = rows + columns - row_ - col_ + line + 1;
	}
	// else return found upper bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
		to = mapping[top] - i_top + line + 1;
}

bool general_pattern::map(const bool backtrack, const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping, const matrix<size_t>& big_matrix)
{
	size_t from, to, last_one = 0;
	bool	know_bounds = false,
			atleast1 = false,
			line_is_row = line < row_,
			line_is_col = line >= row_;

	// go through all elements of "line"
	for (size_t l = 0; l < (line_is_row ? col_ : row_); ++l, ++last_one)
	{
		// real index of l as a line of the big_matrix
		size_t l_index = (line_is_row ? l + row_ : l);

		// if there is a one-entry at the l-th position of added line in the pattern
		if ((lines_[line] >> l) & 1)						
		{
			// and I remember the l-th line
			if ((what_to_remember_[level] >> l_index) & 1)
			{
				size_t index = map_index_[level][l_index];
				know_bounds = false;
				
				// and there is no one-entry at their intersection
				if ((line_is_row && (!big_matrix.at(big_line, mapping[index] - big_matrix.getRow()))) ||
					(line_is_col && (!big_matrix.at(mapping[index], big_line - big_matrix.getRow()))))
					// I can't map the line here
					return false;
			}
			// I don't want to backtrack anymore - I only backtract for line being mapped from avoid algorithm
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
						if (map_approach_ == COMPROMISE || map(false, l_index, level, last_one, mapping, big_matrix))
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
						if (map_approach_ == COMPROMISE || map(false, l_index, level, last_one, mapping, big_matrix))
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

std::vector<size_t> general_pattern::extend(const size_t level, const size_t value, const std::vector<size_t>& mapping)
{
	size_t max = extending_order_[level].size();

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

/* Walking pattern */
walking_pattern::walking_pattern(const matrix<size_t>& pattern, const size_t n)
	: max_walk_part_(n, n, std::pair<size_t, size_t>(0, 0))
{
	// coords of the last visited 1 entry, index for the columns
	size_t last_i = 0, last_j = 0, j;

	// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
	for (size_t sum = 0; sum < pattern.getRow() + pattern.getCol() - 1; ++sum)
	{
		// go through indices for the rows
		for (size_t i = 0; i <= sum; ++i)	
		{
			j = sum - i;

			// I look to the right of the pattern
			if (i >= pattern.getRow())
				break;

			// I look under the pattern
			if (j >= pattern.getCol())
				continue;

			// when I find one-entry or find myself on the last diagonal
			if (pattern.at(i, j) || sum == pattern.getRow() + pattern.getCol() - 2)
			{
				// last visited element is 0 and I did not find any 1 entries
				if (!pattern.at(i, j) && last_i == 0 && last_j == 0 && !pattern.at(last_i, last_j)) {
					assert(!"Pattern has no one entries.");
					throw std::invalid_argument("Pattern has no one entries.");
				}

				// I have found a one-entry somewhere where it shouldn't have been
				if (j < last_j || i < last_i) {
					assert(!"Pattern is not a walking pattern type");
					throw std::invalid_argument("Pattern is not a walking pattern type.");
				}

				// need to find, which elements will be a part of the walk
				// from the previously found one-entry go as for to the bottom as you can
				for (size_t i2 = last_i; i2 < i; ++i2)
				{
					value_.push_back(pattern.at(i2, last_j));
					// vertical
					direction_.push_back(0);					
				}

				// then go to the right and stop right before [i,j]
				for (size_t j2 = last_j; j2 < j; ++j2)
				{
					value_.push_back(pattern.at(i, j2));
					// horizontal
					direction_.push_back(1);					
				}
				// last element of the walk
				last_i = i; last_j = j;							
			}
		}
	}

	// add the last element of the walk
	value_.push_back(pattern.at(last_i, last_j));
}

bool walking_pattern::avoid(const size_t r, const size_t c, const matrix<size_t>& big_matrix)
{
	typedef std::pair<size_t, size_t> pair;
	std::queue<pair> q;						// queue for elements of the matrix that are supposed to be updated
	pair current;							// [x,y] of the currently updated element
	size_t old_c_v, old_c_h, c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of current, c_h of element to the left

	q.push(pair(r, c));
	while (!q.empty())
	{
		current = q.front();
		q.pop();
			
		old_c_v = max_walk_part_.at(current).first;
		old_c_h = max_walk_part_.at(current).second;

		// element on the first row
		if (current.first == 0)		
			c_v_v = 0;
		else
			c_v_v = max_walk_part_.at(current.first - 1, current.second).first;

		// element on the first column
		if (current.second == 0)	
			c_h_h = 0;
		else
			c_h_h = max_walk_part_.at(current.first, current.second - 1).second;
			
	// Initialization - copying those already found walks
		max_walk_part_.at(current).first = c_v_v;
		max_walk_part_.at(current).second = c_h_h;
			
	// Search for longer part of the walk
		// b == 1 or v_{c_v_v + 1} == 0
		if (big_matrix.at(current) || !value_[c_v_v])
		{
			// I found the last element of the walk
			if (c_v_v + 1 == value_.size()) 
				return false;

			// walk continues to the right
			if (direction_[c_v_v]) 
			{
				if (max_walk_part_.at(current).second < c_v_v + 1)
					max_walk_part_.at(current).second = c_v_v + 1;
			}
			// walk continues to the bottom
			else 
			{
				if (max_walk_part_.at(current).first < c_v_v + 1)
					max_walk_part_.at(current).first = c_v_v + 1;
			}
		}

		// N[i,j] == 1
		if (big_matrix.at(current) || !value_[c_h_h]) 
		{
			if (c_h_h + 1 == value_.size())
				return false;
			if (direction_[c_h_h])
			{
				if (max_walk_part_.at(current).second < c_h_h + 1)
					max_walk_part_.at(current).second = c_h_h + 1;
			}
			else
			{
				if (max_walk_part_.at(current).first < c_h_h + 1)
					max_walk_part_.at(current).first = c_h_h + 1;
			}
		}		

		// c_v was changed and there is still an element below the current one
		if (max_walk_part_.at(current).first != old_c_v && current.first + 1 < big_matrix.getRow())
			// if queue is not empty, check whether the element wasn't added before
			if (q.empty() || q.back() != pair(current.first + 1, current.second))	
				q.push(pair(current.first + 1, current.second));

		// c_h was changed and there is still an element to the right
		if (max_walk_part_.at(current).second != old_c_h && current.second + 1 < big_matrix.getCol())
			q.push(pair(current.first, current.second + 1));
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

#endif