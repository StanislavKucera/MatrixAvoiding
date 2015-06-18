#ifndef AvoidanceTest_cpp_
#define AvoidanceTest_cpp_

#include "AvoidanceTests.hpp"
#include <queue>
#include <algorithm>

/* General pattern */
general_pattern::general_pattern(const matrix<size_t>& pattern)
	: k(pattern.getCol()),
	  steps(k << 1),
	  lines_(k << 1),
	  order_(k << 1),
	  what_to_remember_((k << 1) + 1),
	  building_tree_((k << 1) + 1)
{
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			// reading lines of pattern and storing them memory efficiently
			if (pattern.at(i, j))	
			{
				lines_[i] |= 1 << j;
				lines_[j + k] |= 1 << i;
			}
		}
	}

	// finding the best order of line mapping
	//find_DESC_order();
	find_DAG_order();	
	// finding mapped lines I need to remember after each line mapping
	find_what_to_remember();		
}

bool general_pattern::avoid(const matrix<size_t>& N)
{
	// I'm not using the fact I know the position which has changed
	for (size_t i = 0; i <= k << 1; i++)	
	{
		// therefore I start with an empty tree
		building_tree_[i].clear();		
	}	
	// I start with empty mapping - no lines are mapped
	building_tree_[0].push_back(std::vector<size_t>(0));	

	size_t from, to, dont_care;

	// main loop through added lines (loops 2*k times)
	for (size_t i = 0; i < steps; i++)
	{
		// I cannot even map the first (ordered) i lines of the pattern, I definitely cannot map all lines of the pattern
		if (building_tree_[i].size() == 0)	
			return true;

		// loop through the mappings found in the last iteration
		for (size_t m = 0; m < building_tree_[i].size(); m++) 
		{
			// find boundaries of added line:
			find_parallel_bounds(order_[i], i, m, N.getRow(), N.getCol(), from, to, dont_care);

			// map orders_[i] to j-th line of N if possible
			for (size_t j = from; j < to; j++) 
			{	
				// if currenly added line can be mapped to j in mapping m
				if (map(true, order_[i], i, j, m, N)) 
				{
					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line
					extend(i, j, m);
				}
			}
		}
	}

	// after the last important line is mapped, I find out that there is no mapping of the whole pattern - matrix avoids the pattern
	if (building_tree_[steps].empty())
		return true;

	// I have found the mapping of the pattern into the big matrix - it doesn't avoid the pattern
	return false;
}

void general_pattern::find_DESC_order()
{
	// vector of pairs (number of one-entries, line_ID)
	std::vector<std::pair<size_t, size_t> > pairs(k << 1);

	for (size_t i = 0; i < k << 1; i++)
	{
		pairs[i] = std::make_pair((k << 1) - bit_count(lines_[i]), i);
	}

	// sorts the vector according to the first elements of each pair
	// maybe bucket sort instead (n >> k*log k > k)
	std::sort(pairs.begin(), pairs.end());	

	// I just take the ordered vector as a order of mapping
	for (size_t i = 0; i < k << 1; i++)
	{
		order_[i] = pairs[i].second;

		// since it is sorted, as soon as I am adding an empty line, I would be only adding empty lines and there is no reason to do that
		if (lines_[order_[i]] == 0)	
		{
			// so I just reduce the number of lines which need to be mapped
			steps = i;
			break;
		}
	}
}

void general_pattern::find_DAG_order()
{
	// queue for subsets of lines
	std::queue<size_t> q;
	q.push(0);

	// vector of distances for each subset of all lines indices - distance is a number of mapped lines I need to remember throughout the whole algorithm
	std::vector<size_t> distances(1 << (k << 1), (size_t)-1);
	distances[0] = 0;

	// vector for retrieving the order of line mapped after I find the shortest path from 0 lines mapped to all lines mapped
	std::vector<size_t> back_trace(1 << (k << 1), (size_t)-1);

	while (!q.empty())
	{
		size_t current = q.front();
		q.pop();

		// extend current subset by one if possible
		for (size_t i = 0; i < k << 1; i++)
		{
			// this line is already an element of the subset
			if ((current >> i) & 1)
				continue;

			// extend current subset by i-th line
			size_t supset = current | (1 << i);
			// distance from 0 is equal to distance to current + number of elements I need to remember in this step
			size_t count = distances[current] + count_what_to_remember(supset);

			if (count < distances[supset])
			{
				// I haven't seen this subset yet, if I did there is not reason to add it to queue, since algorithm goes through layers of the same number of lines
				if (back_trace[supset] == (size_t)-1)
					q.push(supset);

				// update improved distance
				distances[supset] = count;
				back_trace[supset] = i;
			}
		}
	}

	size_t position = (1 << (k << 1)) - 1;

	// get the best order backtracing found shortest path from 0 lines to all lines
	for (size_t i = 0; i < k << 1; i++)
	{
		order_[(k << 1) - 1 - i] = back_trace[position];
		position = position - (1 << back_trace[position]);
	}
}

size_t general_pattern::count_what_to_remember(size_t current)
{
	bool needed;
	size_t count = bit_count(current);

	// go through all lines
	for (size_t j = 0; j < k << 1; j++)
	{
		// if the line is not an element of current, then do not care about it
		if ((current >> j) & 1)
		{
				// I'm adding the first row or column and I don't remember the second one 
			if (((j == 0 || j == k) && !((current >> (j + 1)) & 1))					
				||
				// I'm adding the last row or column and I don't remember the previous one
				((j == k - 1 || j == (k << 1) - 1) && !((current >> (j - 1)) & 1))	
				||
				// I don't remember either the previous line or the next one
				(((j > 0 && j < k) || (j > k && j < k << 1)) &&
				(!((current >> (j - 1)) & 1) || !((current >> (j + 1)) & 1))))
				// I cannot forget this line
				continue;

			needed = false;

			// go through the line's elements
			for (size_t l = 0; l < k; l++)
			{
				// if I find a one-entry and I don't remember the line, which itersects with computed line at the one-entry
				if (((lines_[j] >> l) & 1) &&
					(((j < k) && !((current >> (l + k)) & 1)) ||
					((j >= k) && !((current >> l) & 1))))
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
	for (size_t i = 1; i <= steps; i++)
	{
		what_do_I_know |= 1 << order_[i - 1];
		what_to_remember_[i] = what_to_remember_[i - 1] | (1 << order_[i - 1]);

		// go through each line and if it is something I believe I need to remember, try if I can forget it
		for (size_t j = 0; j < k << 1; j++)
		{
			// if the line is not an element of current, then do not care about it
			if ((what_to_remember_[i] >> j) & 1)
			{	
					// I'm adding the first row or column and I don't remember the second one 
				if (((j == 0 || j == k) && !((what_do_I_know >> (j + 1)) & 1))					
					||
					// I'm adding the last row or column and I don't remember the previous one
					((j == k - 1 || j == (k << 1) - 1) && !((what_do_I_know >> (j - 1)) & 1))	
					||
					// I don't remember either the previous line or the next one
					(((j > 0 && j < k) || (j > k && j < k << 1)) &&
					(!((what_do_I_know >> (j - 1)) & 1) || !((what_do_I_know >> (j + 1)) & 1))))
					// I cannot forget this line
					continue;

				needed = false;

				// go through the line's elements
				for (size_t l = 0; l < k; l++)
				{
					// if I find a one-entry and I don't remember the line, which itersects with computed line at the one-entry
					if (((lines_[j] >> l) & 1) && 
						(((j < k) && !((what_do_I_know >> (l + k)) & 1)) || 
						 ((j >= k) && !((what_do_I_know >> l) & 1))))
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
				what_to_remember_[i] ^= 1 << j;
			}
		}
	}
}

void general_pattern::find_parallel_bounds(const size_t line, const size_t i, const size_t m, const size_t rows, const size_t columns,
		size_t& from, size_t& to, size_t& top_bound)
{
	size_t	index = 0 - 1,	// index into map m
			bot = 0 - 1,	// index to the lower bound for currently added line in map vector
			top = 0 - 1,	// index to the upper bound for currently added line in map vector
			i_top = 0,		// index of the line of the pattern which is a upper bound of currently added line 
			i_bot = 0 - 1;	// index of the line of the pattern which is a lower bound of currently added line

	// go through all lines and find the lower and upper bound for "line", which I remember in i-th step of the algorithm
	for (i_top = 0; i_top < k << 1; i_top++)
	{
			// iterating through rows and adding a row
		if (((((line << 1) < order_.size()) && ((i_top << 1) < order_.size())) ||	
			// iterating through columns and adding a column
			(((line << 1) >= order_.size()) && ((i_top << 1) >= order_.size()))) &&	
			// I remember this line
			((what_to_remember_[i] >> i_top) & 1))										
		{
			index++;

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
		else if ((what_to_remember_[i] >> i_top) & 1)								
			index++;
		// if there is no upper bound
		else if ((((line << 1) < order_.size()) && (((i_top + 1) << 1) == order_.size())) ||
				 (((line << 1) >= order_.size()) && (i_top + 1 == order_.size())))	
		{
			i_top = 0 - 1;
			break;
		}
	}

	// i have parallel boundaries of 0-1 if there are not any, return correct bounds in both cases:
	// if there is no lower bound
	if (bot == (size_t)(0 - 1))
	{
		// and "line" is a row
		if (line < k)
			from = 0 + line;
		// "line" is a column
		else
			from = rows - k + line;
	}
	// else return found lower bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
		from = building_tree_[i][m][bot] - i_bot + line;

	// if there is no upper bound
	if (top == (size_t)(0 - 1))
	{
		// and "line" is a row
		if (line < k)
		{
			to = rows - k + line + 1;
			top_bound = k;
		}
		// "line" is a column
		else
		{
			to = rows + columns - (k << 1) + line + 1;
			top_bound = k << 1;
		}
	}
	// else return found upper bound with offset, which ensures there at least enough lines in the big matrix to map those from the pattern, which are in between
	else
	{
		to = building_tree_[i][m][top] - i_top + line + 1;
		top_bound = i_top;
	}
}

bool general_pattern::map(const bool backtrack, const size_t line, const size_t i, const size_t j, const size_t m, const matrix<size_t>& N)
{
	size_t from, to, top_bound = 0, last_one;

	// go through all elements of "line"
	for (size_t l = 0; l < k; l++)
	{
		// if there is a one-entry at the l-th position of added line in the pattern
		if ((lines_[line] >> l) & 1)						
		{
			// and I remember the l-th line
			if (((line < k) && ((what_to_remember_[i] >> (l + k)) & 1)) ||
				((line >= k) && ((what_to_remember_[i] >> l) & 1)))		
			{
				size_t index = 0, j2 = 0;

				// go through line and find the l-th one - I do this to find the index of l-th line in the mapping
				for (; j2 < k << 1; j2++)
				{
					// this is the l-th line, index is now pointing to the line of the big matrix, which l-th line was mapped to
					if (((line < k) && (j2 == l + k)) ||
						((line >= k) && (j2 == l)))
						break;

					// I remember this one, so I need to increment the index
					if ((what_to_remember_[i] >> j2) & 1)
						index++;
				}

				// and there is no one-entry at their intersection
				if (((line < k) && (!N.at(j, building_tree_[i][m][index] - N.getRow()))) ||
					((line >= k) && (!N.at(building_tree_[i][m][index], j - N.getRow()))))		
					// I can't map the line here
					return false;
			}
			// I don't want to backtrack anymore - I only backtract for line being mapped from avoid algorithm
			else if (!backtrack)
				continue;
			// I have the bounds from previously found one-entry, need to find another one from the "last one" to "to"
			else if (((line < k) && ((l + k) < top_bound)) ||
					 ((line >= k) && (l < top_bound)))		
			{
				// I don't want to find the same one-entry again, so I increment the index of the line
				last_one++;		
				bool atleast1 = false;

				// go through all possible lines of the big matrix and find a one-entry if there is any
				for (; last_one < to; last_one++)
				{
					// is there a one-entry where I need it?
					if (((line < k) && N.at(j, last_one - N.getRow())) ||
						((line >= k) && N.at(last_one, j - N.getRow())))
					{
						//if (map(false, (line < k) ? (l + k) : l, i, last_one, m, N))
						//{
							atleast1 = true;
							break;
						//}
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
				find_parallel_bounds((line < k) ? l + k : l, i, m, N.getRow(), N.getCol(), from, to, top_bound);
				bool atleast1 = false;

				// go through the lines of a big matrix, where I can map l to and check if there is a one-entry 
				for (last_one = from; last_one < to; last_one++)
				{
					// if there is a one-entry
					if (((line < k) && N.at(j, last_one - N.getRow())) ||
						((line >= k) && N.at(last_one, j - N.getRow())))
					{
						//if (map(false, (line < k) ? (l + k) : l, i, last_one, m, N))
						//{
							atleast1 = true;
							break;
						//}
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

void general_pattern::extend(const size_t i, const size_t j, const size_t m)
{
	// extended mapping - elements are those indices of the big matrix I need to remember in the (i+1)-th step
	std::vector<size_t> extended;
	size_t index = 0;

	// go through all lines
	for (size_t l = 0; l < (k << 1); l++)
	{
		// if I remembered that line in the previous step
		if ((what_to_remember_[i] >> l) & 1)
		{
			// and I want to remember it now as well
			if ((what_to_remember_[i + 1] >> l) & 1)
				// then add it to the extended vector
				extended.push_back(building_tree_[i][m][index]);
			// if I don't need to remember it now, just increase the index to the previous mapping
			index++;
		}
		// this only happens when I need to remember currenly added line
		else if ((what_to_remember_[i + 1] >> l) & 1)	
			extended.push_back(j);
	}

	bool mapped = false;

	// go through all already found mappings in (i+1)-th step and check if extended is not already in there
	for (size_t m2 = 0; m2 < building_tree_[i + 1].size(); m2++)
	{
		mapped = true;

		// go through m2 mapping
		for (size_t l2 = 0; l2 < building_tree_[i + 1][0].size(); l2++)
		{
			// and check whether the index in extended and m2 are the same
			if (extended[l2] != building_tree_[i + 1][m2][l2])
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
		building_tree_[i + 1].push_back(extended);
}

/* Walking pattern */
walking_pattern::walking_pattern(const matrix<size_t>& pattern, const size_t n)
	: max_walk_part_(n, n, std::pair<size_t, size_t>(0, 0))
{
	// coords of the last visited 1 entry, index for the columns
	size_t last_i = 0, last_j = 0, j;

	// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
	for (size_t sum = 0; sum < pattern.getRow() + pattern.getCol() - 1; sum++)
	{
		// go through indices for the rows
		for (size_t i = 0; i <= sum; i++)	
		{
			j = sum - i;

			// I look to the right of the pattern
			if (i >= pattern.getCol())
				break;

			// I look under the pattern
			if (j >= pattern.getRow())
				continue;

			// when I find one-entry or find myself on the last diagonal
			if (pattern.at(i, j) || sum == pattern.getRow() + pattern.getCol() - 2)
			{
				// last visited element is 0 and I did not find any 1 entries
				if (!pattern.at(i, j) && last_i == 0 && last_j == 0 && !pattern.at(last_i, last_j)) 
					throw std::invalid_argument("Pattern has no one entries.");

				// I have found a one-entry somewhere where it shouldn't have been
				if (j < last_j || i < last_i)
					throw std::invalid_argument("Pattern is not a walking pattern type.");

				// need to find, which elements will be a part of the walk
				// from the previously found one-entry go as for to the bottom as you can
				for (size_t i2 = last_i; i2 < i; i2++)
				{
					value_.push_back(pattern.at(i2, last_j));
					// vertical
					direction_.push_back(0);					
				}

				// then go to the right and stop right before [i,j]
				for (size_t j2 = last_j; j2 < j; j2++)
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

bool walking_pattern::avoid(const size_t r, const size_t c, const matrix<size_t>& N)
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
		if (N.at(current) || !value_[c_v_v])
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
		if (N.at(current) || !value_[c_h_h]) 
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
		if (max_walk_part_.at(current).first != old_c_v && current.first + 1 < N.getRow())
			// if queue is not empty, check whether the element wasn't added before
			if (q.empty() || q.back() != pair(current.first + 1, current.second))	
				q.push(pair(current.first + 1, current.second));

		// c_h was changed and there is still an element to the right
		if (max_walk_part_.at(current).second != old_c_h && current.second + 1 < N.getCol())
			q.push(pair(current.first, current.second + 1));
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

#endif