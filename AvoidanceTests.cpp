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
	  orders_(k << 1),
	  what_to_remember_((k << 1) + 1),
	  building_tree_((k << 1) + 1)
{
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			if (pattern.at(i, j))
			{
				lines_[i] |= 1 << j;
				lines_[j + k] |= 1 << i;
			}
		}
	}
	find_DESC_orders();
	// TODO find orders (DAG method)
	find_what_to_remember();
}

bool general_pattern::avoid(const matrix<size_t>& N)
{
	// I'm not using the fact I know the position which has changed
	for (size_t i = 0; i <= k << 1; i++)
	{
		building_tree_[i].clear();
	}	
	size_t from, to;
	building_tree_[0].push_back(std::vector<size_t>(0));

	for (size_t i = 0; i < steps; i++) // main loop through added lines (loops 2*k times)
	{
		if (building_tree_[i].size() == 0)
			return true;
		for (size_t m = 0; m < building_tree_[i].size(); m++) // loop through the mappings found in the last iteration
		{
			// find boundaries of added line:
			find_parallel_bounds(orders_[i], i, m, N.getRow(), N.getCol(), from, to);

			for (size_t j = from; j < to; j++) // map orders_[i] to j-th line of N if possible
			{	
				if (map(i, j, m, N)) // if currenly added line can be mapped to j in mapping m
				{
					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line
					extend(i, j, m);
				}
			}
		}
	}
	if (building_tree_[steps].empty())
		return true;
	return false;
}

void general_pattern::find_DESC_orders()
{
	std::vector<std::pair<size_t, size_t> > pairs(k << 1);
	for (size_t i = 0; i < k << 1; i++)
	{
		pairs[i] = std::make_pair((k << 1) - bit_count(lines_[i]), i);
	}
	std::sort(pairs.begin(), pairs.end());
	// maybe bucket sort instead (n >> k*log k > k)
	for (size_t i = 0; i < k << 1; i++)
	{
		orders_[i] = pairs[i].second;
		if (lines_[orders_[i]] == 0)	// since it is sorted, as soon as I want to add empty line, I would be only adding empty lines and there is no reason to do that
		{
			steps = i;
			break;
		}
	}
}

void general_pattern::find_DAG_orders()
{
	// TODO
}

void general_pattern::find_what_to_remember()
{
	// what to remember for given orders
	bool needed;
	size_t what_do_I_know = 0;
	what_to_remember_[0] = 0;
	for (size_t i = 1; i <= k << 1; i++)
	{
		what_do_I_know |= 1 << orders_[i - 1];
		what_to_remember_[i] = what_to_remember_[i - 1] | (1 << orders_[i - 1]);
		for (size_t j = 0; j < k << 1; j++)
		{
			if ((what_to_remember_[i] >> j) & 1)
			{
				if (((j == 0 || j == k) && !((what_do_I_know >> (j + 1)) & 1))					// I'm adding the first row or column and I don't remember the second one 
					||
					((j == k - 1 || j == (k << 1) - 1) && !((what_do_I_know >> (j - 1)) & 1))	// I'm adding the last row or column and I don't remember the previous one
					||
					(((j > 0 && j < k) || (j > k && j < k << 1)) &&
					(!((what_do_I_know >> (j - 1)) & 1) || !((what_do_I_know >> (j + 1)) & 1))))// I don't remember either the previous line or the next one
					continue;
				needed = false;
				for (size_t l = 0; l < k; l++) // TODO k << 1
				{
					if (((lines_[j] >> l) & 1) && (((j < k) && !((what_do_I_know >> (l + k)) & 1)) || ((j < k) && !((what_do_I_know >> (l + k)) & 1))))
					{
						needed = true;
						break;
					}
				}
				if (needed)
					continue;
				what_to_remember_[i] ^= 1 << j;
			}
		}
	}
}

void general_pattern::find_parallel_bounds(const size_t line, const size_t i, const size_t m, const size_t rows, const size_t columns,
		size_t& from, size_t& to) // TODO kill all bugs
{
	size_t	index = 0 - 1,	// index into map m
		bot = 0 - 1,	// index to the lower bound for currently added line in map vector
		top = 0 - 1,	// index to the upper bound for currently added line in map vector
		i_top = 0,		// index of the line of the pattern which is a upper bound of currently added line 
		i_bot = 0 - 1;	// index of the line of the pattern which is a lower bound of currently added line
	for (i_top = 0; i_top < k << 1; i_top++)
	{
		if (((((line << 1) < orders_.size()) && ((i_top << 1) < orders_.size())) ||		// iterating through rows and adding a row
			(((line << 1) >= orders_.size()) && ((i_top << 1) >= orders_.size()))) &&	// iterating through columns and adding a column
			((what_to_remember_[i] >> i_top) & 1))										// I remember this line
		{
			index++;
			if (i_top < line) // I might change top a few times (up to k times)
			{
				i_bot = i_top;
				bot = index;
			}
			else if (i_top > line)		// I only set top once
			{
				top = index;
				break;
			}
		}
		else if ((what_to_remember_[i] >> i_top) & 1)								// I need to remeber this line, but I'm adding column and going through rows or vice versa
			index++;
		else if ((((line << 1) < orders_.size()) && (((i_top + 1) << 1) == orders_.size())) ||
			(((line << 1) >= orders_.size()) && (i_top + 1 == orders_.size())))	// if there is no upper bound
		{
			i_top = 0 - 1;
			break;
		}
	}
	// i have parallel boundaries
	if (bot == (size_t)(0 - 1))
	{
		if (line < k)
			from = 0 + line;
		else
			from = rows - k + line;
	}
	else
		from = building_tree_[i][m][bot] - i_bot + line;
	if (top == (size_t)(0 - 1))
	{
		if (line < k)
			to = rows - k + line + 1;
		else
			to = rows + columns - (k << 1) + line + 1;
	}
	else
		to = building_tree_[i][m][top] - i_top + line + 1;
	
}

bool general_pattern::map(const size_t i, const size_t j, const size_t m, const matrix<size_t>& N)
{
	for (size_t l = 0; l < k; l++)
	{
		if ((lines_[orders_[i]] >> l) & 1)						// there is a 1 entry at the l-th position of added line in the pattern
		{
			if (((orders_[i] < k) && ((what_to_remember_[i] >> (l + k)) & 1)) ||
				((orders_[i] >= k) && ((what_to_remember_[i] >> l) & 1)))		// I remember the l-th line
			{
				size_t index = 0, j2 = 0;
				for (; j2 < k << 1; j2++)
				{
					if (j2 == l)
						break;
					if ((what_to_remember_[i] >> j2) & 1)
						index++;
				}
				if (((orders_[i] < k) && (!N.at(j, building_tree_[i][m][j2] - N.getRow()))) ||
					((orders_[i] >= k) && (!N.at(building_tree_[i][m][j2], j - N.getRow()))))				// and there is no 1-entry at their intersection
				{
					return false;						// so I can't map the line here
				}
			}
			else
			{
				// TODO check if there is enough 1 entries - now I only check if there is atleast 1, but I might use that one for more than one (could be cutting more)
				size_t from, to;
				find_parallel_bounds((orders_[i] < k) ? l + k : l, i, m, N.getRow(), N.getCol(), from, to);
				bool atleast1 = false;
				for (size_t j2 = from; j2 < to; j2++)
				{
					if (((orders_[i] < k) && N.at(j, j2 - N.getRow())) ||
						((orders_[i] >= k) && N.at(j2, j - N.getRow())))
					{
						atleast1 = true;
						break;
					}
				}
				if (!atleast1)
					return false;
			}
		}
	}
	return true;
}

void general_pattern::extend(const size_t i, const size_t j, const size_t m)
{
	std::vector<size_t> extended;
	size_t index = 0;
	for (size_t l = 0; l < (k << 1); l++)
	{
		if ((what_to_remember_[i] >> l) & 1)
		{
			if ((what_to_remember_[i + 1] >> l) & 1)
				extended.push_back(building_tree_[i][m][index]);
			index++;
		}
		else if ((what_to_remember_[i + 1] >> l) & 1)	// this only happens when I need to remember currenly added line
			extended.push_back(j);
	}
	bool mapped = false;
	for (size_t m2 = 0; m2 < building_tree_[i + 1].size(); m2++)
	{
		mapped = true;
		for (size_t l2 = 0; l2 < building_tree_[i + 1][0].size(); l2++)
		{
			if (extended[l2] != building_tree_[i + 1][m2][l2])
			{
				mapped = false;
				break;
			}
		}
		if (mapped)	// extended has already been added (atleast its different class)
			break;
	}
	if (!mapped)
		building_tree_[i + 1].push_back(extended);
}

/* Walking pattern */
walking_pattern::walking_pattern(const matrix<size_t>& pattern, const size_t n)
	: max_walk_part_(n, n, std::pair<size_t, size_t>(0, 0))
{
	size_t last_i = 0, last_j = 0, j;	// coords of the last visited 1 entry, index for the columns
	// all elements on the same diagonal have the same sum of their coordinates
	for (size_t sum = 0; sum < pattern.getRow() + pattern.getCol() - 1; sum++)
	{
		for (size_t i = 0; i <= sum; i++)	// index for the rows
		{
			j = sum - i;
			if (i >= pattern.getCol())
				break;
			if (j >= pattern.getRow())
				continue;
			// when I find 1 entry or find myself on the last diagonal
			if (pattern.at(i, j) || sum == pattern.getRow() + pattern.getCol() - 2)
			{
				// last visited element is 0 and I did not find any 1 entries
				if (!pattern.at(i, j) && last_i == 0 && last_j == 0 && !pattern.at(last_i, last_j)) 
					throw std::invalid_argument("Pattern has no one entries.");
				if (j < last_j || i < last_i)
					throw std::invalid_argument("Pattern is not a walking pattern type.");
				// need to find, which elements will be a part of the walk
				// from the last found 1 entry go as for to the bottom as you can
				for (size_t i2 = last_i; i2 < i; i2++)
				{
					value_.push_back(pattern.at(i2, last_j));
					direction_.push_back(0);					// vertical
				}
				// then go to the right and stop right before [i,j]
				for (size_t j2 = last_j; j2 < j; j2++)
				{
					value_.push_back(pattern.at(i, j2));
					direction_.push_back(1);					// horizontal
				}
				last_i = i; last_j = j;							// last element of the walk
			}
		}
	}
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

		if (current.first == 0)		// element on the first row
			c_v_v = 0;
		else
			c_v_v = max_walk_part_.at(current.first - 1, current.second).first;
		if (current.second == 0)	// element on the first column
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
			if (c_v_v + 1 == value_.size()) // I found the last element of the walk
				return false;
			if (direction_[c_v_v]) // walk continues to the right
			{
				if (max_walk_part_.at(current).second < c_v_v + 1)
					max_walk_part_.at(current).second = c_v_v + 1;
			}
			else // walk continues to the bottom
			{
				if (max_walk_part_.at(current).first < c_v_v + 1)
					max_walk_part_.at(current).first = c_v_v + 1;
			}
		}
		if (N.at(current) || !value_[c_h_h]) // N[i,j] == 1
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
	return true;
}

#endif