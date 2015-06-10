#ifndef AvoidanceTest_cpp_
#define AvoidanceTest_cpp_

#include "AvoidanceTests.hpp"
#include <queue>

/* General pattern */
general_pattern::general_pattern(const matrix<size_t>& pattern, const size_t N)
	: k(pattern.getCol()),
	  rows_(k),
	  cols_(k),
	  lines_(k << 1),
	  orders_(k << 1),
	  what_to_remember_(k << 1),
	  building_tree_(k << 1)
{
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < k; j++)
		{
			if (pattern.at(i, j))
			{
				rows_[i] |= 1 << j;
				cols_[j] |= 1 << i;
				lines_[i] |= 1 << j;
				lines_[j + k] |= 1 << i;
			}
		}
	}
	for (size_t i = 0; i < orders_.size(); i++)
	{
		orders_[i] = i;
		if (i > 0)
			what_to_remember_[i] = what_to_remember_[i - 1] + i;
	}
	find_DESC_orders();
	// TODO find orders (DAG and DESC)
	find_what_to_remember();
	// TODO what to remember
}

bool general_pattern::avoid(const size_t r, const size_t c, const matrix<size_t>& N)
{
	// I'm not using the fact I know the position which has changed
	// TODO should i become 2k - 1?
	for (size_t i = 0; i < building_tree_.size(); i++) // main loop through added lines (loops 2*k times)
	{
		if (i != 0 && building_tree_[i].size() == 0)
			return true;
		for (size_t m = 0; m < building_tree_[i].size(); m++) // loop through the mappings found in the last iteration
		{
			// find boundaries of added line:
			size_t from, to;
			find_parallel_bounds(i, m, N.getRow(), N.getCol(), from, to);

			size_t l = 0, end = k;
			bool mapped;
			if (orders_[i] < k)
			{
				l = k;
				end = k << 1;
			}
			for (size_t j = from; j <= to; j++) // map orders_[i] to j-th line of N if possible
			{	
				mapped = true;
				for (; l < end; l++)
				{
					if ((lines_[l] >> l) & 1)						// there is a 1 entry at the l-th position of added line in the pattern
					{
						if ((what_to_remember_[i] >> l) & 1)		// I remember the l-th line
						{
							if (!N.at(orders_[i], j))				// and there is no 1-entry at their intersection
							{
								mapped = false;						// so I can't map the line here
								break;
							}
						}
						else
							// TODO find boundaries for j-th line and check if there is enough 1 entries
							;
					}
				}
				if (mapped)
				{
					// extend the mapping - linearly, have to create a new vector, because I might use the current one again for the next line
					std::vector<size_t> extended;
					size_t index = 0;
					for (size_t l = 0; l < (k << 1); l++)
					{
						if ((what_to_remember_[i] >> l) & 1)
						{
							if ((what_to_remember_[i + 1] >> l) & 1)
								extended.push_back(building_tree_[i][m][l]);
							index++;
						}
						else if ((what_to_remember_[i + 1] >> l) & 1)	// this only happens when I need to remember currenly added line
							extended.push_back(j);
					}
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
			}
		}
	}
	if (building_tree_[(k << 1) - 1].empty())
		return true;
	return false;
}

void general_pattern::find_DESC_orders()
{
	// TODO is there an instruction which does this?
	// TODO bucket sort
}

void general_pattern::find_what_to_remember()
{
	// TODO what to remember for given orders
}

void general_pattern::find_parallel_bounds(size_t i, size_t m, size_t rows, size_t columns, size_t& from, size_t& to)
{
	size_t	bot = 0 - 1,	// index to the lower bound for currently added line in map vector
			top = 0 - 1,	// index to the upper bound for currently added line in map vector
			i_top = 0,		// index of the line of the pattern which is a upper bound of currently added line 
			i_bot = 0 - 1;	// index of the line of the pattern which is a lower bound of currently added line
	for (i_top = 0; i_top < orders_.size(); i_top++)
	{
		if ((((orders_[i] << 1) < orders_.size()) && ((i_top << 1) < orders_.size())) ||		// iterating through rows and adding a row
			(((orders_[i] << 1) >= orders_.size()) && ((i_top << 1) >= orders_.size())) &&		// iterating through columns and adding a column
			((what_to_remember_[i] >> i_top) & 1))												// I remember this line
		{
			if (i_top < orders_[i])
			{
				i_bot = i_top;
				bot++;
			}
			else
			{
				top = bot + 1;
				break;
			}
		}
		else if ((what_to_remember_[i] >> i_top) & 1)										// I need to remeber this line, but I'm adding column and going through rows or vice versa
			bot++;
		else if ((((orders_[i] << 1) < orders_.size()) && (((i_top + 1) << 1) == orders_.size())) ||
			(((orders_[i] << 1) >= orders_.size()) && (i_top + 1 == orders_.size())))	// if there is no upper bound
		{
			i_top = 0 - 1;
			break;
		}
	}
	// i have parallel boundaries
	if (bot == 0 - 1)
	{
		if (orders_[i] < k)
			from = 0 + orders_[i];
		else
			from = rows + orders_[i] - k;
	}
	else
		from = building_tree_[i][m][bot] - i_bot + orders_[i];
	if (top == 0 - 1)
	{
		if (orders_[i] < k)
			to = rows - orders_[i];
		else
			to = rows + columns - orders_[i] + k;
	}
	else
		building_tree_[i][m][top] - i_top + orders_[i];
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