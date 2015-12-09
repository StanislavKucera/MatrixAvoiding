#ifndef WalkingPatternFunctions_cpp_
#define WalkingPatternFunctions_cpp_

#include "PatternHeaders.hpp"

#include <queue>
#include <algorithm>
#include <assert.h>

/* Walking pattern */
Walking_pattern::Walking_pattern(const Matrix<size_t>& pattern, size_t n)
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

bool Walking_pattern::avoid(const Matrix<size_t>& big_matrix, std::vector<std::pair<std::pair<size_t, size_t>, size_t> >& /* sizes */, size_t r, size_t c)
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