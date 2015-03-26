#ifndef AvoidanceTest_cpp_
#define AvoidanceTest_cpp_

#include "AvoidanceTests.hpp"
#include <queue>

/* Walking pattern */
	walking_pattern::walking_pattern(const matrix<int>& pattern, const size_t n)
		: max_walk_part_(n, n, std::pair<int, int>(0, 0))
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

	bool walking_pattern::avoid(const size_t r, const size_t c, const matrix<int>& N)
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