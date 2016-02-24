#ifndef WalkingPatternFunctions_cpp_
#define WalkingPatternFunctions_cpp_

#include "PatternHeaders.hpp"

#include <queue>
#include <algorithm>
#include <assert.h>

/* Walking pattern */
Walking_pattern::Walking_pattern(const Matrix<size_t>& pattern, const size_t n)
	: max_walk_part_(n, n, std::pair<size_t, size_t>(0, 0))
{
	// the diagonal walk starts in the top left corner
	top_left = true;
	bool not_left = false;

	// coords of the last visited one-entry, index of the column
	size_t last_i = 0, last_j = 0, j;

	// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
	for (size_t sum = 0; sum < pattern.getRow() + pattern.getCol() - 1; ++sum)
	{
		// go through indices of rows
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
					not_left = true;
					goto top_right;	// yeah goto rules
				}

				// need to find, which elements will be a part of the walk
				// from the previously found one-entry go as far to the bottom as you can
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

	///////////////////////////////////////////////////
top_right:

	if (not_left)
	{
		// the diagonal walk starts in the top right corner
		top_left = false;

		// init containers
		direction_.clear();
		value_.clear();

		// coords of the last visited one-entry, index of the column
		long long last_i = 0, last_j = (long long)pattern.getCol() - 1, j;

		// all elements on the same diagonal have the same difference of their coordinates, go through diagonals
		for (long long diff = 1 - (long long)pattern.getCol(); diff < (long long)pattern.getRow(); ++diff)
		{
			// go through indices of rows
			for (long long i = 0; i < (long long)pattern.getRow(); ++i)
			{
				j = i - diff;

				// I look to the left of the pattern
				if (i < diff)
					continue;

				// I look to the right of the pattern
				if (j >= (long long)pattern.getCol())
					break;

				// when I find one-entry or find myself on the last diagonal
				if (pattern.at(i, j) || diff == (long long)pattern.getRow() - 1)
				{
					// last visited element is 0 and I did not find any 1 entries
					if (!pattern.at(i, j) && last_i == 0 && last_j == (long long)pattern.getCol() - 1 && !pattern.at(last_i, last_j)) {
						assert(!"Pattern has no one entries.");
						throw std::invalid_argument("Pattern has no one entries.");
					}

					// I have found a one-entry somewhere where it shouldn't have been
					if (j > last_j || i < last_i) {
						assert(!"Pattern is not a walking pattern type");
						throw std::invalid_argument("Pattern is not a walking pattern type.");
					}

					// need to find, which elements will be a part of the walk
					// from the previously found one-entry go as far to the bottom as you can
					for (long long i2 = last_i; i2 < i; ++i2)
					{
						value_.push_back(pattern.at(i2, last_j));
						// vertical
						direction_.push_back(0);
					}

					// then go to the left and stop right before [i,j]
					for (long long j2 = last_j; j2 > j; --j2)
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
}

bool Walking_pattern::avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& /* sizes */, const size_t r, const size_t c, const size_t& force_end)
{
	typedef std::pair<size_t, size_t> pair;
	std::queue<pair> q;						// queue for elements of the matrix that are supposed to be updated
	pair current;							// [x,y] of the currently updated element
	size_t old_c_v, old_c_h, c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of current, c_h of element to the left

	q.push(pair(r, c));
	while (!q.empty())
	{
		// the function is forced to end from outside
		if (force_end == 1)
			return false;

		current = q.front();
		q.pop();
			
		old_c_v = max_walk_part_.at(current).first;
		old_c_h = max_walk_part_.at(current).second;

		// element on the first row
		if (current.first == 0)		
			c_v_v = 0;
		else
			c_v_v = max_walk_part_.at(current.first - 1, current.second).first;

		if (top_left)
		{
			// element on the first column
			if (current.second == 0)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second - 1).second;
		}
		else
		{
			// element on the last column
			if (current.second == max_walk_part_.getCol() - 1)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second + 1).second;
		}
			
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

			// walk continues to the right/left
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

		if (top_left)
		{
			// c_h was changed and there is still an element to the right
			if (max_walk_part_.at(current).second != old_c_h && current.second + 1 < big_matrix.getCol())
				q.push(pair(current.first, current.second + 1));
		}
		else
		{
			// c_h was changed and there is still an element to the left
			if (max_walk_part_.at(current).second != old_c_h && current.second > 0)
				q.push(pair(current.first, current.second - 1));
		}
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

bool Walking_pattern::parallel_avoid(const size_t /* threads_count */, const Matrix<size_t>& big_matrix, std::vector<Counter>& /* sizes */, const size_t r, const size_t c, const size_t& force_end)
{
	if (top_left)
	{
		// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
		for (size_t sum = r + c; sum < max_walk_part_.getRow() + max_walk_part_.getCol() - 1; ++sum)
		{
			// the function is forced to end from outside
			if (force_end == 1)
				return false;
//#pragma omp parallel for //num_threads(threads_count)
			// go through indices of rows
			for (long long temp = r; temp <= (long long)sum; ++temp)
			{
				size_t i = (size_t)temp;
				size_t j = sum - i;
				size_t c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of i, j, c_h of element to the left

				// I look under the pattern
				if (i >= max_walk_part_.getRow())
					break;

				// I look to the right of the pattern
				if (j >= max_walk_part_.getCol() || j < c)
					continue;

				// element on the first row
				if (i == 0)
					c_v_v = 0;
				else
					c_v_v = max_walk_part_.at(i - 1, j).first;

				// element on the first column
				if (j == 0)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(i, j - 1).second;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).first = c_v_v;
				max_walk_part_.at(i, j).second = c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == value_.size())
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if (max_walk_part_.at(i, j).second < c_v_v + 1)
							max_walk_part_.at(i, j).second = c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if (max_walk_part_.at(i, j).first < c_v_v + 1)
							max_walk_part_.at(i, j).first = c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == value_.size())
						return false;
					if (direction_[c_h_h])
					{
						if (max_walk_part_.at(i, j).second < c_h_h + 1)
							max_walk_part_.at(i, j).second = c_h_h + 1;
					}
					else
					{
						if (max_walk_part_.at(i, j).first < c_h_h + 1)
							max_walk_part_.at(i, j).first = c_h_h + 1;
					}
				}
			}
		}
	}
	else
	{
		// all elements on the same diagonal have the same difference of their coordinates, go through diagonals
		for (long long diff = 1 - (long long)max_walk_part_.getCol(); diff < (long long)max_walk_part_.getRow(); ++diff)
		{
			// the function is forced to end from outside
			if (force_end == 1)
				return false;
//#pragma omp parallel for //num_threads(threads_count)
			// go through indices of rows
			for (long long i = 0; i < (long long)max_walk_part_.getRow(); ++i)
			{
				long long j = i - diff;
				long long c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of i, j, c_h of element to the left

				// I look to the left of the pattern
				if (i < diff)
					continue;

				// I look to the right of the pattern
				if (j >= (long long)max_walk_part_.getCol())
					break;

				// element on the first row
				if (i == 0)
					c_v_v = 0;
				else
					c_v_v = (long long)max_walk_part_.at(i - 1, j).first;


				// element on the last column
				if (j == (long long)max_walk_part_.getCol() - 1)
					c_h_h = 0;
				else
					c_h_h = (long long)max_walk_part_.at(i, j + 1).second;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).first = (size_t)c_v_v;
				max_walk_part_.at(i, j).second = (size_t)c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == (long long)value_.size())
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if ((long long)max_walk_part_.at(i, j).second < c_v_v + 1)
							max_walk_part_.at(i, j).second = (size_t)c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if ((long long)max_walk_part_.at(i, j).first < c_v_v + 1)
							max_walk_part_.at(i, j).first = (size_t)c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == (long long)value_.size())
						return false;
					if (direction_[c_h_h])
					{
						if ((long long)max_walk_part_.at(i, j).second < c_h_h + 1)
							max_walk_part_.at(i, j).second = (size_t)c_h_h + 1;
					}
					else
					{
						if ((long long)max_walk_part_.at(i, j).first < c_h_h + 1)
							max_walk_part_.at(i, j).first = (size_t)c_h_h + 1;
					}
				}
			}
		}
	}
	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

#endif