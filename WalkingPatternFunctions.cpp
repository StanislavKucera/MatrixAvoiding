#ifndef WalkingPatternFunctions_cpp_
#define WalkingPatternFunctions_cpp_

#include "PatternHeaders.hpp"

#include <queue>
#include <algorithm>
#include <assert.h>

/* Walking pattern */
Walking_pattern::Walking_pattern(const Matrix<bool>& pattern, const int n)
	: max_walk_part_(n, n), entries_in_queue_(n, n), direction_(pattern.getRow() + pattern.getCol() - 2), value_(pattern.getRow() + pattern.getCol() - 1),
	entries_to_compute_(2 * n), diagonals_in_queue_(2 * n, false), row_(n), col_(n)
{	
	// the diagonal walk starts in the top left corner
	top_left = true;
	bool not_left = false;
	int index = -1;

	// coords of the last visited one-entry, index of the column
	int last_i = 0, last_j = 0, j;

	// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
	for (int sum = 0; sum < pattern.getRow() + pattern.getCol() - 1; ++sum)
	{
		// go through indices of rows
		for (int i = 0; i <= sum; ++i)	
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
				for (int i2 = last_i; i2 < i; ++i2)
				{
					value_[++index] = pattern.at(i2, last_j);
					// vertical
					direction_[index] = 0;
				}

				// then go to the right and stop right before [i,j]
				for (int j2 = last_j; j2 < j; ++j2)
				{
					value_[++index] = pattern.at(i, j2);
					// horizontal
					direction_[index] = 1;					
				}

				// last element of the walk
				last_i = i; last_j = j;							
			}
		}
	}

	// add the last element of the walk
	value_[++index] = pattern.at(last_i, last_j);

	///////////////////////////////////////////////////
top_right:

	index = -1;

	if (not_left)
	{
		// the diagonal walk starts in the top right corner
		top_left = false;

		// coords of the last visited one-entry, index of the column
		int last_i = 0, last_j = pattern.getCol() - 1, j;

		// all elements on the same diagonal have the same difference of their coordinates, go through diagonals
		for (int diff = 1 - pattern.getCol(); diff < pattern.getRow(); ++diff)
		{
			// go through indices of rows
			for (int i = 0; i < pattern.getRow(); ++i)
			{
				j = i - diff;

				// I look to the left of the pattern
				if (i < diff)
					continue;

				// I look to the right of the pattern
				if (j >= pattern.getCol())
					break;

				// when I find one-entry or find myself on the last diagonal
				if (pattern.at(i, j) || diff == pattern.getRow() - 1)
				{
					// last visited element is 0 and I did not find any 1 entries
					if (!pattern.at(i, j) && last_i == 0 && last_j == pattern.getCol() - 1 && !pattern.at(last_i, last_j)) {
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
					for (int i2 = last_i; i2 < i; ++i2)
					{
						value_[++index] = pattern.at(i2, last_j);
						// vertical
						direction_[index] = 0;
					}

					// then go to the left and stop right before [i,j]
					for (int j2 = last_j; j2 > j; --j2)
					{
						value_[++index] = pattern.at(i, j2);
						// horizontal
						direction_[index] = 1;
					}

					// last element of the walk
					last_i = i; last_j = j;
				}
			}
		}

		// add the last element of the walk
		value_[++index] = pattern.at(last_i, last_j);
	}

	size_ = (int)value_.size();
}

bool Walking_pattern::avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& /* sizes */, const std::atomic<bool>& force_end)
{
	typedef std::pair<int, int> pair;
	std::queue<pair> q;						// queue for elements of the matrix that are supposed to be updated
	pair current;							// [x,y] of the currently updated element
	int old_c_v, old_c_h, c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of current, c_h of element to the left

	q.emplace(r, c);
	while (!q.empty())
	{
		// the function is forced to end from outside
		if (force_end)
			return false;

		current = q.front();
		q.pop();
			
		old_c_v = max_walk_part_.at(current).c_v;
		old_c_h = max_walk_part_.at(current).c_h;

		// element on the first row
		if (current.first == 0)		
			c_v_v = 0;
		else
			c_v_v = max_walk_part_.at(current.first - 1, current.second).c_v;

		if (top_left)
		{
			// element on the first column
			if (current.second == 0)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second - 1).c_h;
		}
		else
		{
			// element on the last column
			if (current.second == col_ - 1)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second + 1).c_h;
		}
			
	// Initialization - copying those already found walks
		max_walk_part_.at(current).c_v = c_v_v;
		max_walk_part_.at(current).c_h = c_h_h;
			
	// Search for longer part of the walk
		// b == 1 or v_{c_v_v + 1} == 0
		if (big_matrix.at(current) || !value_[c_v_v])
		{
			// I found the last element of the walk
			if (c_v_v + 1 == size_)
				return false;

			// walk continues to the right/left
			if (direction_[c_v_v]) 
			{
				if (max_walk_part_.at(current).c_h < c_v_v + 1)
					max_walk_part_.at(current).c_h = c_v_v + 1;
			}
			// walk continues to the bottom
			else 
			{
				if (max_walk_part_.at(current).c_v < c_v_v + 1)
					max_walk_part_.at(current).c_v = c_v_v + 1;
			}
		}

		// N[i,j] == 1
		if (big_matrix.at(current) || !value_[c_h_h]) 
		{
			if (c_h_h + 1 == size_)
				return false;

			if (direction_[c_h_h])
			{
				if (max_walk_part_.at(current).c_h < c_h_h + 1)
					max_walk_part_.at(current).c_h = c_h_h + 1;
			}
			else
			{
				if (max_walk_part_.at(current).c_v < c_h_h + 1)
					max_walk_part_.at(current).c_v = c_h_h + 1;
			}
		}		

		// c_v was changed and there is still an element below the current one
		if (max_walk_part_.at(current).c_v != old_c_v && current.first + 1 < row_)
			// if queue is not empty, check whether the element wasn't added before
			if (q.empty() || q.back() != pair(current.first + 1, current.second))	
				q.emplace(current.first + 1, current.second);

		if (top_left)
		{
			// c_h was changed and there is still an element to the right
			if (max_walk_part_.at(current).c_h != old_c_h && current.second + 1 < col_)
				q.emplace(current.first, current.second + 1);
		}
		else
		{
			// c_h was changed and there is still an element to the left
			if (max_walk_part_.at(current).c_h != old_c_h && current.second > 0)
				q.emplace(current.first, current.second - 1);
		}
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

bool Walking_pattern::initial_avoid(const Matrix<bool>& big_matrix)
{
	if (top_left)
	{
		// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
		for (int sum = 0; sum < max_walk_part_.getRow() + max_walk_part_.getCol() - 1; ++sum)
		{
			// go through indices of rows
			for (int i = 0; i <= sum; ++i)
			{
				int j = sum - i;
				int c_v_v, c_h_h;	// c_v of element to the top of i, j, c_h of element to the left

				// I look under the pattern
				if (i >= max_walk_part_.getRow())
					break;

				// I look to the right of the pattern
				if (j >= max_walk_part_.getCol())
					continue;

				// element on the first row
				if (i == 0)
					c_v_v = 0;
				else
					c_v_v = max_walk_part_.at(i - 1, j).c_v;

				// element on the first column
				if (j == 0)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(i, j - 1).c_h;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).c_v = c_v_v;
				max_walk_part_.at(i, j).c_h = c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == size_)
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if (max_walk_part_.at(i, j).c_h < c_v_v + 1)
							max_walk_part_.at(i, j).c_h = c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_v_v + 1)
							max_walk_part_.at(i, j).c_v = c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == size_)
						return false;

					if (direction_[c_h_h])
					{
						if (max_walk_part_.at(i, j).c_h < c_h_h + 1)
							max_walk_part_.at(i, j).c_h = c_h_h + 1;
					}
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_h_h + 1)
							max_walk_part_.at(i, j).c_v = c_h_h + 1;
					}
				}
			}
		}
	}
	else
	{
		// all elements on the same diagonal have the same difference of their coordinates, go through diagonals
		for (int diff = 1 - max_walk_part_.getCol(); diff < max_walk_part_.getRow(); ++diff)
		{
			// go through indices of rows
			for (int i = 0; i < max_walk_part_.getRow(); ++i)
			{
				int j = i - diff;
				int c_v_v, c_h_h;	// c_v of element to the top of i, j, c_h of element to the left
				
				// I look to the right of the pattern
				if (j >= max_walk_part_.getCol())
					break;

				// I look to the left of the pattern
				if (j < 0)
					continue;

				// element on the first row
				if (i == 0)
					c_v_v = 0;
				else
					c_v_v = max_walk_part_.at(i - 1, j).c_v;

				// element on the last column
				if (j == max_walk_part_.getCol() - 1)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(i, j + 1).c_h;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).c_v = c_v_v;
				max_walk_part_.at(i, j).c_h = c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == size_)
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if (max_walk_part_.at(i, j).c_h < c_v_v + 1)
							max_walk_part_.at(i, j).c_h = c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_v_v + 1)
							max_walk_part_.at(i, j).c_v = c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == size_)
						return false;

					if (direction_[c_h_h])
					{
						if (max_walk_part_.at(i, j).c_h < c_h_h + 1)
							max_walk_part_.at(i, j).c_h = c_h_h + 1;
					}
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_h_h + 1)
							max_walk_part_.at(i, j).c_v = c_h_h + 1;
					}
				}
			}
		}
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

//bool Walking_pattern::lazy_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& /* sizes */, const std::atomic<bool>& force_end)
/*{
	typedef std::pair<int, int> pair;
	pair current, last = std::make_pair(-1, -1);							// [x,y] of the currently updated element
	int old_c_v, old_c_h, c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of current, c_h of element to the left
	
	if (top_left)
		ql_.emplace(r, c);
	else
		qr_.emplace(r, c);

	while (!ql_.empty() || !qr_.empty())
	{
		// the function is forced to end from outside
		if (force_end)
			return false;

		if (top_left)
		{
			current = ql_.top();
			ql_.pop();
		}
		else
		{
			current = qr_.top();
			qr_.pop();
		}

		if (current == last)
			continue;

		last = current;

		old_c_v = max_walk_part_.at(current).c_v;
		old_c_h = max_walk_part_.at(current).c_h;

		// element on the first row
		if (current.first == 0)
			c_v_v = 0;
		else
			c_v_v = max_walk_part_.at(current.first - 1, current.second).c_v;

		if (top_left)
		{
			// element on the first column
			if (current.second == 0)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second - 1).c_h;
		}
		else
		{
			// element on the last column
			if (current.second == max_walk_part_.getCol() - 1)
				c_h_h = 0;
			else
				c_h_h = max_walk_part_.at(current.first, current.second + 1).c_h;
		}

		// Initialization - copying those already found walks
		max_walk_part_.at(current).c_v = c_v_v;
		max_walk_part_.at(current).c_h = c_h_h;

		// Search for longer part of the walk
		// b == 1 or v_{c_v_v + 1} == 0
		if (big_matrix.at(current) || !value_[c_v_v])
		{
			// I found the last element of the walk
			if (c_v_v + 1 == size_)
			{
				if (top_left)
				{
					ql_.emplace(current);

					if (current.first + 1 < big_matrix.getRow())
						ql_.emplace(current.first + 1, current.second);

					if (current.second + 1 < big_matrix.getCol())
						ql_.emplace(current.first, current.second + 1);
				}
				else
				{
					qr_.emplace(current);

					if (current.first + 1 < big_matrix.getRow())
						qr_.emplace(current.first + 1, current.second);

					if (current.second > 0)
						qr_.emplace(current.first, current.second - 1);
				}

				return false;
			}

			// walk continues to the right/left
			if (direction_[c_v_v])
			{
				if (max_walk_part_.at(current).c_h < c_v_v + 1)
					max_walk_part_.at(current).c_h = c_v_v + 1;
			}
			// walk continues to the bottom
			else
			{
				if (max_walk_part_.at(current).c_v < c_v_v + 1)
					max_walk_part_.at(current).c_v = c_v_v + 1;
			}
		}

		// N[i,j] == 1
		if (big_matrix.at(current) || !value_[c_h_h])
		{
			if (c_h_h + 1 == size_)
			{
				if (top_left)
				{
					ql_.emplace(current);

					if (current.first + 1 < big_matrix.getRow())
						ql_.emplace(current.first + 1, current.second);

					if (current.second + 1 < big_matrix.getCol())
						ql_.emplace(current.first, current.second + 1);
				}
				else
				{
					qr_.emplace(current);

					if (current.first + 1 < big_matrix.getRow())
						qr_.emplace(current.first + 1, current.second);
					
					if (current.second > 0)
						qr_.emplace(current.first, current.second - 1);
				}

				return false;
			}


			if (direction_[c_h_h])
			{
				if (max_walk_part_.at(current).c_h < c_h_h + 1)
					max_walk_part_.at(current).c_h = c_h_h + 1;
			}
			else
			{
				if (max_walk_part_.at(current).c_v < c_h_h + 1)
					max_walk_part_.at(current).c_v = c_h_h + 1;
			}
		}		

		if (top_left)
		{
			// c_v was changed and there is still an element below the current one
			if (max_walk_part_.at(current).c_v != old_c_v && current.first + 1 < big_matrix.getRow())
				ql_.emplace(current.first + 1, current.second);
		}
		else
		{
			// c_v was changed and there is still an element below the current one
			if (max_walk_part_.at(current).c_v != old_c_v && current.first + 1 < big_matrix.getRow())
				qr_.emplace(current.first + 1, current.second);
		}

		if (top_left)
		{
			// c_h was changed and there is still an element to the right
			if (max_walk_part_.at(current).c_h != old_c_h && current.second + 1 < big_matrix.getCol())
				ql_.emplace(current.first, current.second + 1);
		}
		else
		{
			// c_h was changed and there is still an element to the left
			if (max_walk_part_.at(current).c_h != old_c_h && current.second > 0)
				qr_.emplace(current.first, current.second - 1);
		}
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}*/

bool Walking_pattern::lazy_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& /* sizes */, const std::atomic<bool>& force_end)
{
	std::pair<int, int> current;		// [x,y] of the currently updated element
	int old_c_v, old_c_h, c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of current, c_h of element to the left

	// add [r, c] to the structures if needed
	if (top_left)
	{
		if (!diagonals_in_queue_[r + c])
		{
			diagonals_in_queue_[r + c] = true;
			diagonals_to_compute_.emplace(r + c);
		}

		if (!entries_in_queue_.at(r, c))
		{
			entries_in_queue_.at(r, c) = true;
			entries_to_compute_[r + c].emplace_back(r, c);
		}
	}
	else
	{
		if (!diagonals_in_queue_[r - c + col_])
		{
			diagonals_in_queue_[r - c + col_] = true;
			diagonals_to_compute_.emplace(r - c + col_);
		}

		if (!entries_in_queue_.at(r, c))
		{
			entries_in_queue_.at(r, c) = true;
			entries_to_compute_[r - c + col_].emplace_back(r, c);
		}
	}
	
	while (!diagonals_to_compute_.empty())
	{
		const int diagonal = diagonals_to_compute_.top();

		while (!entries_to_compute_[diagonal].empty())
		{
			// the function is forced to end from outside
			if (force_end)
				return false;

			current = entries_to_compute_[diagonal].back();

			old_c_v = max_walk_part_.at(current).c_v;
			old_c_h = max_walk_part_.at(current).c_h;

			// element on the first row
			if (current.first == 0)
				c_v_v = 0;
			else
				c_v_v = max_walk_part_.at(current.first - 1, current.second).c_v;

			if (top_left)
			{
				// element on the first column
				if (current.second == 0)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(current.first, current.second - 1).c_h;
			}
			else
			{
				// element on the last column
				if (current.second == col_ - 1)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(current.first, current.second + 1).c_h;
			}

			// Initialization - copying those already found walks
			max_walk_part_.at(current).c_v = c_v_v;
			max_walk_part_.at(current).c_h = c_h_h;

			// Search for longer part of the walk
			// b == 1 or v_{c_v_v + 1} == 0
			if (big_matrix.at(current) || !value_[c_v_v])
			{
				// I found the last element of the walk
				if (c_v_v + 1 == size_)
				{
					const int r_ = current.first;
					const int c_ = current.second;

					if (top_left)
					{
						if (r_ + c_ + 1 < row_ + col_ - 1 && !diagonals_in_queue_[r_ + c_ + 1])
						{
							diagonals_in_queue_[r_ + c_ + 1] = true;
							diagonals_to_compute_.emplace(r_ + c_ + 1);
						}

						if (r_ + 1 < row_ && !entries_in_queue_.at(r_ + 1, c_))
						{
							entries_in_queue_.at(r_ + 1, c_) = true;
							entries_to_compute_[r_ + c_ + 1].emplace_back(r_ + 1, c_);
						}

						if (c_ + 1 < col_ && !entries_in_queue_.at(r_, c_ + 1))
						{
							entries_in_queue_.at(r_, c_ + 1) = true;
							entries_to_compute_[r_ + c_ + 1].emplace_back(r_, c_ + 1);
						}
					}
					else
					{
						if (r_ - c_ + 1 < row_ && !diagonals_in_queue_[r_ - c_ + col_ + 1])
						{
							diagonals_in_queue_[r_ - c_ + col_ + 1] = true;
							diagonals_to_compute_.emplace(r_ - c_ + col_ + 1);
						}

						if (r_ + 1 < row_ && !entries_in_queue_.at(r_ + 1, c_))
						{
							entries_in_queue_.at(r_ + 1, c_) = true;
							entries_to_compute_[r_ - c_ + col_ + 1].emplace_back(r_ + 1, c_);
						}

						if (c_ > 0 && !entries_in_queue_.at(r_, c_ - 1))
						{
							entries_in_queue_.at(r_, c_ - 1) = true;
							entries_to_compute_[r_ - c_ + col_ + 1].emplace_back(r_, c_ - 1);
						}
					}

					return false;
				}

				// walk continues to the right/left
				if (direction_[c_v_v])
				{
					if (max_walk_part_.at(current).c_h < c_v_v + 1)
						max_walk_part_.at(current).c_h = c_v_v + 1;
				}
				// walk continues to the bottom
				else
				{
					if (max_walk_part_.at(current).c_v < c_v_v + 1)
						max_walk_part_.at(current).c_v = c_v_v + 1;
				}
			}

			// N[i,j] == 1
			if (big_matrix.at(current) || !value_[c_h_h])
			{
				if (c_h_h + 1 == size_)
				{
					const int r_ = current.first;
					const int c_ = current.second;

					if (top_left)
					{
						if (r_ + c_ + 1 < row_ + col_ - 1 && !diagonals_in_queue_[r_ + c_ + 1])
						{
							diagonals_in_queue_[r_ + c_ + 1] = true;
							diagonals_to_compute_.emplace(r_ + c_ + 1);
						}

						if (r_ + 1 < row_ && !entries_in_queue_.at(r_ + 1, c_))
						{
							entries_in_queue_.at(r_ + 1, c_) = true;
							entries_to_compute_[r_ + c_ + 1].emplace_back(r_ + 1, c_);
						}

						if (c_ + 1 < col_ && !entries_in_queue_.at(r_, c_ + 1))
						{
							entries_in_queue_.at(r_, c_ + 1) = true;
							entries_to_compute_[r_ + c_ + 1].emplace_back(r_, c_ + 1);
						}
					}
					else
					{
						if (r_ - c_ + 1 < row_ && !diagonals_in_queue_[r_ - c_ + col_ + 1])
						{
							diagonals_in_queue_[r_ - c_ + col_ + 1] = true;
							diagonals_to_compute_.emplace(r_ - c_ + col_ + 1);
						}

						if (r_ + 1 < row_ && !entries_in_queue_.at(r_ + 1, c_))
						{
							entries_in_queue_.at(r_ + 1, c_) = true;
							entries_to_compute_[r_ - c_ + col_ + 1].emplace_back(r_ + 1, c_);
						}

						if (c_ > 0 && !entries_in_queue_.at(r_, c_ - 1))
						{
							entries_in_queue_.at(r_, c_ - 1) = true;
							entries_to_compute_[r_ - c_ + col_ + 1].emplace_back(r_, c_ - 1);
						}
					}

					return false;
				}

				if (direction_[c_h_h])
				{
					if (max_walk_part_.at(current).c_h < c_h_h + 1)
						max_walk_part_.at(current).c_h = c_h_h + 1;
				}
				else
				{
					if (max_walk_part_.at(current).c_v < c_h_h + 1)
						max_walk_part_.at(current).c_v = c_h_h + 1;
				}
			}

			if (top_left)
			{
				// c_v was changed and there is still an element below the current one
				if (current.first + 1 < row_ && max_walk_part_.at(current).c_v != old_c_v && !entries_in_queue_.at(current.first + 1, current.second))
				{
					entries_in_queue_.at(current.first + 1, current.second) = true;
					entries_to_compute_[current.first + current.second + 1].emplace_back(current.first + 1, current.second);

					if (!diagonals_in_queue_[current.first + current.second + 1])
					{
						diagonals_in_queue_[current.first + current.second + 1] = true;
						diagonals_to_compute_.emplace(current.first + current.second + 1);
					}
				}

				// c_h was changed and there is still an element to the right
				if (current.second + 1 < col_ && max_walk_part_.at(current).c_h != old_c_h && !entries_in_queue_.at(current.first, current.second + 1))
				{
					entries_in_queue_.at(current.first, current.second + 1) = true;
					entries_to_compute_[current.first + current.second + 1].emplace_back(current.first, current.second + 1);

					if (!diagonals_in_queue_[current.first + current.second + 1])
					{
						diagonals_in_queue_[current.first + current.second + 1] = true;
						diagonals_to_compute_.emplace(current.first + current.second + 1);
					}
				}
			}
			else
			{
				// c_v was changed and there is still an element below the current one
				if (current.first + 1 < row_ && max_walk_part_.at(current).c_v != old_c_v && !entries_in_queue_.at(current.first + 1, current.second))
				{
					entries_in_queue_.at(current.first + 1, current.second) = true;
					entries_to_compute_[current.first - current.second + col_ + 1].emplace_back(current.first + 1, current.second);

					if (!diagonals_in_queue_[current.first - current.second + col_ + 1])
					{
						diagonals_in_queue_[current.first - current.second + col_ + 1] = true;
						diagonals_to_compute_.emplace(current.first - current.second + col_ + 1);
					}
				}

				// c_h was changed and there is still an element to the left
				if (current.second > 0 && max_walk_part_.at(current).c_h != old_c_h && !entries_in_queue_.at(current.first, current.second - 1))
				{
					entries_in_queue_.at(current.first, current.second - 1) = true;
					entries_to_compute_[current.first - current.second + col_ + 1].emplace_back(current.first, current.second - 1);

					if (!diagonals_in_queue_[current.first - current.second + col_ + 1])
					{
						diagonals_in_queue_[current.first - current.second + col_ + 1] = true;
						diagonals_to_compute_.emplace(current.first - current.second + col_ + 1);
					}
				}
			}

			entries_to_compute_[diagonal].pop_back();
			entries_in_queue_.at(current) = false;
		}

		diagonals_to_compute_.pop();
		diagonals_in_queue_[diagonal] = false;
	}

	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

bool Walking_pattern::parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& /* sizes */, const int /* threads_count */, const std::atomic<bool>& force_end)
{
	// r and c are always valid in parallel_avoid

	if (top_left)
	{
		// all elements on the same diagonal have the same sum of their coordinates, go through diagonals
		for (int sum = r + c; sum < max_walk_part_.getRow() + max_walk_part_.getCol() - 1; ++sum)
		{
			// the function is forced to end from outside
			if (force_end)
				return false;

//#pragma omp parallel for //num_threads(threads_count)
			// go through indices of rows
			for (int i = r; i <= sum; ++i)
			{
				int j = sum - i;
				int c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of i, j, c_h of element to the left

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
					c_v_v = max_walk_part_.at(i - 1, j).c_v;

				// element on the first column
				if (j == 0)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(i, j - 1).c_h;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).c_v = c_v_v;
				max_walk_part_.at(i, j).c_h = c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == size_)
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if (max_walk_part_.at(i, j).c_h < c_v_v + 1)
							max_walk_part_.at(i, j).c_h = c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_v_v + 1)
							max_walk_part_.at(i, j).c_v = c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == size_)
						return false;

					if (direction_[c_h_h])
					{
						if (max_walk_part_.at(i, j).c_h < c_h_h + 1)
							max_walk_part_.at(i, j).c_h = c_h_h + 1;
					}
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_h_h + 1)
							max_walk_part_.at(i, j).c_v = c_h_h + 1;
					}
				}
			}
		}
	}
	else
	{
		// all elements on the same diagonal have the same difference of their coordinates, go through diagonals
		for (int diff = r - c; diff < max_walk_part_.getRow(); ++diff)
		{
			// the function is forced to end from outside
			if (force_end)
				return false;

//#pragma omp parallel for //num_threads(threads_count)
			// go through indices of rows
			for (int i = 0; i < max_walk_part_.getRow(); ++i)
			{
				int j = i - diff;
				int c_v_v, c_h_h;	// c_v and c_h before an update, c_v of element to the top of i, j, c_h of element to the left

				// I look to the left of the pattern
				if (i < diff)
					continue;

				// I look to the right of the pattern
				if (j >= max_walk_part_.getCol())
					break;

				// element on the first row
				if (i == 0)
					c_v_v = 0;
				else
					c_v_v = max_walk_part_.at(i - 1, j).c_v;


				// element on the last column
				if (j == max_walk_part_.getCol() - 1)
					c_h_h = 0;
				else
					c_h_h = max_walk_part_.at(i, j + 1).c_h;

				// Initialization - copying those already found walks
				max_walk_part_.at(i, j).c_v = c_v_v;
				max_walk_part_.at(i, j).c_h = c_h_h;

				// Search for longer part of the walk
				// b == 1 or v_{c_v_v + 1} == 0
				if (big_matrix.at(i, j) || !value_[c_v_v])
				{
					// I found the last element of the walk
					if (c_v_v + 1 == size_)
						return false;

					// walk continues to the right/left
					if (direction_[c_v_v])
					{
						if (max_walk_part_.at(i, j).c_h < c_v_v + 1)
							max_walk_part_.at(i, j).c_h = c_v_v + 1;
					}
					// walk continues to the bottom
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_v_v + 1)
							max_walk_part_.at(i, j).c_v = c_v_v + 1;
					}
				}

				// N[i,j] == 1
				if (big_matrix.at(i, j) || !value_[c_h_h])
				{
					if (c_h_h + 1 == size_)
						return false;

					if (direction_[c_h_h])
					{
						if (max_walk_part_.at(i, j).c_h < c_h_h + 1)
							max_walk_part_.at(i, j).c_h = c_h_h + 1;
					}
					else
					{
						if (max_walk_part_.at(i, j).c_v < c_h_h + 1)
							max_walk_part_.at(i, j).c_v = c_h_h + 1;
					}
				}
			}
		}
	}
	// I haven't mapped the last element of the walk - matrix avoids the pattern
	return true;
}

#endif