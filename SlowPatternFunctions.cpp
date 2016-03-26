#ifndef SlowPatternFunctions_cpp_
#define SlowPatternFunctions_cpp_

#include "PatternHeaders.hpp"

bool Slow_pattern::parallel_avoid(const Matrix<bool>& big_matrix, const int /* r */, const int /* c */, std::vector<Counter>& /* sizes */, const int /* threads_count */, const std::atomic<bool>& force_end)
{
	//////////////////////////
	// not parallel for now //
	//////////////////////////

	done_ = false;
	// goes through all subsets of rows and columns of the right cardinality and tests whether the pattern can be mapped to that subset
	test_all_subsets(0ll, 0ll, rows_, cols_, (int)big_matrix.getRow(), (int)big_matrix.getCol(), big_matrix, force_end);

	if (done_)
		return false;

	return true;
}

void Slow_pattern::test_all_subsets(int v_map, int h_map, int v_ones, int h_ones, int v_vals, int h_vals, const Matrix<bool>& big_matrix, const std::atomic<bool>& force_end)
{
	// the function is forced to end from outside
	if (force_end)
		done_ = true;

	// I have already found a mapping
	if (done_)
		return;

	// There is no hope - I have used too many one-entries or the set is too big, or there is less lines left than one to be used
	if (v_ones < 0 || h_ones < 0 || v_vals < 0 || h_vals < 0 || v_ones > v_vals || h_ones > h_vals)
		return;

	// There is enough one-entries already, so I fill the rest with zeros
	if (v_ones == 0 && v_vals > 0)
	{
		v_map <<= v_vals;
		v_vals = 0;
	}
	// There is exactly as many one-entries missing as how many entries I still need to add
	else if (v_ones == v_vals && v_ones != 0)
	{
		for (int i = 0; i < v_vals; ++i)
		{
			v_map *= 2;
			v_map += 1;
		}

		v_vals = 0;
		v_ones = 0;
	}
	
	// There is enough one-entries already, so I fill the rest with zeros
	if (h_ones == 0 && h_vals > 0)
	{
		h_map <<= h_vals;
		h_vals = 0;
	}
	// There is exactly as many one-entries missing as how many entries I still need to add
	else if (h_ones == h_vals && h_ones != 0)
	{
		for (int i = 0; i < h_vals; ++i)
		{
			h_map *= 2;
			h_map += 1;
		}

		h_vals = 0;
		h_ones = 0;
	}

	// The subset of rows and columns was chosen
	if (v_vals == 0 && v_ones == 0 && h_vals == 0 && h_ones == 0)
	{
		std::vector<size_t> rows, cols;

		for (int i = 0; i < big_matrix.getRow(); ++i)
		{
			if (v_map & 1)
				rows.emplace_back(i);
			v_map /= 2;
		}

		for (int i = 0; i < big_matrix.getCol(); ++i)
		{
			if (h_map & 1)
				cols.emplace_back(i);
			h_map /= 2;
		}

		for (auto& one : one_entries_)
		{
			if (!big_matrix.at(rows[one.first], cols[one.second]))
				return;
		}

		done_ = true;
		return;
	}

	// I either extend the subset by 1 or by 0 or I don't change it -> 3^2 possibilities - 1 because I won't call myself on the same intance
	if (v_ones > 0 && h_ones > 0 && v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2 + 1, h_map * 2 + 1, v_ones - 1, h_ones - 1, v_vals - 1, h_vals - 1, big_matrix, force_end);
	if (v_vals > v_ones && h_ones > 0 && h_vals > h_ones)
		test_all_subsets(v_map * 2, h_map * 2 + 1, v_ones, h_ones - 1, v_vals - 1, h_vals - 1, big_matrix, force_end);
	if (v_ones > 0 && v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2 + 1, h_map * 2, v_ones - 1, h_ones, v_vals - 1, h_vals - 1, big_matrix, force_end);
	if (v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2, h_map * 2, v_ones, h_ones, v_vals - 1, h_vals - 1, big_matrix, force_end);
	if (h_vals == 0 && v_ones > 0 && v_vals > v_ones)
		test_all_subsets(v_map * 2 + 1, h_map, v_ones - 1, h_ones, v_vals - 1, h_vals, big_matrix, force_end);
	if (h_vals == 0 && v_vals > v_ones)
		test_all_subsets(v_map * 2, h_map, v_ones, h_ones, v_vals - 1, h_vals, big_matrix, force_end);
	if (v_vals == 0 && h_ones > 0 && h_vals > h_ones)
		test_all_subsets(v_map, h_map * 2 + 1, v_ones, h_ones - 1, v_vals, h_vals - 1, big_matrix, force_end);
	if (v_vals == 0 && h_vals > h_ones)
		test_all_subsets(v_map, h_map * 2, v_ones, h_ones, v_vals, h_vals - 1, big_matrix, force_end);
}

#endif