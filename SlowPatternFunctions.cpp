#ifndef SlowPatternFunctions_cpp_
#define SlowPatternFunctions_cpp_

#include "PatternHeaders.hpp"

void Slow_pattern::test_all_subsets(long long v_map, long long h_map, long long v_ones, long long h_ones, long long v_vals, long long h_vals, const Matrix<size_t>& big_matrix)
{
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
		for (size_t i = 0; i < v_vals; ++i)
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
		for (size_t i = 0; i < h_vals; ++i)
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

		for (size_t i = 0; i < big_matrix.getRow(); ++i)
		{
			if (v_map & 1)
				rows.push_back(i);
			v_map /= 2;
		}

		for (size_t i = 0; i < big_matrix.getCol(); ++i)
		{
			if (h_map & 1)
				cols.push_back(i);
			h_map /= 2;
		}

		for (auto& one : one_entries_)
		{
			if (big_matrix.at(rows[one.first], cols[one.second]) == 0)
				return;
		}

		done_ = true;
		return;
	}

	// I either extend the subset by 1 or by 0 or I don't change it -> 3^2 possibilities - 1 because I won't call myself on the same intance
	if (v_ones > 0 && h_ones > 0 && v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2 + 1, h_map * 2 + 1, v_ones - 1, h_ones - 1, v_vals - 1, h_vals - 1, big_matrix);
	if (v_vals > v_ones && h_ones > 0 && h_vals > h_ones)
		test_all_subsets(v_map * 2, h_map * 2 + 1, v_ones, h_ones - 1, v_vals - 1, h_vals - 1, big_matrix);
	if (v_ones > 0 && v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2 + 1, h_map * 2, v_ones - 1, h_ones, v_vals - 1, h_vals - 1, big_matrix);
	if (v_vals > v_ones && h_vals > h_ones)
		test_all_subsets(v_map * 2, h_map * 2, v_ones, h_ones, v_vals - 1, h_vals - 1, big_matrix);
	if (h_vals == 0 && v_ones > 0 && v_vals > v_ones)
		test_all_subsets(v_map * 2 + 1, h_map, v_ones - 1, h_ones, v_vals - 1, h_vals, big_matrix);
	if (h_vals == 0 && v_vals > v_ones)
		test_all_subsets(v_map * 2, h_map, v_ones, h_ones, v_vals - 1, h_vals, big_matrix);
	if (v_vals == 0 && h_ones > 0 && h_vals > h_ones)
		test_all_subsets(v_map, h_map * 2 + 1, v_ones, h_ones - 1, v_vals, h_vals - 1, big_matrix);
	if (v_vals == 0 && h_vals > h_ones)
		test_all_subsets(v_map, h_map * 2, v_ones, h_ones, v_vals, h_vals - 1, big_matrix);
}

#endif