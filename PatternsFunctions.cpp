#ifndef PatternsFunctions_cpp_
#define PatternsFunctions_cpp_

#include "PatternHeaders.hpp"

bool Patterns::avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<std::vector<Counter> >& sizes)
{
	if (changed_ == (int)patterns_.size() - 1)
		changed_ = -1;

	if (changed_ != -1)
		throw my_exception("A pattern is not in a valid state!");

	sizes.resize(patterns_.size());

	for (auto& pattern : patterns_)
	{
		++changed_;

		if (!pattern->avoid(big_matrix, r, c, sizes[changed_]))
		{
			sizes.resize(changed_ + 1);
			return false;
		}
	}

	return true;
}

bool Patterns::revert(const Matrix<bool>& big_matrix, const int r, const int c)
{
	for (; changed_ != -1; --changed_)
		if (!patterns_[changed_]->revert(big_matrix, r, c))
			throw my_exception("Matrix after reverting contains the pattern!");

	return true;
}

bool Patterns::parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<std::vector<Counter> >& sizes, const int threads_count)
{
	if (changed_ == (int)patterns_.size() - 1)
		changed_ = -1;

	if (changed_ != -1)
		throw my_exception("A pattern is not in a valid state!");

	sizes.resize(patterns_.size());

	for (auto& pattern : patterns_)
	{
		++changed_;

		if (!pattern->parallel_avoid(big_matrix, r, c, sizes[changed_], threads_count))
		{
			sizes.resize(changed_ + 1);
			return false;
		}
	}

	return true;
}

bool Patterns::parallel_revert(const Matrix<bool>& big_matrix, const int r, const int c, const int threads_count)
{
	for (; changed_ != -1; --changed_)
		if (!patterns_[changed_]->parallel_revert(big_matrix, r, c, threads_count))
			throw my_exception("Matrix after reverting contains the pattern!");

	return true;
}

std::vector<std::vector<int> > Patterns::get_order() const
{
	std::vector<std::vector<int> > orders(patterns_.size());

	for (size_t i = 0; i < patterns_.size(); ++i)
		orders[i] = std::move(patterns_[i]->get_order());

	return orders;
}

bool Patterns::avoid(const int r, const int c, std::vector<std::vector<Counter> >& sizes, const std::atomic_bool& forced_end)
{
	big_matrix_.flip(r, c);

	//////////////////////// why don't I use changed_?
	int changed = -1;

	sizes.resize(patterns_.size());

	for (auto& pattern : patterns_)
	{
		++changed;

		// forcing abort from outside or the matrix doesn't avoid the pattern
		if (forced_end || !pattern->avoid(big_matrix_, r, c, sizes[changed], forced_end))
		{
			// flip the bit back
			big_matrix_.flip(r, c);
			//std::cout << "Avoid [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " (" << forced_end << ")fail" << std::endl;

			// and revert pattern structures if needed
			for (; changed != -1; --changed)
				if (!patterns_[changed]->revert(big_matrix_, r, c))
					throw my_exception("Matrix after reverting contains the pattern!");

			return false;
		}
	}

	//std::cout << "Avoid [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " (" << forced_end << ")success" << std::endl;
	//changes.push_back(std::make_pair(r, c));
	return true;
}

bool Patterns::revert(const int r, const int c)//, const Matrix<bool>& mat)
{
	// flip the bit back
	big_matrix_.flip(r, c);
	//changes.push_back(std::make_pair(r, c));
	//std::cout << "Revert [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " ()success" << std::endl;

	//std::vector<Counter> sizes;

	// and revert pattern structures if needed
	for (auto& pattern : patterns_)
	{
		if (!pattern->revert(big_matrix_, r, c)) {
			//if (!pattern->avoid(big_matrix_, sizes, r, c))
			//if (check_matrix(mat))
				throw my_exception("Matrix after reverting contains the pattern!");
			throw my_exception("Matrix after reverting contains the pattern!");
		}
	}

	return true;
}

bool Patterns::check_matrix(const Matrix<bool>& mat)
{
	bool diff = false;

	for (int i = 0; i != big_matrix_.getRow(); ++i)
		for (int j = 0; j != big_matrix_.getCol(); ++j)
			if (mat.at(i, j) != big_matrix_.at(i, j))
				diff = true;

	return diff;
}

void Patterns::construct_threads(const Matrix<bool>& big_matrix)
{
	for (auto& pattern : patterns_)
		pattern->construct_threads(big_matrix);
}

void Patterns::destruct_threads()
{
	for (auto& pattern : patterns_)
		pattern->destruct_threads();
}

#endif