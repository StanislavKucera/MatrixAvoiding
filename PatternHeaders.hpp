#ifndef PatternHeaders_hpp_
#define PatternHeaders_hpp_

#include "HelpFunctionsAndStructures.hpp"

#include <algorithm>
#include <assert.h>

#include <iostream>

/// For the purposes of complexity, let \theta(k) be the number of lines (rows and columns) of the pattern
/// and \theta(n) be the number of lines of the big matrix.

class Pattern
{
public:
	virtual bool avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& sizes, const size_t r, const size_t c) = 0;
	virtual bool revert(const Matrix<size_t>& big_matrix, const size_t r, const size_t c) = 0;
	virtual std::vector<size_t> get_order() const = 0;
	virtual Pattern* get_new_instance() const = 0;
};

class Slow_pattern
	: public Pattern
{
public:
	Slow_pattern(const Matrix<size_t>& pattern) : one_entries_(0), rows_(pattern.getRow()), cols_(pattern.getCol())
	{
		for (size_t i = 0; i < pattern.getRow(); ++i)
			for (size_t j = 0; j < pattern.getCol(); ++j)
				if (pattern.at(i, j) == 1)
					one_entries_.push_back(std::make_pair(i, j));
	}

	bool avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& /* sizes */, const size_t /* r */ = (size_t)-1, const size_t /* c */ = (size_t)-1)
	{
		done_ = false;
		// goes through all subsets of rows and columns of the right cardinality and tests whether the pattern can be mapped to that subset
		test_all_subsets(0ll, 0ll, rows_, cols_, (long long)big_matrix.getRow(), (long long)big_matrix.getCol(), big_matrix);

		if (done_)
			return false;

		return true;
	}
	bool revert(const Matrix<size_t>& /* big_matrix */, const size_t /* r */, const size_t /* c */) { return true; }
	std::vector<size_t> get_order() const { return std::vector<size_t>(); }
	Pattern* get_new_instance() const { return new Slow_pattern(*this); }
private:
	std::vector<std::pair<size_t, size_t> > one_entries_;	// list of all one entries of the pattern
	const long long rows_, cols_;							// size of the pattern
	bool done_;												// indicator whether the avoidance testing has failed (the matrix does not avoid the pattern)

	void test_all_subsets(long long v_map, long long h_map, long long v_ones, long long h_ones, long long v_vals, long long h_vals, const Matrix<size_t>& big_matrix);
};

template<typename T>
class General_pattern
	: public Pattern
{
public:
	/// <summary>
	/// Constructor of the pattern which stores the lines memory efficiently, identifies empty ones, computes the order of line mapping,
	/// precomputes which lines it needs to remember in each step, how to find parallel bound and ho to extend previous mapping.
	/// </summary>
	/// <param name="pattern">Binary matrix which will form the pattern.</param>
	/// <param name="order">Enum determining which function will be used for line ordering.</param>
	/// <param name="map">Enum determining what conditions will map function check.</param>
	/// <param name="custom_order">Order of lines given by user in case order is set to CUSTOM.</param>
	General_pattern(const Matrix<size_t>& pattern, const Order order = DESC, const Map map_approach = SUPERACTIVE, std::vector<size_t>&& custom_order = std::vector<size_t>());

	/// <summary>
	/// Tests if the pattern avoids given matrix as a submatrix.
	/// Returns true if it does, false if the matrix contains the pattern.
	/// The program takes one line of the pattern after another and tries to map them to every possible line of the resulting matrix. 
	/// It is, as it sounds, a brute force method (O(n^k)). To make it more efficient, the mappings, which have "important lines" mapped to
	/// the same lines of big matrix are shrinked into one mapping.
	/// </summary>
	/// <param name="big_matrix">Matrix for which is tested whether it avoids the pattern.</param>
	/// <param name="r">Row of the big matrix that has been changed.</param>
	/// <param name="c">Column of the big matrix that has been changed.</param>
	/// <param name="sizes">Vector of numbers of found mappings on each level.</param>
	bool avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& sizes, const size_t r = (size_t)-1, const size_t c = (size_t)-1);
	bool revert(const Matrix<size_t>& /* big_matrix */, const size_t /* r */, const size_t /* c */) { return true; }
	std::vector<size_t> get_order() const { return order_; }
	Pattern* get_new_instance() const { return new General_pattern(*this); }
private:
	const size_t	row_,								// number of rows of the pattern
					col_;								// number of columns of the pattern
	std::vector<size_t> lines_,							// binary number for each line of a pattern having one at i-th position if the pattern has one-entry there
														// lines_[i] = (1011)_2 ... i-th line of the pattern has one-enty at 0th, 1st and 3rd position
						order_,							// order of lines in which I am going to be mapping them
														// order_[i] = j ... in i-th step I'm going to map j-th line if the pattern
						what_to_remember_;				// for each adding line I know which of them I still need to remember for next mapping
														// what_to_remember[i] = (001010)_2 ... in i-th step I remember where I mapped the 1st and 3rd line of the pattern
	std::vector<std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t> > > > parallel_bound_indices_;
		// vector through levels - vector through lines - pair of pairs - pair of lower and upper bounds
		// parallel_bound_indices_[i][j] = ((bot, top),(i_bot, i_top)) ... in i-th step, j-th line of the pattern is bounded by i_bot line of the pattern from the bottom
		//	and by i_top line of the pattern from the top; bot and top are indices to mapping structure - mapping[bot] = b ... bot-th line is mapped to b line
	std::vector<std::vector<size_t> > extending_order_;	// vector through levels - vector of indices of the mapping which are needed for the extended one
														// extending_order_[i][j] = k ... in i-th step, j-th linewill be emplaced at k-th position of the mapping
	std::vector<std::vector<size_t> > map_index_;		// vector through levels - vector of indices of lines in the mapping
														// map_index_[i][j] = k ... in i-th step, j-th line is on the k-th position in the mapping
	std::vector<Container<T> > building_tree_;			// container for found mapping at each level
	size_t	steps_,										// number of steps I'm going to do = number of lines I need to map (excluding empty lines)
			empty_lines_;								// binary number of lines with no one-entries
	const Map map_approach_;							// choosen way of mapping algorithm - use recursion for nonmapped lines or not

	/// <summary>
	/// For given line of the pattern computes lines of the big matrix, which bound its mapping.
	/// The bounds are stored to from and to variables meaning that "line" can be mapped to lines [from, to).
	/// Takes constant time since it knows where to look, because indices to the vectors are precalculated.
	/// </summary>
	/// <param name="line">Index of the line of the pattern for which bounds are calculated.</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="rows">Number of rows of the big matrix.</param>
	/// <param name="columns">Number of columns of the big matrix.</param>
	/// <param name="from">Index of the line of the big matrix which bounds "line" from the bottom.</param>
	/// <param name="to">Index of the line of the big matrix which bounds "line" from the top.</param>
	/// <param name="r">Row of the entry that was changed in the last iteration if I know it.</param>
	/// <param name="c">Column of the entry that was changed in the last iteration if I know it.</param>
	void find_parallel_bounds(const size_t line, const size_t level, const std::vector<size_t>& mapping, const size_t rows, const size_t columns,
		size_t& from, size_t& to, const size_t r = (size_t)-1, const size_t c = (size_t)-1) const;
	
	/// <summary>
	/// For given line of the pattern computes lines of the big matrix, which bound its mapping.
	/// The bounds are stored to from and to variables meaning that "line" can be mapped to lines [from, to).
	/// Takes constant time since it knows where to look, because indices to the vectors are precalculated.
	/// </summary>
	/// <param name="line">Index of the line of the pattern for which bounds are calculated.</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="rows">Number of rows of the big matrix.</param>
	/// <param name="columns">Number of columns of the big matrix.</param>
	/// <param name="from">Index of the line of the big matrix which bounds "line" from the bottom.</param>
	/// <param name="to">Index of the line of the big matrix which bounds "line" from the top.</param>
	/// <param name="r">Row of the entry that was changed in the last iteration if I know it.</param>
	/// <param name="c">Column of the entry that was changed in the last iteration if I know it.</param>
	bool check_orthogonal_bounds(const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping,
		const size_t orthogonal_line, const size_t big_orthognal_line, const Matrix<size_t>& big_matrix) const;

	/// <summary>
	/// Orders lines of the pattern according to the number of one-entries descendingly.
	/// Takes k*log(k) time because it needs to sort the lines at the end.
	/// If it encounters a line with no one-entries, it reduces the number of steps of the whole algorithm,
	/// since empty line can be mapped anywhere.
	/// </summary>
	void find_DESC_order();

	/// <summary>
	/// Orders lines of the pattern so that there is the smallest number of lines it needs to remember throughout the whole algorithm.
	/// Smallest in this case means smallest number as a sum of all numbers.
	/// Takes 2^k time because it needs to try all the subsets of lines to find out the best one.
	/// </summary>
	void find_SUM_order();

	/// <summary>
	/// Orders lines of the pattern so that there is the smallest number of lines it needs to remember throughout the whole algorithm.
	/// Smallest in this case means smallest number in the worst case.
	/// Takes 2^k time because it needs to try all the subsets of lines to find out the best one.
	/// </summary>
	void find_MAX_order();
	
	/// <summary>
	/// For given subset returns the number of lines it needs to remember (excluding those, which are not needed).
	/// Takes k^2 time since for each line of the subset it checks whether line - 1 and line + 1 are in the subset
	/// and if they are it checks the same condition for all the lines that intersect the line in a one-entry.
	/// </summary>
	/// <param name="current">Given subset of lines for which I calculate how many lines I need to remember.</param>
	size_t count_what_to_remember(const size_t current) const;
	
	/// <summary>
	/// For given order computes, which already mapped lines need to be stored and which can be forgotten
	/// Takes k^3 time since it k times for each mapped line checks whether line - 1 and line + 1 have been mapped as well
	/// and if they have it checks the same condition for all the lines that intersect the line in a one-entry.
	/// </summary>
	void find_what_to_remember();

	/// <summary>
	/// Precalculates parallel bounds (indices to the mapping, from which I will take the lines) for all lines
	/// for which it makes sense - I will be either adding them or checking them in order to get lower number of mappings.
	/// Takes k^3 time since there is 2*k steps and in which it k times calls find_bound_indices.
	/// </summary>
	void find_parralel_bound_indices();

	/// <summary>
	/// Precalculates parallel bounds (indices to the mapping, from which I will take the lines) for given line.
	/// Takes linear time according to k. It goes through what_to_remeber_ and finds the nearest lines.
	/// </summary>
	/// <param name="line">Given line for which bounds are being precalculated.</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	void find_bound_indices(const size_t line, const size_t level);

	/// <summary>
	/// Precomputes which values of mapping I need to store in the one which is one step forward.
	/// Takes linear time according to k. It goes through the lines in previous step and decides whether to remember them.
	/// </summary>
	void find_extending_order();
	
	/// <summary>
	/// Checks if it is possible to map given line of the pattern to given big_line of the big_matrix.
	/// It goes through line entries and if it finds one-entry, it either checks there is a one-entry in the big_matrix
	/// if the crossing line is already mapped or checks if there is enough one-entries in the big matrix.
	/// Moreover for the crossing line, which is not mapped, it checks if it can be mapped to the line with the found one-entry.
	/// Returns true if the line can be mapped to the big_line and the mapping makes sense according to those lines, which have already been mapped.
	/// Takes up to n^2 time since it goes through the entries of the big_line and of each of them it can go through the crossing line.
	/// </summary>
	/// <param name="backtrack">Indicator whether we want to recursively check mapping possibility for other lines.</param>
	/// <param name="line">Index of the line of the pattern which I am trying to map.</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="big_line">Index of the line of the big matrix which I am trying to map the line to.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="big_matrix">Reference to the big matrix for which I test pattern avoiding.</param>
	bool map(const bool backtrack, const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping, const Matrix<size_t>& big_matrix);
	
	/// <summary>
	/// Extends previous mapping after deciding to which big line the line should be mapped.
	/// It is precomputed which elements of previous mapping are needed for the new one and where to put the new element.
	/// Takes k time, for extending.
	/// </summary>
	/// <param name="return">The extended mapping</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="big_line">Index of the line of the big matrix which I mapped the line to.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	std::vector<size_t> extend(const size_t level, const size_t big_line, const std::vector<size_t>& mapping) const;
};

/// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all one-entries.
/// It may contain zero-entries as well, but no one-entry can be not included in the walk.
class Walking_pattern
	: public Pattern
{
public:
	Walking_pattern(const Matrix<size_t>& pattern, const size_t n);
	
	/// <summary>
	/// Tests if the pattern avoids given matrix as a submatrix.
	/// Returns true if it does, false if the matrix contains the pattern.
	/// Takes n^2 time. It looks at position [r,c], recalculates c_v and c_h (in constant time) and if it changes, recalculates those values
	/// that might be inflicted by the change.
	/// </summary>
	/// <param name="big_matrix">Matrix for which is tested whether it avoids the pattern.</param>
	/// <param name="r">Row of the big matrix that has been changed.</param>
	/// <param name="c">Column of the big matrix that has been changed.</param>
	/// <param name="sizes">Vector of numbers of found mappings on each level.</param>
	bool avoid(const Matrix<size_t>& big_matrix, std::vector<Counter>& sizes, const size_t r = (size_t)-1, const size_t c = (size_t)-1);
	
	// reverts changes in max_walk_part matrix after an unsuccessful change of the big matrix
	bool revert(const Matrix<size_t>& big_matrix, const size_t r, const size_t c) { std::vector<Counter> sizes; return avoid(big_matrix, sizes, r, c); }
	std::vector<size_t> get_order() const { return std::vector<size_t>(); }
	Pattern* get_new_instance() const { return new Walking_pattern(*this); }
private:
	Matrix<std::pair<size_t, size_t> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical) and value of v_i.
	std::vector<size_t> direction_, value_;
	bool top_left;
};

class Patterns
{
public:
	Patterns() : patterns_(0), big_matrix_(), changed_(-1) {}
	Patterns(const Patterns& copy) : big_matrix_(copy.big_matrix_), changed_(-1)
	{
		for (auto& pattern : copy.patterns_)
			patterns_.push_back(pattern->get_new_instance());
	}

	bool avoid(const Matrix<size_t>& big_matrix, std::vector<std::vector<Counter> >& sizes, const size_t r, const size_t c)
	{
		if (changed_ == patterns_.size() - 1)
			changed_ = -1;

		if (changed_ != -1)
		{
			assert(!"A pattern is not in a valid state!");
			throw my_exception("A pattern is not in a valid state!");
		}

		sizes.resize(patterns_.size());

		for (auto& pattern : patterns_)
		{
			++changed_;

			if (!pattern->avoid(big_matrix, sizes[changed_], r, c))
				return false;
		}

		return true;
	}
	bool revert(const Matrix<size_t>& big_matrix, const size_t r, const size_t c)
	{
		for (; changed_ != -1; --changed_)
			if (!patterns_[changed_]->revert(big_matrix, r, c))
			{
				assert(!"Matrix after reverting contains the pattern!");
				throw my_exception("Matrix after reverting contains the pattern!");
			}

		return true;
	}
	std::vector<std::vector<size_t> > get_order() const
	{
		std::vector<std::vector<size_t> > orders(patterns_.size());

		for (size_t i = 0; i < patterns_.size(); ++i)
			orders[i] = std::move(patterns_[i]->get_order());

		return orders;
	}

	bool avoid(std::vector<std::vector<Counter> >& sizes, const size_t r, const size_t c, const size_t& forced_end)
	{
		big_matrix_.flip(r, c);

		long long changed = -1;

		sizes.resize(patterns_.size());

		for (auto& pattern : patterns_)
		{
			++changed;

			// forcing abort from outside or the matrix doesn't avoid the pattern
			if (forced_end || !pattern->avoid(big_matrix_, sizes[changed], r, c))
			{
				// flip the bit back
				big_matrix_.flip(r, c);
				//std::cout << "Avoid [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " (" << forced_end << ")fail" << std::endl;

				// and revert pattern structures if needed
				for (; changed != -1; --changed)
					if (!patterns_[changed]->revert(big_matrix_, r, c))
					{
						assert(!"Matrix after reverting contains the pattern!");
						throw my_exception("Matrix after reverting contains the pattern!");
					}

				return false;
			}
		}

		//std::cout << "Avoid [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " (" << forced_end << ")success" << std::endl;
		changes.push_back(std::make_pair(r, c));
		return true;
	}
	bool revert(const size_t r, const size_t c, const Matrix<size_t>& mat)
	{
		// flip the bit back
		big_matrix_.flip(r, c);
		changes.push_back(std::make_pair(r, c));
		//std::cout << "Revert [" << r << "," << c << "] = " << big_matrix_.at(r, c) << " ()success" << std::endl;

		//std::vector<Counter> sizes;

		// and revert pattern structures if needed
		for (auto& pattern : patterns_)
		{
			if (!pattern->revert(big_matrix_, r, c))
			//if (!pattern->avoid(big_matrix_, sizes, r, c))
			{
				if (check_matrix(mat))
					assert(!"Matrix after reverting contains the pattern!");
				throw my_exception("Matrix after reverting contains the pattern!");
			}
		}

		return true;
	}
	bool check_matrix(const Matrix<size_t>& mat)
	{
		bool diff = false;

		for (size_t i = 0; i != big_matrix_.getRow(); ++i)
			for (size_t j = 0; j != big_matrix_.getCol(); ++j)
				if (mat.at(i, j) != big_matrix_.at(i, j))
				{
					diff = true;
					std::cout << "[" << i << "," << j << "]" << std::endl;
				}

		return diff;
	}

	void add(Pattern* pattern) { patterns_.push_back(pattern); }
	void set_matrix(const Matrix<size_t>& big_matrix) { big_matrix_ = big_matrix; }
private:
	std::vector<Pattern*> patterns_;
	Matrix<size_t> big_matrix_;
	long long changed_;
	std::vector<std::pair<size_t, size_t> > changes;
};

#endif