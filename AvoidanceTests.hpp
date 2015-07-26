#ifndef AvoidanceTest_hpp_
#define AvoidanceTest_hpp_

#include "Matrix.hpp"
#include <set>

/// For the purposes of complexity, let \theta(k) be the number of lines (rows and columns) of the pattern
/// and \theta(n) be the number of lines of the big matrix.

enum Type { GENERAL, WALKING };

// Enum for the line ordering functions
enum Order { DESC, SUM, MAX, AUTO, CUSTOM };

enum Map { RECURSION, COMPROMISE, NORECURSION };

enum Map_container { VECTOR, SET };

class general_pattern
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
	general_pattern(const matrix<size_t>& pattern, Order order = DESC, Map map_approach = RECURSION, std::vector<size_t>&& custom_order = std::vector<size_t>());
protected:
	Map map_approach_;												// choosen way of mapping algorithm - use recursion for nonmapped lines or not
	size_t	row_,													// number of rows of the pattern
			col_,													// number of columns of the pattern
			steps_,													// number of steps I'm going to do
			empty_lines_;											// binary number of lines with no one-entries
	std::vector<size_t> lines_,										// binary number for each line of a pattern having one at i-th position if the pattern has one-entry there
						order_,										// order of lines in which I am going to be mapping them
						what_to_remember_;							// for each adding line I know which of them I still need to remember for next mapping
	std::vector<std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t> > > > parallel_bound_indices_;
																	// vector through levels - vector through lines - pair of pairs - pair of lower and upper bounds
	std::vector<std::vector<size_t> > extending_order_;				// vector through levels - vector of indices of the mapping which are needed for the extended one
	std::vector<std::vector<size_t> > map_index_;					// vector through levels - vector of indices of lines in the mapping

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
	void find_parallel_bounds(const size_t line, const size_t level, const std::vector<size_t>& mapping, const size_t rows, const size_t columns,
		size_t& from, size_t& to);
	
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
	size_t count_what_to_remember(const size_t current);
	
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
	bool map(const bool backtrack, const size_t line, const size_t level, const size_t big_line, const std::vector<size_t>& mapping, const matrix<size_t>& big_matrix);
	
	/// <summary>
	/// Extends previous mapping after deciding to which big line the line should be mapped.
	/// It is precomputed which elements of previous mapping are needed for the new one and where to put the new element.
	/// Takes k time, for extending.
	/// </summary>
	/// <param name="return">The extended mapping</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="big_line">Index of the line of the big matrix which I mapped the line to.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	std::vector<size_t> extend(const size_t level, const size_t big_line, const std::vector<size_t>& mapping);
};

// general pattern having std::vector as a container for found mappings
class general_vector_pattern
	: public general_pattern
{
public:
	general_vector_pattern(const matrix<size_t>& pattern, Order order = DESC, Map map_approach = RECURSION, std::vector<size_t>&& custom_order = std::vector<size_t>())
		: general_pattern(pattern, order, map_approach, std::move(custom_order)), building_tree_(2)
	{}

	/// <summary>
	/// Tests if the pattern avoids given matrix as a submatrix.
	/// Returns true if it does, false if the matrix contains the pattern.
	/// The program takes one line of the pattern after another and tries to map them to every possible line of the resulting matrix. 
	/// It is, as it sounds, a brute force method (O(n^k)). To make it more efficient, the mappings, which have "important lines" mapped to
	/// the same lines of big matrix are shrinked into one mapping.
	/// </summary>
	/// <param name="N">Matrix for which is tested whether it avoids the pattern.</param>
	bool avoid(const matrix<size_t>& N);
private:
	std::vector<std::vector<std::vector<size_t> > > building_tree_;
};

// general pattern having std::set as a container for found mappings
class general_set_pattern
	: public general_pattern
{
public:
	general_set_pattern(const matrix<size_t>& pattern, Order order = DESC, Map map_approach = RECURSION, std::vector<size_t>&& custom_order = std::vector<size_t>())
		: general_pattern(pattern, order, map_approach, std::move(custom_order)), building_tree_(2)
	{}

	/// <summary>
	/// Tests if the pattern avoids given matrix as a submatrix.
	/// Returns true if it does, false if the matrix contains the pattern.
	/// The program takes one line of the pattern after another and tries to map them to every possible line of the resulting matrix. 
	/// It is, as it sounds, a brute force method (O(n^k)). To make it more efficient, the mappings, which have "important lines" mapped to
	/// the same lines of big matrix are shrinked into one mapping.
	/// </summary>
	/// <param name="N">Matrix for which is tested whether it avoids the pattern.</param>
	bool avoid(const matrix<size_t>& N);
private:
	std::vector<std::set<std::vector<size_t> > > building_tree_;
};


/// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all one-entries.
/// It may contain zero-entries as well, but no one-entry can be not included in the walk.
class walking_pattern
{
public:
	walking_pattern(const matrix<size_t>& pattern, const size_t n);
	
	/// <summary>
	/// Tests if the pattern avoids given matrix as a submatrix.
	/// Returns true if it does, false if the matrix contains the pattern.
	/// Takes n^2 time. It looks at position [r,c], recalculates c_v and c_h (in constant time) and if it changes, recalculates those values
	/// that might be inflicted by the change.
	/// </summary>
	/// <param name="r">Row of the big matrix that has been changed.</param>
	/// <param name="c">Column of the big matrix that has been changed.</param>
	/// <param name="big_matrix">Matrix for which is tested whether it avoids the pattern.</param>
	virtual bool avoid(const size_t r, const size_t c, const matrix<size_t>& big_matrix);
	
	// reverts changes in max_walk_part matrix after an unsuccessful change of the big matrix
	virtual bool revert(const size_t r, const size_t c, const matrix<size_t>& big_matrix) { return avoid(r, c, big_matrix); }
private:
	matrix<std::pair<size_t, size_t> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical) and value of v_i.
	std::vector<size_t> direction_, value_;
};

/// <summary>
/// Computes the number of one-entries of a binary number.
/// Takes constant time and space.
/// </summary>
/// <param name="n">A binary number for which number of bits is computed.</param>
inline size_t bit_count(size_t n)	// I have used a function from the internet: -http://blogs.msdn.com/b/jeuge/archive/2005/06/08/hakmem-bit-count.aspx
{
	size_t uCount = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
	return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

#endif