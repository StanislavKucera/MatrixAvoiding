#ifndef AvoidanceTest_hpp_
#define AvoidanceTest_hpp_

#include "Matrix.hpp"

// General matrix pattern
class general_pattern
{
public:
	explicit general_pattern(const matrix<size_t>& pattern);
	// returns true if matrix N avoids general pattern this as a submatrix
	bool avoid(const matrix<size_t>& N);
private:
	size_t	k,														// size of a matrix (k x k)
			steps;													// number of steps I'm going to do
	std::vector<size_t> lines_,										// binary number for each line of a pattern having one at i-th position if the pattern has one-entry there
						order_,									// order of lines in which I am going to be mapping them
						what_to_remember_;							// for each adding line I know which of them I still need to remember for next mapping
	std::vector<std::vector<std::vector<size_t> > > building_tree_;	// vector through layers - vector through mappings on each layer - vector through indices of mapped lines

	// for given index of mapped line returns (last two arguments) indices of big matrix which line can be mapped in [from, to)
	void find_parallel_bounds(const size_t, const size_t, const size_t, const size_t, const size_t, size_t&, size_t&, size_t&);
	// orders lines of the pattern according to the number of one-entries descendingly
	void find_DESC_order();
	void find_DAG_order();
	// for given order computes, which already mapped lines need to be stored and which can be forgotten
	void find_what_to_remember();
	// returns true if given line of the pattern can be mapped into a given line of the big matrix
	bool map(const bool, const size_t, const size_t, const size_t, const size_t, const matrix<size_t>&);
	// adds given mapping to the vector of all mappings if it is not already in there
	void extend(const size_t, const size_t, const size_t);
};

// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all 1 entries (may contain 0 entries too)
class walking_pattern
{
public:
	walking_pattern(const matrix<size_t>& pattern, const size_t n);
	// returns true if matrix N avoids walking pattern, on which function is called. Using already computed table by updating it from [r,c] and beyond.
	virtual bool avoid(const size_t r, const size_t c, const matrix<size_t>& N);
	// reverts changes in max_walk_part matrix after an unsuccessful change of the big matrix
	virtual bool revert(const size_t r, const size_t c, const matrix<size_t>& N) { return avoid(r, c, N); }
private:
	matrix<std::pair<size_t, size_t> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical) and value of v_i.
	std::vector<size_t> direction_, value_;
};

// computes the number of one-entries of a binary number in constant time and space
inline size_t bit_count(size_t n)	// I have used a function from the internet: -http://blogs.msdn.com/b/jeuge/archive/2005/06/08/hakmem-bit-count.aspx
{
	size_t uCount = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
	return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

#endif