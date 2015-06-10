#ifndef AvoidanceTest_hpp_
#define AvoidanceTest_hpp_

#include "Matrix.hpp"

class grandfather_pattern
{
public:
	virtual bool avoid(const size_t r, const size_t c, const matrix<size_t>& N) { return 0; }
	virtual bool revert(const size_t r, const size_t c, const matrix<size_t>& N) { return 0; }
};

// General matrix pattern
class general_pattern
	: public grandfather_pattern
{
public:
	explicit general_pattern(const matrix<size_t>& pattern, const size_t n);

	virtual bool avoid(const size_t r, const size_t c, const matrix<size_t>& N);
private:
	size_t k;
	std::vector<size_t> rows_, cols_, lines_, orders_, what_to_remember_;
	std::vector<std::vector<std::vector<size_t> > > building_tree_;		// vector through layers - vector through mappings on each layer - vector through indices of mapped lines

	void find_parallel_bounds(size_t i, size_t m, size_t rows, size_t columns, size_t& from, size_t& to);
	void find_DESC_orders();
	void find_what_to_remember();
};

// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all 1 entries (may contain 0 entries too)
class walking_pattern
	: public grandfather_pattern
{
public:
	walking_pattern(const matrix<size_t>& pattern, const size_t n);

	// Returns true if matrix N avoids walking pattern, on which function is called. Using already computed table by updating it from [r,c] and beyond.
	virtual bool avoid(const size_t r, const size_t c, const matrix<size_t>& N);
	virtual bool revert(const size_t r, const size_t c, const matrix<size_t>& N) { return avoid(r, c, N); }
private:
	matrix<std::pair<size_t, size_t> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical) and value of v_i.
	std::vector<size_t> direction_, value_;
};

#endif