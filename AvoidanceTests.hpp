#ifndef AvoidanceTest_hpp_
#define AvoidanceTest_hpp_

#include "Matrix.hpp"

class grandfather_pattern
{
public:
	virtual bool avoid(const size_t r, const size_t c, const matrix<int>& N) { return 0; }
};

// General matrix pattern
class general_pattern
	: public grandfather_pattern
{
public:
	explicit general_pattern(const matrix<int>& p) : pattern_(p) {}
	explicit general_pattern(matrix<int>&& p) : pattern_(std::move(p)) {}

	virtual bool avoid(const size_t r, const size_t c, const matrix<int>& N); // not implemented yet - hopefully I will find more efficient way than brute force
private:
	matrix<int> pattern_;
};

// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all 1 entries (may contain 0 entries too)
class walking_pattern
	: public grandfather_pattern
{
public:
	walking_pattern(const matrix<int>& pattern, const size_t n);

	// Returns true if matrix N avoids walking pattern, on which function is called. Using already computed table by updating it from [r,c] and beyond.
	virtual bool avoid(const size_t r, const size_t c, const matrix<int>& N);
private:
	matrix<std::pair<size_t, size_t> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical) and value of v_i.
	std::vector<int> direction_, value_;
};

#endif