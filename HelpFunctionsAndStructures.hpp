#ifndef HelpFunctions_hpp_
#define HelpFunctions_hpp_

#include "Matrix.hpp"

#include <set>
#include <unordered_set>

enum Type { GENERAL, WALKING };

// Enum for the line ordering functions
enum Order { DESC, SUM, MAX, AUTO, CUSTOM };

enum Map { RECURSION, COMPROMISE, NORECURSION };

enum Map_container { VECTOR, SET, HASH };

// hash function for a vector of size_t
class size_t_vector_hasher {
public:
	size_t operator()(const std::vector<size_t>& vec) const {
		size_t seed = 0;
		for (auto& i : vec) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

struct my_exception : std::exception
{
	my_exception(const char* message) : message_(message) {}

	const char* what() const
	{
		return message_;
	}
private:
	const char* message_;
};

template<class T>
class Container
{
public:
	typedef T::iterator	iterator;
	typedef T::const_iterator const_iterator;

	/// <summary>
	/// Initializes the container by erasing all item inside and filling it with mapping of size 0.
	/// Takes constant time and space.
	/// </summary>
	void init();

	/// <summary>
	/// Inserts a new mapping into a container only if it is not in the container yet.
	/// Takes time depending on the container T.
	/// </summary>
	void insertWithoutDuplicates(std::vector<size_t> mapping);

	iterator begin()		{ return T.begin(); }
	const_iterator cbegin()	{ return T.cbegin(); }
	iterator end()			{ return T.end(); }
	const_iterator cend()	{ return T.cend(); }
private:
	T container_;
};

template<>
inline void Container<std::vector<std::vector<size_t> > >::init()
{
	container_.clear();
	container_.push_back(std::vector<size_t>(0));
}

template<>
inline void Container<std::set<std::vector<size_t> > >::init()
{
	container_.clear();
	container_.insert(std::vector<size_t>(0));
}

template<>
inline void Container<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> >::init()
{
	container_.clear();
	container_.insert(std::vector<size_t>(0));
}

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