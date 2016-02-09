#ifndef HelpFunctions_hpp_
#define HelpFunctions_hpp_

#include "Matrix.hpp"

#include <set>
#include <unordered_set>

enum Type { GENERAL, WALKING, SLOW };

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

	const char* what() throw() { return message_; }
private:
	const char* message_;
};

template<typename T>
class Container
{
public:
	typedef typename T::iterator		iterator;
	typedef typename T::const_iterator	const_iterator;

	/// <summary>
	/// Initializes the container by erasing all item inside and filling it with mapping of size 0.
	/// Takes constant time and space.
	/// </summary>
	void init();

	/// <summary>
	/// Inserts a new mapping into a container only if it is not in the container yet.
	/// Takes time depending on the container T.
	/// </summary>
	void insert_without_duplicates(std::vector<size_t>&& mapping);

	void clear()					{ container_.clear(); }
	size_t size() const				{ return container_.size(); }

	iterator begin()				{ return container_.begin(); }
	const_iterator cbegin() const	{ return container_.cbegin(); }
	iterator end()					{ return container_.end(); }
	const_iterator cend() const		{ return container_.cend(); }
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

template<>
inline void Container<std::vector<std::vector<size_t> > >::insert_without_duplicates(std::vector<size_t>&& mapping)
{
	// go through all already found mappings in (i+1)-th step and check if extended is not already in there
	for (auto& mapping2 : container_)
	{
		// extended has already been added (atleast its different class) - I won't add it for the second time
		if (mapping == mapping2)
			return;
	}

	// if extended is not yet an element, add it to the tree
	container_.push_back(mapping);
}

template<>
inline void Container<std::set<std::vector<size_t> > >::insert_without_duplicates(std::vector<size_t>&& mapping)
{
	// duplicates are dealt with automagically
	container_.insert(mapping);
}

template<>
inline void Container<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> >::insert_without_duplicates(std::vector<size_t>&& mapping)
{
	// duplicates are dealt with automagically
	container_.insert(mapping);
}

struct Counter
{
	Counter() : tries(0), maps(0), uniques(0) {}

	size_t tries,	// number of attempts to map a line
		maps,		// number of successful attempts to map a line
		uniques;	// number of unique mappings (two mapping are the same if the important mapped lines are the same)
};

/// <summary>
/// Computes the number of one-entries of a binary number.
/// Takes constant time and space.
/// </summary>
/// <param name="n">A binary number for which number of bits is computed.</param>
inline size_t bit_count(const size_t n)	// I have used a function from the internet: -http://blogs.msdn.com/b/jeuge/archive/2005/06/08/hakmem-bit-count.aspx
{
	const size_t uCount = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
	return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

#endif