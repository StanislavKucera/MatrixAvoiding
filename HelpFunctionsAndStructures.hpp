#ifndef HelpFunctions_hpp_
#define HelpFunctions_hpp_

#include "Matrix.hpp"

#include <set>
#include <unordered_set>
#include <mutex>
#include <string>
#include <sstream>
//#include <concurrent_unordered_set.h>

enum Type { GENERAL, WALKING, SLOW };

// Enum for the line ordering functions
enum Order { DESC, SUM, MAX, AUTO, CUSTOM };

struct Map 
{ 
	Map() : enough_entries(true), recursion(true), orthogonal_bounds(true) {}
	Map(const int map)
	{
		switch (map)
		{
		case 0:
			enough_entries = recursion = orthogonal_bounds = false;
			break;
		case 1:
			enough_entries = true;
			recursion = orthogonal_bounds = false;
			break;
		case 2:
			orthogonal_bounds = true;
			enough_entries = recursion = false;
			break;
		case 3:
			enough_entries = orthogonal_bounds = true;
			recursion = false;
			break;
		case 4:
			enough_entries = recursion = true;
			orthogonal_bounds = false;
			break;
		case 5:
			enough_entries = recursion = orthogonal_bounds = true;
			break;
		default:
			enough_entries = recursion = orthogonal_bounds = true;
			break;
		}
	}

	bool enough_entries,
		recursion,
		orthogonal_bounds;
};

enum Map_container { VECTOR, SET, HASH };

enum Parallel_mode { SERIAL, MCMC, MCMC2, MAP };

// hash function for a vector of size_t
class int_vector_hasher {
public:
	int operator()(const std::vector<int>& vec) const
	{
		int seed = 0;

		for (auto& i : vec)
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);

		return seed;
	}
};

struct my_exception : std::exception
{
	my_exception(const char* message) : message_(message) {}

	const char* what() const throw() { return message_; }
private:
	const char* message_;
};

template<typename T>
class Container
{
public:
	typedef typename T::iterator		iterator;
	typedef typename T::const_iterator	const_iterator;

	Container() : container_(), write_() {}
	Container(const Container<T>& cont) : container_(cont.container_), write_() {}

	/// <summary>
	/// Initializes the container by erasing all item inside and filling it with mapping of size 0.
	/// Takes constant time and space.
	/// </summary>
	void init();

	/// <summary>
	/// Inserts a new mapping into a container only if it is not in the container yet.
	/// Takes time depending on the container T.
	/// </summary>
	void insert_without_duplicates(std::vector<int>&& mapping);
	void insert_without_duplicates(const std::vector<int>& mapping);
	void insert_without_duplicates(const Container<T>& mappings);

	void parallel_insert_without_duplicates(std::vector<int>&& mapping);

	void clear()					{ container_.clear(); }
	int size() const				{ return (int)container_.size(); }
	bool empty() const				{ return container_.empty(); }

	iterator begin()				{ return container_.begin(); }
	const_iterator begin() const	{ return container_.cbegin(); }
	const_iterator cbegin() const	{ return container_.cbegin(); }
	iterator end()					{ return container_.end(); }
	const_iterator end() const		{ return container_.cend(); }
	const_iterator cend() const		{ return container_.cend(); }
private:
	T container_;
	std::mutex write_;
};

template<>
inline void Container<std::vector<std::vector<int> > >::init()
{
	container_.clear();
	container_.emplace_back(0);
}

template<>
inline void Container<std::set<std::vector<int> > >::init()
{
	container_.clear();
	container_.emplace(0);
}

template<>
inline void Container<std::unordered_set<std::vector<int>, int_vector_hasher> >::init()
{
	container_.clear();
	container_.emplace(0);
}

template<>
inline void Container<std::vector<std::vector<int> > >::insert_without_duplicates(std::vector<int>&& mapping)
{
	// go through all already found mappings in (i+1)-th step and check if extended is not already in there
	for (auto& mapping2 : container_)
	{
		// extended has already been added (atleast its different class) - I won't add it for the second time
		if (mapping == mapping2)
			return;
	}

	// if extended is not yet an element, add it to the tree
	container_.emplace_back(std::move(mapping));
}

template<>
inline void Container<std::vector<std::vector<int> > >::insert_without_duplicates(const std::vector<int>& mapping)
{
	// go through all already found mappings in (i+1)-th step and check if extended is not already in there
	for (auto& mapping2 : container_)
	{
		// extended has already been added (atleast its different class) - I won't add it for the second time
		if (mapping == mapping2)
			return;
	}

	// if extended is not yet an element, add it to the tree
	container_.emplace_back(mapping);
}

template<>
inline void Container<std::vector<std::vector<int> > >::insert_without_duplicates(const Container<std::vector<std::vector<int> > >& mappings)
{
	for (const std::vector<int>& mapping : mappings)
	{
		// go through all already found mappings in (i+1)-th step and check if extended is not already in there
		for (const auto& mapping2 : container_)
		{
			// extended has already been added (atleast its different class) - I won't add it for the second time
			if (mapping == mapping2)
				return;
		}

		// if extended is not yet an element, add it to the tree
		container_.emplace_back(mapping);
	}
}

template<>
inline void Container<std::set<std::vector<int> > >::insert_without_duplicates(std::vector<int>&& mapping)
{
	// duplicates are dealt with automagically
	container_.emplace(mapping);
}

template<>
inline void Container<std::set<std::vector<int> > >::insert_without_duplicates(const std::vector<int>& mapping)
{
	// duplicates are dealt with automagically
	container_.emplace(mapping);
}

template<>
inline void Container<std::set<std::vector<int> > >::insert_without_duplicates(const Container<std::set<std::vector<int> > >& mappings)
{
	// duplicates are dealt with automagically
	container_.insert(mappings.cbegin(), mappings.cend());
}

template<>
inline void Container<std::unordered_set<std::vector<int>, int_vector_hasher> >::insert_without_duplicates(std::vector<int>&& mapping)
{
	// duplicates are dealt with automagically
	container_.emplace(mapping);
}

template<>
inline void Container<std::unordered_set<std::vector<int>, int_vector_hasher> >::insert_without_duplicates(const std::vector<int>& mapping)
{
	// duplicates are dealt with automagically
	container_.emplace(mapping);
}

template<>
inline void Container<std::unordered_set<std::vector<int>, int_vector_hasher> >::insert_without_duplicates(const Container<std::unordered_set<std::vector<int>, int_vector_hasher> >& mappings)
{
	// duplicates are dealt with automagically
	container_.insert(mappings.cbegin(), mappings.cend());
}

template<>
inline void Container<std::vector<std::vector<int> > >::parallel_insert_without_duplicates(std::vector<int>&& mapping)
{
	// I will play with this later
	write_.lock();
	// go through all already found mappings in (i+1)-th step and check if extended is not already in there
	for (auto& mapping2 : container_)
	{
		// extended has already been added (atleast its different class) - I won't add it for the second time
		if (mapping == mapping2)
			return;
	}

	// if extended is not yet an element, add it to the tree
	container_.emplace_back(std::move(mapping));

	write_.unlock();
}

template<>
inline void Container<std::set<std::vector<int> > >::parallel_insert_without_duplicates(std::vector<int>&& mapping)
{
	std::unique_lock<std::mutex> lck(write_);
	container_.emplace(mapping);
}

template<>
inline void Container<std::unordered_set<std::vector<int>, int_vector_hasher> >::parallel_insert_without_duplicates(std::vector<int>&& mapping)
{
	std::unique_lock<std::mutex> lck(write_);
	container_.emplace(mapping);
}

struct Counter
{
	Counter() : tries(0), maps(0), uniques(0) {}

	int tries,		// number of attempts to map a line
		maps,		// number of successful attempts to map a line
		uniques;	// number of unique mappings (two mapping are the same if the important mapped lines are the same)
};

struct Job
{
	Job() : r(-1), c(-1), avoid(true) {}
	Job(const int row, const int col, bool a) : r(row), c(col), avoid(a) {}
	Job(const Job& j) : r(j.r), c(j.c), avoid(j.avoid) {}
	
	int r, c;
	bool avoid;
};

inline bool operator<(const Job&, const Job&) { return false; }

struct Task
{
	Task() : job(), id(-1), next_id(-1), returned(false), synced(false) {}
	Task(const Job& j, const int i, const int ni, const bool r, const bool s) : job(j), id(i), next_id(ni), returned(r), synced(s) {}

	Job job;
	int id,
		next_id;
	bool returned,
		synced;
};

/// <summary>
/// Computes the number of one-entries of a binary number.
/// Takes constant time and space.
/// </summary>
/// <param name="n">A binary number for which number of bits is computed.</param>
inline int bit_count(const int n)	// I have used a function from the internet: -http://blogs.msdn.com/b/jeuge/archive/2005/06/08/hakmem-bit-count.aspx
{
	const int uCount = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
	return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

inline int my_stoi(const std::string& s) {
	std::istringstream str(s);
	int ret;
	str >> ret;
	return ret;
}

#endif