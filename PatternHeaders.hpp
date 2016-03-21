#ifndef PatternHeaders_hpp_
#define PatternHeaders_hpp_

#include "HelpFunctionsAndStructures.hpp"

#include <atomic>
#include <queue>
#include <condition_variable>
#include <thread>

/// For the purposes of complexity, let \theta(k) be the number of lines (rows and columns) of the pattern
/// and \theta(n) be the number of lines of the big matrix.

class Pattern
{
public:
	virtual bool avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const std::atomic_bool& force_end = false) = 0;
	virtual bool revert(const Matrix<bool>& big_matrix, const int r, const int c) = 0;
	virtual bool parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const int threads_count, const std::atomic_bool& force_end = false) = 0;
	virtual bool parallel_revert(const Matrix<bool>& big_matrix, const int r, const int c, const int threads_count) = 0;
	virtual std::vector<int> get_order() const = 0;
	virtual Pattern* get_new_instance() const = 0;
	virtual void construct_threads(const Matrix<bool>& big_matrix) = 0;
	virtual void destruct_threads() = 0;
};

class Slow_pattern
	: public Pattern
{
public:
	Slow_pattern(const Matrix<bool>& pattern) : one_entries_(0), rows_(pattern.getRow()), cols_(pattern.getCol())
	{
		for (int i = 0; i < pattern.getRow(); ++i)
			for (int j = 0; j < pattern.getCol(); ++j)
				if (pattern.at(i, j) == 1)
					one_entries_.push_back(std::make_pair(i, j));
	}

	bool avoid(const Matrix<bool>& big_matrix, const int /* r */, const int /* c */, std::vector<Counter>& /* sizes */, const std::atomic_bool& force_end = 0)
	{
		done_ = false;
		// goes through all subsets of rows and columns of the right cardinality and tests whether the pattern can be mapped to that subset
		test_all_subsets(0ll, 0ll, rows_, cols_, big_matrix.getRow(), big_matrix.getCol(), big_matrix, force_end);

		if (done_)
			return false;

		return true;
	}
	bool revert(const Matrix<bool>& /* big_matrix */, const int /* r */, const int /* c */) { return true; }
	bool parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const int threads_count, const std::atomic_bool& force_end = 0);
	bool parallel_revert(const Matrix<bool>& /* big_matrix */, const int /* r */, const int /* c */, const int /* threads_count */) { return true; }
	std::vector<int> get_order() const { return std::vector<int>(); }
	Pattern* get_new_instance() const { return new Slow_pattern(*this); }
	void construct_threads(const Matrix<bool>& /* big_matrix */) {}
	void destruct_threads() {}
private:
	std::vector<std::pair<int, int> > one_entries_;			// list of all one entries of the pattern
	const int rows_, cols_;									// size of the pattern
	bool done_;												// indicator whether the avoidance testing has failed (the matrix does not avoid the pattern)

	void test_all_subsets(int v_map, int h_map, int v_ones, int h_ones, int v_vals, int h_vals, const Matrix<bool>& big_matrix, const std::atomic_bool& force_end);
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
	General_pattern(const Matrix<bool>& pattern, const int threads_count, const Order order = DESC, const Map map_approach = SUPERACTIVE,  std::vector<int>&& custom_order = std::vector<int>());
	General_pattern(const General_pattern<T>& copy) : row_(copy.row_), col_(copy.col_), lines_(copy.lines_), order_(copy.order_), what_to_remember_(copy.what_to_remember_),
		parallel_bound_indices_(copy.parallel_bound_indices_), extending_order_(copy.extending_order_), map_index_(copy.map_index_), building_tree_(2), steps_(copy.steps_),
		empty_lines_(copy.empty_lines_), map_approach_(copy.map_approach_) {}

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
	bool avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const std::atomic_bool& force_end = false);
	bool revert(const Matrix<bool>& /* big_matrix */, const int /* r */, const int /* c */) { return true; }
	bool parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const int threads_count, const std::atomic_bool& force_end = false);
	bool parallel_revert(const Matrix<bool>& /* big_matrix */, const int /* r */, const int /* c */, const int /* threads_count */) { return true; }
	std::vector<int> get_order() const { return order_; }
	Pattern* get_new_instance() const { return new General_pattern(*this); }
	void construct_threads(const Matrix<bool>& big_matrix)
	{
		for (int index = 0; index != (int)threads_.size(); ++index)
			threads_[index] = std::thread(&General_pattern::worker, this, index, std::ref(big_matrix));
	}
	void destruct_threads()
	{
		end_ = true;

		for (int index = 0; index < (int)threads_.size(); ++index)
		{
			{
				std::unique_lock<std::mutex> lck(mtxs_[index]);
				cvs_[index].notify_one();
			}

			threads_[index].join();
		}
	}
private:
	const int	row_,									// number of rows of the pattern
				col_;									// number of columns of the pattern
	std::vector<int>	lines_,							// binary number for each line of a pattern having one at i-th position if the pattern has one-entry there
														// lines_[i] = (1011)_2 ... i-th line of the pattern has one-enty at 0th, 1st and 3rd position
						order_,							// order of lines in which I am going to be mapping them
														// order_[i] = j ... in i-th step I'm going to map j-th line if the pattern
						what_to_remember_;				// for each adding line I know which of them I still need to remember for next mapping
														// what_to_remember[i] = (001010)_2 ... in i-th step I remember where I mapped the 1st and 3rd line of the pattern
	std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > > parallel_bound_indices_;
		// vector through levels - vector through lines - pair of pairs - pair of lower and upper bounds
		// parallel_bound_indices_[i][j] = ((bot, top),(i_bot, i_top)) ... in i-th step, j-th line of the pattern is bounded by i_bot line of the pattern from the bottom
		//	and by i_top line of the pattern from the top; bot and top are indices to mapping structure - mapping[bot] = b ... bot-th line is mapped to b line
	std::vector<std::vector<int> > extending_order_;	// vector through levels - vector of indices of the mapping which are needed for the extended one
														// extending_order_[i][j] = k ... in i-th step, j-th linewill be emplaced at k-th position of the mapping
	std::vector<std::vector<int> > map_index_;			// vector through levels - vector of indices of lines in the mapping
														// map_index_[i][j] = k ... in i-th step, j-th line is on the k-th position in the mapping
	std::vector<Container<T> > building_tree_;			// container for found mappings at each level
	int	steps_,											// number of steps I'm going to do = number of lines I need to map (excluding empty lines)
		empty_lines_,									// binary number of lines with no one-entries
		level_;											// the level of the computation - index of the line being mapped
	const Map map_approach_;							// choosen way of mapping algorithm - use recursion for nonmapped lines or not

	/// <summary>
	/// For given line of the pattern computes lines of the big matrix, which bound its mapping.
	/// The bounds are stored to from and to variables meaning that "line" can be mapped to lines [from, to).
	/// Takes constant time since it knows where to look, because indices to the vectors are precalculated.
	/// </summary>
	/// <param name="line">Index of the line of the pattern for which bounds are calculated.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="rows">Number of rows of the big matrix.</param>
	/// <param name="columns">Number of columns of the big matrix.</param>
	/// <param name="from">Index of the line of the big matrix which bounds "line" from the bottom.</param>
	/// <param name="to">Index of the line of the big matrix which bounds "line" from the top.</param>
	/// <param name="r">Row of the entry that was changed in the last iteration if I know it.</param>
	/// <param name="c">Column of the entry that was changed in the last iteration if I know it.</param>
	void find_parallel_bounds(const int line, const std::vector<int>& mapping, const int rows, const int columns, int& from, int& to, const int r = -1, const int c = -1) const;
	
	/// <summary>
	/// For given line of the pattern computes lines of the big matrix, which bound its mapping.
	/// The bounds are stored to from and to variables meaning that "line" can be mapped to lines [from, to).
	/// Takes constant time since it knows where to look, because indices to the vectors are precalculated.
	/// </summary>
	/// <param name="line">Index of the line of the pattern for which bounds are calculated.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="rows">Number of rows of the big matrix.</param>
	/// <param name="columns">Number of columns of the big matrix.</param>
	/// <param name="from">Index of the line of the big matrix which bounds "line" from the bottom.</param>
	/// <param name="to">Index of the line of the big matrix which bounds "line" from the top.</param>
	/// <param name="r">Row of the entry that was changed in the last iteration if I know it.</param>
	/// <param name="c">Column of the entry that was changed in the last iteration if I know it.</param>
	bool check_orthogonal_bounds(const int line, const int big_line, const std::vector<int>& mapping, const int orthogonal_line, const int big_orthognal_line, const Matrix<bool>& big_matrix) const;

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
	int count_what_to_remember(const int current) const;
	
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
	void find_bound_indices(const int line, const int level);

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
	/// <param name="big_line">Index of the line of the big matrix which I am trying to map the line to.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	/// <param name="big_matrix">Reference to the big matrix for which I test pattern avoiding.</param>
	bool map(const bool backtrack, const int line, const int big_line, const std::vector<int>& mapping, const Matrix<bool>& big_matrix);
	
	/// <summary>
	/// Extends previous mapping after deciding to which big line the line should be mapped.
	/// It is precomputed which elements of previous mapping are needed for the new one and where to put the new element.
	/// Takes k time, for extending.
	/// </summary>
	/// <param name="return">The extended mapping</param>
	/// <param name="level">The level I am at - how many lines I have mapped already.</param>
	/// <param name="big_line">Index of the line of the big matrix which I mapped the line to.</param>
	/// <param name="mapping">The mapping I am extending.</param>
	std::vector<int> extend(const int big_line, const std::vector<int>& mapping) const;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	// parallel avoid stuff //
	//////////////////////////

	// thread pool
	std::vector<std::thread> threads_;
	// queue of found mappings for each thread - stores all mapping found by its worker and those are then taken by the main thread and added to building_tree_ without duplicates
	// accessed by exactly two threads - push by its worker; front, pop by the main thread - mutex needed
	std::vector<std::queue<std::vector<int> > > qs_;
	// mutex for each queue
	std::vector<std::mutex> mutexes_;
	// condition variable for each thread - if there is nothing to compute for a worker it waits for the condition variable
	// accessed by exactly two threads - wait by its worker; notify_one by the main thread
	std::vector<std::condition_variable> cvs_;
	// mutex for each condition variable
	std::vector<std::mutex> mtxs_;
	// indicator whether a thread waits for its condition variable or is computing something
	// accessed by exactly two threads - write by its worker; read by the main thread - mutex not needed because of atomicity
	std::vector<std::atomic_bool> sleeps_;
	// condition variable for the main thread - if workers calculate new mappings and none of them is done yet, the main threads waits
	// accessed by all the threads - notify_one by all workers; wait by the main thread
	std::condition_variable cv_;
	// mutex for the condition variable
	std::mutex mtx_;
	// pointer to a mapping that is being extended
	// accessed by all threads - read by all workers; write by the main thread - mutex not needed because it is writen only when workers sleep
	typename T::const_pointer mapping_ptr_ = nullptr;
	// index of a line which currently mapped line is being mapped to
	// accessed by all threads - read and write by all workers; write by the main thread - mutex not needed because of atomicity
	std::atomic_int big_line_;
	// index of a line behind the last one that is possible to map the currently mapped line to. If big_line = big_line_to tests are done
	// accessed by all threads - read by all workers; write by the main thread - mutex not needed because of atomicity (this really doesn't need to be atomic)
	std::atomic_int big_line_to_;
	// are all possible mappings of the current line tested? It is TRUE when there is no more work for a thread (big_line = big_line_to)
	// accessed by all threads - read and write by all workers; read and write by the main thread - mutex not needed because of atomicity
	std::atomic_bool done_;
	// indicator that the MCMCgenerator ends - it need to destroy all workers
	std::atomic_bool end_;
	// indicator that there is some work for the main thread and it wasn't just woke up randomly
	std::atomic_bool something_is_mapped_or_done_;

	void worker(const int index, const Matrix<bool>& big_matrix);
};

/// A matrix pattern in which exists a walk from left-upper corner to right-bottom corner, which contains all one-entries.
/// It may contain zero-entries as well, but no one-entry can be not included in the walk.
class Walking_pattern
	: public Pattern
{
public:
	Walking_pattern(const Matrix<bool>& pattern, const int n);
	
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
	bool avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const std::atomic_bool& force_end = false);	
	// reverts changes in max_walk_part matrix after an unsuccessful change of the big matrix
	bool revert(const Matrix<bool>& /* big_matrix */, const int r, const int c) { changes_.emplace_back(std::make_pair(r, c)); return true; }
	bool parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<Counter>& sizes, const int threads_count, const std::atomic_bool& force_end = false);
	bool parallel_revert(const Matrix<bool>& big_matrix, const int r, const int c, const int threads_count)
	{ 
		std::vector<Counter> sizes;
		return parallel_avoid(big_matrix, r, c, sizes, threads_count);
	}
	std::vector<int> get_order() const { return std::vector<int>(); }
	Pattern* get_new_instance() const { return new Walking_pattern(*this); }
	virtual void construct_threads(const Matrix<bool>& /* big_matrix */) {}
	void destruct_threads() {}
private:
	Matrix<std::pair<int, int> > max_walk_part_;	// table of calculated [c_v,c_h] for all elements
	
	// indexed by index of v_i, the element of the walk, gives the direction of the next element (0 for vertical)
	std::vector<bool> direction_;
	// indexed by index of v_i, the element of the walk, gives the value of v_i.
	std::vector<int> value_;
	std::vector<std::pair<int, int> > changes_;
	bool top_left;
	int size_;
};

class Patterns
{
public:
	Patterns() : patterns_(), big_matrix_(), changed_(-1) {}
	Patterns(const Patterns& copy) : big_matrix_(copy.big_matrix_), changed_(-1)
	{
		for (auto& pattern : copy.patterns_)
			patterns_.push_back(pattern->get_new_instance());
	}

	bool avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<std::vector<Counter> >& sizes);
	bool revert(const Matrix<bool>& big_matrix, const int r, const int c);
	bool parallel_avoid(const Matrix<bool>& big_matrix, const int r, const int c, std::vector<std::vector<Counter> >& sizes, const int threads_count);
	bool parallel_revert(const Matrix<bool>& big_matrix, const int r, const int c, const int threads_count);
	std::vector<std::vector<int> > get_order() const;

	bool avoid(const int r, const int c, std::vector<std::vector<Counter> >& sizes, const std::atomic_bool& forced_end);
	bool revert(const int r, const int c);// , const Matrix<bool>& mat);
	bool check_matrix(const Matrix<bool>& mat);

	void add(Pattern* pattern) { patterns_.push_back(pattern); }
	void set_matrix(const Matrix<bool>& big_matrix) { big_matrix_ = big_matrix; }
	void construct_threads(const Matrix<bool>& big_matrix);
	void destruct_threads();
private:
	std::vector<Pattern*> patterns_;
	Matrix<bool> big_matrix_;
	int changed_;
	std::vector<std::pair<int, int> > changes;
};

#endif