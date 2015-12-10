#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include "HelpFunctionsAndStructures.hpp"

#include <ostream>
#include <time.h>

// matrix statistics

// performance statistics - is size_t big enough?
class Performance_Statistics
{
public:
	Performance_Statistics() : success_sizes(0), fail_sizes(0), success_counter(0), success_levels(0), success_time(0), fail_counter(0), fail_levels(0), fail_time(0) {}

	void addData(size_t iter, bool success, size_t time, const std::vector<Counter>& sizes)
	{
		// maybe use iter for something? Data may be very different for iter in [0,1000] and iter in [1M,100M]
		iter = iter;

		if (sizes.size() > success_sizes.size())
		{
			success_sizes.resize(sizes.size());
			fail_sizes.resize(sizes.size());
		}

		if (success)
		{
			++success_counter;
			success_time += time;
			success_levels += sizes.size();

			for (size_t i = 0; i < sizes.size(); ++i)
			{
				success_sizes[i].tries += sizes[i].tries;
				success_sizes[i].maps += sizes[i].maps;
				success_sizes[i].uniques += sizes[i].uniques;
			}
		}
		else
		{
			++fail_counter;
			fail_time += time;
			fail_levels += sizes.size();

			for (size_t i = 0; i < sizes.size(); ++i)
			{
				fail_sizes[i].tries += sizes[i].tries;
				fail_sizes[i].maps += sizes[i].maps;
				fail_sizes[i].uniques += sizes[i].uniques;
			}
		}
	}
	void printData(std::ostream& output)
	{
		const size_t iter = success_counter + fail_counter;
		const double s_counter = (double)success_counter;
		const double f_counter = (double)fail_counter;

		output << "Successful calls of avoid (matrix avoids the pattern):\n";
		output << "Count: " << success_counter << " - " << s_counter / iter * 100 << "%\n";
		output << "Average time per call: " << success_time / s_counter / CLOCKS_PER_SEC << " sec\n";
		output << "Average number of lines mapped: " << success_levels / s_counter << "\n";
		output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";
		for (size_t i = 0; i < success_sizes.size(); ++i)
		{
			output << i + 1 << ": "
				<< success_sizes[i].tries / s_counter << " - "
				<< success_sizes[i].maps / s_counter << " - "
				<< success_sizes[i].uniques / s_counter << "\n";
		}

		output << "\nUnsuccessful calls of avoid (matrix doesn't avoid the pattern):\n";
		output << "Count: " << fail_counter << " - " << f_counter / iter * 100 << "%\n";
		output << "Average time per call: " << fail_time / f_counter / CLOCKS_PER_SEC << " sec\n";
		output << "Average number of lines mapped: " << fail_levels / f_counter << "\n";
		output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";
		for (size_t i = 0; i < fail_sizes.size(); ++i)
		{
			output << i + 1 << ": "
				<< fail_sizes[i].tries / f_counter << " - "
				<< fail_sizes[i].maps / f_counter << " - "
				<< fail_sizes[i].uniques / f_counter << "\n";
		}
	}
private:
	std::vector<Counter> success_sizes,		// for each mapped line, there is a sum of all mappings found during the whole MCMC generation process while mapping the line
		fail_sizes;
	size_t success_counter,					// count of the successful calls of avoid(), successful mean the matrix did avoid the pattern
		success_levels,						// sum of the numbers of lines mapped until the end of avoid
		success_time,						// total time (in processor cycles) spent in successful calls of avoid
		fail_counter,
		fail_levels,
		fail_time;
};

#endif