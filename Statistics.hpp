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
	Performance_Statistics(const size_t p, const size_t i) : success_sizes(0), fail_sizes(0), max_success_sizes(0), max_fail_sizes(0),
		success_counter(0), success_levels(0), success_time(0), fail_counter(0), fail_levels(0), fail_time(0), mod(i / p), iter(i) {}

	void addData(size_t iter, bool success, size_t time, const std::vector<Counter>& sizes)
	{
		const size_t index = iter / mod;

        if (success_sizes.size() <= index)
        {
            success_sizes.push_back(std::vector<Counter>());
			fail_sizes.push_back(std::vector<Counter>());
			max_success_sizes.push_back(std::vector<Counter>());
			max_fail_sizes.push_back(std::vector<Counter>());
			success_counter.push_back(0);
			success_levels.push_back(0);
			success_time.push_back(0);
			fail_counter.push_back(0);
			fail_levels.push_back(0);
			fail_time.push_back(0);
        }

		if (sizes.size() > success_sizes[index].size())
		{
			size_t size = success_sizes[index].size();

			success_sizes[index].resize(sizes.size());
			fail_sizes[index].resize(sizes.size());
			max_success_sizes[index].resize(sizes.size());
			max_fail_sizes[index].resize(sizes.size());

			for (; size < sizes.size(); ++size)
			{
				max_success_sizes[index][size].tries = 0;
				max_success_sizes[index][size].maps = 0;
				max_success_sizes[index][size].uniques = 0;
				max_fail_sizes[index][size].tries = 0;
				max_fail_sizes[index][size].maps = 0;
				max_fail_sizes[index][size].uniques = 0;
			}
		}

		if (success)
		{
			++success_counter[index];
			success_time[index] += time;
			success_levels[index] += sizes.size();

			for (size_t i = 0; i < sizes.size(); ++i)
			{
				success_sizes[index][i].tries += sizes[i].tries;
				if (sizes[i].tries > max_success_sizes[index][i].tries)
					max_success_sizes[index][i].tries = sizes[i].tries;
				success_sizes[index][i].maps += sizes[i].maps;
				if (sizes[i].maps > max_success_sizes[index][i].maps)
					max_success_sizes[index][i].maps = sizes[i].maps;
				success_sizes[index][i].uniques += sizes[i].uniques;
				if (sizes[i].uniques > max_success_sizes[index][i].uniques)
					max_success_sizes[index][i].uniques = sizes[i].uniques;
			}
		}
		else
		{
			++fail_counter[index];
			fail_time[index] += time;
			fail_levels[index] += sizes.size();

			for (size_t i = 0; i < sizes.size(); ++i)
			{
				fail_sizes[index][i].tries += sizes[i].tries;
				if (sizes[i].tries > max_fail_sizes[index][i].tries)
					max_fail_sizes[index][i].tries = sizes[i].tries;
				fail_sizes[index][i].maps += sizes[i].maps;
				if (sizes[i].maps > max_fail_sizes[index][i].maps)
					max_fail_sizes[index][i].maps = sizes[i].maps;
				fail_sizes[index][i].uniques += sizes[i].uniques;
				if (sizes[i].uniques > max_fail_sizes[index][i].uniques)
					max_fail_sizes[index][i].uniques = sizes[i].uniques;
			}
		}
	}
	void printData(std::ostream& output)
	{
		output << "Successful calls of avoid (matrix avoids the pattern):\n";
		for (size_t m = 0; m < iter / mod; ++m)
		{
			output << "Iteration: " << m << "\n";
			const double s_counter = (double)success_counter[m];

			output << "Count: " << success_counter[m] << " - " << s_counter / (success_counter[m] + fail_counter[m]) * 100 << "%\n";
			if (success_counter[m] != 0)
			{
				output << "Average time per call: " << success_time[m] / s_counter / CLOCKS_PER_SEC << " sec\n";
				output << "Average number of lines mapped: " << success_levels[m] / s_counter << "\n";
				output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";
				for (size_t i = 0; i < success_sizes[m].size(); ++i)
				{
					output << i + 1 << ": "
						<< success_sizes[m][i].tries / s_counter << " (" << max_success_sizes[m][i].tries << ") - "
						<< success_sizes[m][i].maps / s_counter << " (" << max_success_sizes[m][i].maps << ") - "
						<< success_sizes[m][i].uniques / s_counter << " (" << max_success_sizes[m][i].uniques << ")\n";
				}
			}
			output << "\n";
		}

		output << "Unsuccessful calls of avoid (matrix doesn't avoid the pattern):\n";
		for (size_t m = 0; m < iter / mod; ++m)
		{
			output << "Iteration: " << m << "\n";
			const double f_counter = (double)fail_counter[m];

			output << "Count: " << fail_counter[m] << " - " << f_counter / (success_counter[m] + fail_counter[m]) * 100 << "%\n";
			if (fail_counter[m] != 0)
			{
				output << "Average time per call: " << fail_time[m] / f_counter / CLOCKS_PER_SEC << " sec\n";
				output << "Average number of lines mapped: " << fail_levels[m] / f_counter << "\n";
				output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";
				for (size_t i = 0; i < fail_sizes[m].size(); ++i)
				{
					output << i + 1 << ": "
						<< fail_sizes[m][i].tries / f_counter << " (" << max_fail_sizes[m][i].tries << ") - "
						<< fail_sizes[m][i].maps / f_counter << " (" << max_fail_sizes[m][i].maps << ") - "
						<< fail_sizes[m][i].uniques / f_counter << " (" << max_fail_sizes[m][i].uniques << ")\n";
				}
			}
			output << "\n";
        }
	}
	void printCsv(std::ostream& output)
	{
		output << "Success\n";
		output << "from;to;time;count;ratio;average call time;average lines mapped;line number;average map calls;average maps found;average unique maps;max map calls;max maps found;max unique maps";
		output << "\n";

		for (size_t m = 0; m < iter / mod; ++m)
		{
			const double s_counter = (double)success_counter[m];

			output << m * mod + 1 << ";" << (m + 1) * mod << ";";														// from;to;
			output << success_time[m] / CLOCKS_PER_SEC << " sec;";														// time;
			output << success_counter[m] << ";" << s_counter / (success_counter[m] + fail_counter[m]) * 100 << " %;";	// count;ratio;
			if (success_counter[m] != 0)
			{
				output << success_time[m] / s_counter / CLOCKS_PER_SEC << " sec;";											// average call time;
				output << success_levels[m] / s_counter;																	// average lines mapped
				for (size_t i = 0; i < success_sizes[m].size(); ++i)
				{
					output << ";" << i + 1 << ";"
						<< success_sizes[m][i].tries / s_counter << ";"
						<< success_sizes[m][i].maps / s_counter << ";"
						<< success_sizes[m][i].uniques / s_counter << ";"
						<< max_success_sizes[m][i].tries << ";"
						<< max_success_sizes[m][i].maps << ";"
						<< max_success_sizes[m][i].uniques;
				}
			}
			output << "\n";
		}

		output << "\nFail\n";
		output << "from;to;time;count;ratio;average call time;average lines mapped;line number;average map calls;average maps found;average unique maps;max map calls;max maps found;max unique maps";
		output << "\n";

		for (size_t m = 0; m < iter / mod; ++m)
		{
			const double f_counter = (double)fail_counter[m];

			output << m * mod + 1 << ";" << (m + 1) * mod << ";";													// from;to;
			output << fail_time[m] / CLOCKS_PER_SEC << " sec;";														// time;
			output << fail_counter[m] << ";" << f_counter / (success_counter[m] + fail_counter[m]) * 100 << " %;";	// count;ratio;
			if (fail_counter[m] != 0)
			{
				output << fail_time[m] / f_counter / CLOCKS_PER_SEC << " sec;";											// average call time;
				output << fail_levels[m] / f_counter;																	// average lines mapped
				for (size_t i = 0; i < fail_sizes[m].size(); ++i)
				{
					output << ";" << i + 1 << ";"
						<< fail_sizes[m][i].tries / f_counter << ";"
						<< fail_sizes[m][i].maps / f_counter << ";"
						<< fail_sizes[m][i].uniques / f_counter << ";"
						<< max_fail_sizes[m][i].tries << ";"
						<< max_fail_sizes[m][i].maps << ";"
						<< max_fail_sizes[m][i].uniques;
				}
			}
			output << "\n";
		}
	}
private:
	std::vector<std::vector<Counter> > success_sizes,		// for each mapped line, there is a sum of all mappings found during the whole MCMC generation process while mapping the line
		fail_sizes;
	std::vector<std::vector<Counter> > max_success_sizes,		// for each mapped line, there is a sum of all mappings found during the whole MCMC generation process while mapping the line
		max_fail_sizes;
	std::vector<size_t> success_counter,					// count of the successful calls of avoid(), successful mean the matrix did avoid the pattern
		success_levels,						// sum of the numbers of lines mapped until the end of avoid
		success_time,						// total time (in processor cycles) spent in successful calls of avoid
		fail_counter,
		fail_levels,
		fail_time;
	size_t mod,
		iter;
};

#endif
