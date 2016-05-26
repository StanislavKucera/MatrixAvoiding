#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include "HelpFunctionsAndStructures.hpp"
#include "EasyBMP/EasyBMP.h"

#include <ostream>

// matrix statistics
class Matrix_Statistics
{
public:
	Matrix_Statistics(const int from, const int to, const int n, const int freq, const bool ones) : all_entries(n, n), max_all_entries(n, n),
		ones_count(0), max_ones_count(0), hist_count(0), iter_from(from), iter_to(to), frequency(freq), count_max_ones(ones) {}

	void add_data(const int iter, const int ones, const Matrix<bool>& big_matrix)
	{
		if (count_max_ones && ones > max_ones_count)
		{
			max_ones_count = ones;
			max_all_entries = big_matrix;
		}

		if (iter < iter_from || iter > iter_to || frequency == 0 || iter % frequency == 0)
			return;

		for (int i = 0; i < big_matrix.getRow(); ++i)
			for (int j = 0; j < big_matrix.getCol(); ++j)
				if (big_matrix.at(i, j))
					++all_entries.at(i, j);

		ones_count += ones;
		++hist_count;
	}
	void print_bmp_histogram(const char* output) const
	{
		BMP hist;
		hist.SetSize(all_entries.getRow(), all_entries.getCol());
		hist.SetBitDepth(8);
		CreateGrayscaleColorTable(hist);

		for (int i = 0; i < all_entries.getRow(); ++i)
			for (int j = 0; j < all_entries.getCol(); ++j)
			{
				const ebmpBYTE color = (ebmpBYTE)((1 - (all_entries.at(i, j) / (double)hist_count)) * 255);
				hist(i, j)->Red = color;
				hist(i, j)->Green = color;
				hist(i, j)->Blue = color;
			}

		hist.WriteToFile(output);
	}
	void print_text_histogram(const char* output) const
	{
		std::ofstream oFile(output);
		oFile << all_entries.Print() << std::endl << hist_count;
		oFile.close();
	}
	void print_histogram(std::ostream& output) const
	{
		for (int i = 0; i < all_entries.getRow(); i++)
		{
			for (int j = 0; j < all_entries.getCol(); j++)
			{
				if (j != 0)
					output << " ";

				output << all_entries.at(i, j) / (double)hist_count;
			}

			output << "\n";
		}
	}
	void print_bmp_max_ones(const char* output) const
	{
		BMP matrix;
		matrix.SetSize(max_all_entries.getRow(), max_all_entries.getCol());
		matrix.SetBitDepth(1);
		CreateGrayscaleColorTable(matrix);

		for (int i = 0; i < max_all_entries.getRow(); ++i)
			for (int j = 0; j < max_all_entries.getCol(); ++j)
			{
				const ebmpBYTE color = (ebmpBYTE)((1 - max_all_entries.at(i, j)) * 255);
				matrix(i, j)->Red = color;
				matrix(i, j)->Green = color;
				matrix(i, j)->Blue = color;
			}

		matrix.WriteToFile(output);
	}
	void print_text_max_ones(const char* output) const
	{
		std::ofstream oFile(output);
		oFile << max_all_entries.Print();
		oFile.close();
	}
	void print_max_ones(std::ostream& output) const { output << max_all_entries.Print(); }
private:
	Matrix<long long> all_entries;
	Matrix<bool> max_all_entries;
	long long ones_count,
		max_ones_count,
		hist_count;
	const int iter_from,
		iter_to,
		frequency;
	const bool count_max_ones;
};

// performance statistics - is size_t big enough?
class Performance_Statistics
{
public:
	Performance_Statistics(const int p, const int i) : success_sizes(0), fail_sizes(0), max_success_sizes(0), max_fail_sizes(0),
		orders(0), success_time(0), fail_time(0), success_counter(0), success_levels(0), fail_counter(0), fail_levels(0), mod(i / p), iter(i) {}

	void add_data(const int iter, const bool success, const double time, const std::vector<std::vector<Counter> >& sizes)
	{
		const size_t index = iter / mod;

        if (success_sizes.size() <= index)
        {
            success_sizes.emplace_back();
			fail_sizes.emplace_back();
			max_success_sizes.emplace_back();
			max_fail_sizes.emplace_back();
			success_counter.emplace_back(0);
			success_levels.emplace_back(0);
			success_time.emplace_back(0);
			fail_counter.emplace_back(0);
			fail_levels.emplace_back(0);
			fail_time.emplace_back(0);
        }

		if (success_sizes[index].size() < sizes.size())
		{
			success_sizes[index].resize(sizes.size());
			fail_sizes[index].resize(sizes.size());
			max_success_sizes[index].resize(sizes.size());
			max_fail_sizes[index].resize(sizes.size());
		}

		if (success)
		{
			++success_counter[index];
			success_time[index] += time;
			
			if (sizes.empty())
				return;

			success_levels[index] += sizes[0].size();
		}
		else
		{
			++fail_counter[index];
			fail_time[index] += time;

			if (sizes.empty())
				return;

			fail_levels[index] += sizes[0].size();
		}

		for (size_t ind = 0; ind < sizes.size(); ++ind)
		{
			if (sizes[ind].size() > success_sizes[index][ind].size())
			{
				success_sizes[index][ind].resize(sizes[ind].size());
				fail_sizes[index][ind].resize(sizes[ind].size());
				max_success_sizes[index][ind].resize(sizes[ind].size());
				max_fail_sizes[index][ind].resize(sizes[ind].size());
			}

			if (success)
			{
				for (size_t i = 0; i < sizes[ind].size(); ++i)
				{
					success_sizes[index][ind][i].tries += sizes[ind][i].tries;
					if (sizes[ind][i].tries > max_success_sizes[index][ind][i].tries)
						max_success_sizes[index][ind][i].tries = sizes[ind][i].tries;

					success_sizes[index][ind][i].maps += sizes[ind][i].maps;
					if (sizes[ind][i].maps > max_success_sizes[index][ind][i].maps)
						max_success_sizes[index][ind][i].maps = sizes[ind][i].maps;

					success_sizes[index][ind][i].uniques += sizes[ind][i].uniques;
					if (sizes[ind][i].uniques > max_success_sizes[index][ind][i].uniques)
						max_success_sizes[index][ind][i].uniques = sizes[ind][i].uniques;
				}
			}
			else
			{
				for (size_t i = 0; i < sizes[ind].size(); ++i)
				{
					fail_sizes[index][ind][i].tries += sizes[ind][i].tries;
					if (sizes[ind][i].tries > max_fail_sizes[index][ind][i].tries)
						max_fail_sizes[index][ind][i].tries = sizes[ind][i].tries;

					fail_sizes[index][ind][i].maps += sizes[ind][i].maps;
					if (sizes[ind][i].maps > max_fail_sizes[index][ind][i].maps)
						max_fail_sizes[index][ind][i].maps = sizes[ind][i].maps;

					fail_sizes[index][ind][i].uniques += sizes[ind][i].uniques;
					if (sizes[ind][i].uniques > max_fail_sizes[index][ind][i].uniques)
						max_fail_sizes[index][ind][i].uniques = sizes[ind][i].uniques;
				}
			}
		}
	}
	void set_order(std::vector<std::vector<int> >&& o, bool parallel = false) { orders = std::move(o); parallel_ = parallel; }
	void print_data(std::ostream& output) const
	{
		output << "Used orders:\n";

		for (auto& order : orders)
		{
			if (!order.empty())
				for (size_t i = 0; i < order.size(); ++i)
					output << " " << order[i];
			else
				output << " No order.";

			output << "\n";
		}

		if (parallel_)
			return;

		output << "\nSuccessful calls of avoid (matrix avoids the pattern):\n";

		for (size_t ind = 0; ind < orders.size(); ++ind)
		{
			output << "Pattern " << ind + 1 << ":\n";

			for (int m = 0; m < iter / mod; ++m)
			{
				const double s_counter = (double)success_counter[m];

				output << "Iteration: " << m << "\n";
				output << "Count: " << success_counter[m] << " - " << s_counter / (success_counter[m] + fail_counter[m]) * 100 << "%\n";

				if (success_counter[m] != 0)
				{
					output << "Average time per call: " << success_time[m] / s_counter / CLOCKS_PER_SEC << " sec\n";

					// not a general pattern
					if (orders[ind].empty())
					{
						output << "\n";
						continue;
					}

					output << "Average number of lines mapped: " << success_levels[m] / s_counter << "\n";
					output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";

					for (size_t i = 0; i < success_sizes[m][ind].size(); ++i)
					{
						output << i + 1 << ": "
							<< success_sizes[m][ind][i].tries / s_counter << " (" << max_success_sizes[m][ind][i].tries << ") - "
							<< success_sizes[m][ind][i].maps / s_counter << " (" << max_success_sizes[m][ind][i].maps << ") - "
							<< success_sizes[m][ind][i].uniques / s_counter << " (" << max_success_sizes[m][ind][i].uniques << ")\n";
					}
				}

				output << "\n";
			}
		}

		output << "Unsuccessful calls of avoid (matrix doesn't avoid the pattern):\n";

		for (size_t ind = 0; ind < orders.size(); ++ind)
		{
			output << "Pattern " << ind + 1 << ":\n";

			for (int m = 0; m < iter / mod; ++m)
			{
				const double f_counter = (double)fail_counter[m];

				output << "Iteration: " << m << "\n";
				output << "Count: " << fail_counter[m] << " - " << f_counter / (success_counter[m] + fail_counter[m]) * 100 << "%\n";

				if (fail_counter[m] != 0)
				{
					output << "Average time per call: " << fail_time[m] / f_counter / CLOCKS_PER_SEC << " sec\n";

					// not a general pattern
					if (orders[ind].empty())
					{
						output << "\n";
						continue;
					}

					output << "Average number of lines mapped: " << fail_levels[m] / f_counter << "\n";
					output << "Average number of mapping attempts - successful mappings - unique mappings for each mapped line:\n";

					for (size_t i = 0; i < fail_sizes[m][ind].size(); ++i)
					{
						output << i + 1 << ": "
							<< fail_sizes[m][ind][i].tries / f_counter << " (" << max_fail_sizes[m][ind][i].tries << ") - "
							<< fail_sizes[m][ind][i].maps / f_counter << " (" << max_fail_sizes[m][ind][i].maps << ") - "
							<< fail_sizes[m][ind][i].uniques / f_counter << " (" << max_fail_sizes[m][ind][i].uniques << ")\n";
					}
				}

				output << "\n";
			}
		}
	}
	void print_csv(std::ostream& output) const
	{
		output << "Used orders:\n";

		for (auto& order : orders)
		{
			if (!order.empty())
				for (size_t i = 0; i < order.size(); ++i)
					output << ";" << order[i];
			else
				output << ";No order.";

			output << "\n";
		}

		if (parallel_)
			return;

		output << "Success\n";
		output << "from;to;time;count;ratio;average call time;average lines mapped;line number;average map calls;average maps found;average unique maps;max map calls;max maps found;max unique maps";
		output << "\n";

		for (size_t ind = 0; ind < orders.size(); ++ind)
		{
			output << "Pattern " << ind + 1 << ";";

			for (int m = 0; m < iter / mod; ++m)
			{
				const double s_counter = (double)success_counter[m];

				output << m * mod + 1 << ";" << (m + 1) * mod << ";";														// from;to;
				output << success_time[m] / CLOCKS_PER_SEC << " sec;";														// time;
				output << success_counter[m] << ";" << s_counter / (success_counter[m] + fail_counter[m]) * 100 << " %;";	// count;ratio;

				if (success_counter[m] != 0)
				{
					output << success_time[m] / s_counter / CLOCKS_PER_SEC << " sec;";											// average call time;

					// not a general pattern
					if (orders[ind].empty())
					{
						output << "\n";
						continue;
					}

					output << success_levels[m] / s_counter;																	// average lines mapped

					for (size_t i = 0; i < success_sizes[m][ind].size(); ++i)
					{
						output << ";" << i + 1 << ";"
							<< success_sizes[m][ind][i].tries / s_counter << ";"
							<< success_sizes[m][ind][i].maps / s_counter << ";"
							<< success_sizes[m][ind][i].uniques / s_counter << ";"
							<< max_success_sizes[m][ind][i].tries << ";"
							<< max_success_sizes[m][ind][i].maps << ";"
							<< max_success_sizes[m][ind][i].uniques;
					}
				}

				output << "\n";
			}
		}

		output << "\nFail\n";
		output << "from;to;time;count;ratio;average call time;average lines mapped;line number;average map calls;average maps found;average unique maps;max map calls;max maps found;max unique maps";
		output << "\n";

		for (size_t ind = 0; ind < orders.size(); ++ind)
		{
			output << "Pattern " << ind + 1 << ";";

			for (int m = 0; m < iter / mod; ++m)
			{
				const double f_counter = (double)fail_counter[m];

				output << m * mod + 1 << ";" << (m + 1) * mod << ";";													// from;to;
				output << fail_time[m] / CLOCKS_PER_SEC << " sec;";														// time;
				output << fail_counter[m] << ";" << f_counter / (success_counter[m] + fail_counter[m]) * 100 << " %;";	// count;ratio;

				if (fail_counter[m] != 0)
				{
					output << fail_time[m] / f_counter / CLOCKS_PER_SEC << " sec;";											// average call time;

					// not a general pattern
					if (orders[ind].empty())
					{
						output << "\n";
						continue;
					}

					output << fail_levels[m] / f_counter;																	// average lines mapped

					for (size_t i = 0; i < fail_sizes[m][ind].size(); ++i)
					{
						output << ";" << i + 1 << ";"
							<< fail_sizes[m][ind][i].tries / f_counter << ";"
							<< fail_sizes[m][ind][i].maps / f_counter << ";"
							<< fail_sizes[m][ind][i].uniques / f_counter << ";"
							<< max_fail_sizes[m][ind][i].tries << ";"
							<< max_fail_sizes[m][ind][i].maps << ";"
							<< max_fail_sizes[m][ind][i].uniques;
					}
				}
			}

			output << "\n";
		}
	}
private:
	std::vector<std::vector<std::vector<Counter> > > success_sizes,		// for each mapped line, there is a sum of all mappings found during the whole MCMC generation process while mapping the line
													 fail_sizes;
	std::vector<std::vector<std::vector<Counter> > > max_success_sizes,	// for each mapped line, there is a sum of all mappings found during the whole MCMC generation process while mapping the line
													 max_fail_sizes;
	std::vector<std::vector<int> > orders;								// order of the lines when using general pattern
	std::vector<double> success_time,									// total time (in processor cycles) spent in successful calls of avoid
						fail_time;
	std::vector<int>	success_counter,								// count of the successful calls of avoid(), successful mean the matrix did avoid the pattern
						success_levels,									// sum of the numbers of lines mapped until the end of avoid
						fail_counter,
						fail_levels;
	const int	mod,
				iter;
	bool parallel_;
};

#endif
