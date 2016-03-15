#ifdef _MSC_VER 
// MS Visual Studio gives warnings when using fopen.
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef API_cpp_
#define API_cpp_

#include "API.hpp"
#include "GeneralPatternFunctions.hpp"

int main()
{
	size_t N(100),
		iter(1000),
		threads_count(1);
	size_t hist_from(0),
		hist_to(1000),
		hist_freq(100);							// matrix statistics settings
	std::string bmp_file("notset"),
		hist_file("notset"),
		max_ones_file("notset"),
		csv_file("notset"),
		perf_file("notset"),
		init_matrix("zero");
	bool console_time(false),
		console_pattern(false),
		console_matrix(false),
		console_perf(false),
		console_csv(false),
		console_hist(false),
		console_max_ones(false);
	Parallel_mode parallel_mode(SERIAL);

	std::ifstream config("config.txt");

	std::vector<Pattern_info> pattern_info = parse_config(config, N, iter, hist_from, hist_to, hist_freq, bmp_file, hist_file, max_ones_file, csv_file, perf_file,
		init_matrix, console_time, console_pattern, console_matrix, console_perf, console_csv, console_hist, console_max_ones, threads_count, parallel_mode);

	std::chrono::system_clock::time_point start, end;
	Patterns patterns;
	Matrix<size_t> result;

	if (init_matrix == "zero")
		result = Matrix<size_t>(N, N);
	else
	{
		try
		{
			result = Matrix<size_t>(N, N, init_matrix);
		}
		catch (...)
		{
			std::cerr << "Cannot read the initial matrix from \"" << init_matrix << "\"." << std::endl;
		}
	}

	set_patterns(patterns, pattern_info, result, init_matrix != "zero", N, threads_count);

	Matrix_Statistics matrix_stats(hist_from, hist_to, N, hist_freq);
	Performance_Statistics perf_stats(5, iter);

	//////////////////////////////////////////////////////
	start = std::chrono::system_clock::now();

	if (parallel_mode == MCMC)
		parallelMCMCgenerator(iter, patterns, result, perf_stats, matrix_stats, threads_count - 1);
	else if (parallel_mode == MCMC2)
		parallelMCMCgenerator2(iter, patterns, result, perf_stats, matrix_stats, threads_count);
	else //if (parallel_mode == MAP)
		MCMCgenerator(iter, patterns, result, perf_stats, matrix_stats, (parallel_mode == SERIAL ? 1 : threads_count));

	end = std::chrono::system_clock::now();
	//////////////////////////////////////////////////////

	// if output file is specified
	if (bmp_file != "notset") {
		BMP matrix;
		matrix.SetSize(N, N);
		matrix.SetBitDepth(1);
		CreateGrayscaleColorTable(matrix);

		for (size_t i = 0; i < N; ++i)
			for (size_t j = 0; j < N; ++j)
			{
				matrix(i, j)->Red = (ebmpBYTE)((1 - result.at(i, j)) * 255);
				matrix(i, j)->Green = (ebmpBYTE)((1 - result.at(i, j)) * 255);
				matrix(i, j)->Blue = (ebmpBYTE)((1 - result.at(i, j)) * 255);
			}

		matrix.WriteToFile(bmp_file.c_str());
	}

	// if performance stats file is specified
	if (perf_file != "notset") {
		std::ofstream opFile(perf_file);
		opFile << "Total running time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.\n\n";
		perf_stats.print_data(opFile);
		opFile.close();
	}

	// if performance stats csv file is specified
	if (csv_file != "notset") {
		std::ofstream opcFile(csv_file);
		opcFile << "Total running time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.\n\n";
		perf_stats.print_csv(opcFile);
		opcFile.close();
	}

	// if histogram file is specified
	if (hist_file != "notset")
		matrix_stats.print_histogram(hist_file.c_str());

	// if max ones file is specified
	if (max_ones_file != "notset")
		matrix_stats.print_max_ones(max_ones_file.c_str());

	if (console_matrix)
		std::cout << result.Print();
	if (console_pattern)
	{
		std::cout << "\nAvoiding patterns:\n\n";

		for (const auto& pat : pattern_info)
		{
			size_t row, col, val;
			std::ifstream iFile(pat.pattern_file);
			iFile >> row >> col;

			for (size_t i = 0; i < row; ++i)
			{
				for (size_t j = 0; j < col; ++j)
				{
					iFile >> val;

					if (j != 0)
						std::cout << " ";

					std::cout << val;
				}

				std::cout << std::endl;
			}

			std::cout << std::endl;

			iFile.close();
		}
	}
	if (console_time)
		std::cout << "\nRunning time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.\n";
	if (console_perf)
	{
		std::cout << "\nPerformance statistics:\n";
		perf_stats.print_data(std::cout);
	}
	if (console_csv)
	{
		std::cout << "\nPerformance csv statistics:\n";
		perf_stats.print_csv(std::cout);
	}
	if (console_hist)
	{
		std::cout << "\nHistogram:\n";
		matrix_stats.print_histogram(std::cout);
	}
	if (console_max_ones)
	{
		std::cout << "\nMatrix with the maximum number of one-entries:\n";
		matrix_stats.print_max_ones(std::cout);
	}

	getchar();
	return 0;
}

#endif
