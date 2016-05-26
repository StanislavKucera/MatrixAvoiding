#ifdef _MSC_VER 
// MS Visual Studio gives warnings when using fopen.
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef API_cpp_
#define API_cpp_

#include "API.hpp"
#include "GeneralPatternFunctions.hpp"

int main(int argc, char **argv)
{
	int N(100),
		iter(1000),
		threads_count(1),
		random_seed(-1);
	int hist_from(0),
		hist_to(1000),
		hist_freq(100);							// matrix statistics settings
	std::vector<Console_output> console_outputs;
	std::vector<std::string> output_files,
		hist_files,
		max_ones_files,
		csv_files,
		perf_files;
	std::string	init_matrix("zero"),
		config_file("config.txt");
	Parallel_mode parallel_mode(SERIAL);

	if (argc > 1)
		config_file = argv[1];

	{
		FILE* test = fopen(config_file.c_str(), "r");

		if (!test)
		{
			std::cerr << "Cannot open configuration file \"" << config_file << "\". Please check that all directories are created and accessible.\n";
			getchar();
			return 1;
		}
		else
			fclose(test);
	}

	std::ifstream config(config_file);

	std::vector<Pattern_info> pattern_info = parse_config(config, N, iter, random_seed, hist_from, hist_to, hist_freq, console_outputs,
		output_files, hist_files, max_ones_files, csv_files, perf_files, init_matrix, threads_count, parallel_mode);

	config.close();

	std::chrono::system_clock::time_point start, end;
	Patterns patterns;
	Matrix<bool> result;

	if (init_matrix == "zero")
		result = Matrix<bool>(N, N);
	else
	{
		try
		{
			result = Matrix<bool>(N, N, init_matrix);
		}
		catch (...)
		{
			std::cerr << "Cannot read the initial matrix from \"" << init_matrix << "\"." << std::endl;
		}
	}

	set_patterns(patterns, pattern_info, result, init_matrix != "zero", N, threads_count);

	Matrix_Statistics matrix_stats(hist_from, (hist_to == -1) ? iter : hist_to, N, (hist_files.empty()) ? 0 : hist_freq, !max_ones_files.empty());
	Performance_Statistics perf_stats(5, iter);

	//////////////////////////////////////////////////////
	start = std::chrono::system_clock::now();

	if (parallel_mode == MCMC)
		parallelMCMCgenerator(iter, patterns, result, perf_stats, matrix_stats, threads_count - 1, random_seed);
	else if (parallel_mode == MCMC2)
		parallelMCMCgenerator2(iter, patterns, result, perf_stats, matrix_stats, threads_count, random_seed);
	else //if (parallel_mode == MAP)
		MCMCgenerator(iter, patterns, result, perf_stats, matrix_stats, (parallel_mode == SERIAL ? 1 : threads_count), random_seed);

	end = std::chrono::system_clock::now();
	//////////////////////////////////////////////////////

	// if output file is specified
	for (const auto& output_file : output_files)
	{
		if (output_file == "console")
			continue;
		else if (output_file.substr(output_file.size() - 4) == ".bmp")
		{
			BMP matrix;
			matrix.SetSize(N, N);
			matrix.SetBitDepth(1);
			CreateGrayscaleColorTable(matrix);

			for (int i = 0; i != N; ++i)
				for (int j = 0; j != N; ++j)
				{
					matrix(i, j)->Red = (ebmpBYTE)((1 - result.at(i, j)) * 255);
					matrix(i, j)->Green = (ebmpBYTE)((1 - result.at(i, j)) * 255);
					matrix(i, j)->Blue = (ebmpBYTE)((1 - result.at(i, j)) * 255);
				}

			matrix.WriteToFile(output_file.c_str());
		}
		else
		{
			std::ofstream outFile(output_file);
			outFile << result.Print();
			outFile.close();
		}
	}

	// if performance stats file is specified
	for (const auto& perf_file : perf_files)
	{
		if (perf_file == "console")
			continue;

		std::ofstream opFile(perf_file);
		opFile << "Total running time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.\n\n";
		perf_stats.print_data(opFile);
		opFile.close();
	}

	// if performance stats csv file is specified
	for (const auto& csv_file : csv_files)
	{
		if (csv_file == "console")
			continue;

		std::ofstream opcFile(csv_file);
		opcFile << "Total running time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.\n\n";
		perf_stats.print_csv(opcFile);
		opcFile.close();
	}

	// if histogram file is specified
	for (const auto& hist_file : hist_files)
	{
		if (hist_file == "console")
			continue;
		else if (hist_file.substr(hist_file.size() - 4) == ".bmp")
			matrix_stats.print_bmp_histogram(hist_file.c_str());
		else
			matrix_stats.print_text_histogram(hist_file.c_str());
	}

	// if max ones file is specified
	for (const auto& max_ones_file : max_ones_files)
	{
		if (max_ones_file == "console")
			continue;
		else if (max_ones_file.substr(max_ones_file.size() - 4) == ".bmp")
			matrix_stats.print_bmp_max_ones(max_ones_file.c_str());
		else
			matrix_stats.print_text_max_ones(max_ones_file.c_str());
	}

	// console output
	for (const auto& console_output : console_outputs)
	{
		switch (console_output)
		{
		case MATRIX:
			std::cout << result.Print();
			break;
		case TIME:
			std::cout << "Running time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0 << " sec.";
			break;
		case PATTERNS:
			std::cout << "Avoiding patterns:\n\n";

			for (const auto& pat : pattern_info)
			{
				int row, col, val;
				std::ifstream iFile(pat.pattern_file);
				iFile >> row >> col;

				for (int i = 0; i < row; ++i)
				{
					for (int j = 0; j < col; ++j)
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

			break;
		case PERF:
			std::cout << "Performance statistics:\n";
			perf_stats.print_data(std::cout);
			break;
		case CSV:
			std::cout << "Performance csv statistics:\n";
			perf_stats.print_csv(std::cout);
			break;
		case HIST:
			std::cout << "\nHistogram:\n";
			matrix_stats.print_histogram(std::cout);
			break;
		case MAX_ONES:
			std::cout << "\nMatrix with the maximum number of one-entries:\n";
			matrix_stats.print_max_ones(std::cout);
			break;
		default:
			break;
		}

		std::cout << std::endl << std::endl;
	}

	getchar();
	return 0;
}

#endif