#ifdef _MSC_VER 
// MS Visual Studio gives warnings when using fopen.
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef API_cpp_
#define API_cpp_

#include "GeneralPatternFunctions.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"

#include <set>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <string>																							//								^
#include <sstream>																							//								|
#include <assert.h>																							//     <---------> horizontal	| vertical	
#include <exception>																						//								|
#include <stdexcept>																						//								v
#include <time.h>


inline Type getType(const std::string& type)
{
	if (type == "general" || type == "GENERAL" || type == "G" || type == "g")
		return GENERAL;
	else if (type == "walking" || type == "WALKING" || type == "W" || type == "w")
		return WALKING;
	else if (type == "slow" || type == "SLOW" || type == "S" || type == "s")
		return SLOW;
	else {
		assert(!"Pattern type not supported.");
		throw my_exception("Pattern type not supported. Choose WALKING or GENERAL.");
		return GENERAL;
	}
}

inline Map getMap(const std::string& map)
{
	if (map == "recursion" || map == "RECURSION" || map == "R" || map == "r" || map == "REC" || map == "Rec" || map == "rec")
		return RECURSION;
	else if (map == "compromise" || map == "COMPROMISE" || map == "C" || map == "c" || map == "COM" || map == "Com" || map == "com")
		return COMPROMISE;
	else if (map == "norecursion" || map == "NORECURSION" || map == "N" || map == "n" || map == "NOR" || map == "Nor" || map == "nor")
		return NORECURSION;
	else {
		assert(!"Mapping approach not supported.");
		throw my_exception("Mapping approach not supported. Choose RECURSION, NORECURSION or COMPROMISE.");
		return RECURSION;
	}
}

inline Map_container getContainer(const std::string& container)
{
	if (container == "vector" || container == "VECTOR" || container == "V" || container == "v" || container == "std::vector")
		return VECTOR;
	else if (container == "set" || container == "SET" || container == "S" || container == "s" || container == "std::set")
		return SET;
	else if (container == "hash" || container == "HASH" || container == "H" || container == "h" || container == "std::unordered_set")
		return HASH;
	else {
		assert(!"Map container not supported.");
		throw my_exception("Map container not supported. Choose VECTOR or SET.");
		return VECTOR;
	}
}

inline Order getOrder(const std::string& order)
{
	if (order == "desc" || order == "DESC" || order == "D" || order == "d")
		return DESC;
	else if (order == "max" || order == "MAX" || order == "M" || order == "m")
		return MAX;
	else if (order == "sum" || order == "SUM" || order == "S" || order == "s")
		return SUM;
	else if (order == "custom" || order == "CUSTOM" || order == "C" || order == "c")
		return CUSTOM;
	else if (order == "auto" || order == "AUTO" || order == "A" || order == "a" || "automatic")
		return AUTO;
	else {
		assert(!"Line ordering not supported.");
		throw my_exception("Line ordering not supported. Choose DESC, MAX, SUM, AUTO or CUSTOM.");
		return DESC;
	}
}

inline bool getBool(const std::string& write)
{
	if (write == "yes" || write == "YES" || write == "Y" || write == "y")
		return true;
	else if (write == "no" || write == "NO" || write == "N" || write == "n")
		return false;
	else {
		assert(!"Not sure what to write into the console.");
		throw my_exception("Not sure what to write into the console. Choose yes or no in config.txt.");
		return false;
	}
}

inline std::vector<size_t> readOrder(const std::string& custom_order)
{
	std::vector<size_t> ret;
	std::ifstream iFile(custom_order);
	size_t val;

	while (!iFile.eof())
	{
		iFile >> val;
		ret.push_back(val);
	}

	return ret;
}

int main()
{
	size_t N, iter;
	size_t from, to, freq;							// matrix statistics settings
	clock_t t;
	std::string junk, param;
	bool initialized = false;

	std::ifstream config("config.txt");

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Program input:
	for (size_t i = 0; i < 5; ++i) config >> junk;	// Size of the generated matrix:
	config >> N;									// #size#

	for (size_t i = 0; i < 6; ++i) config >> junk;	// Number of iterations of the generator:
	config >> iter;									// #iterations#

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Pattern file:
	config >> param;								// #pattern#
	Matrix<size_t> pattern(param);

	for (size_t i = 0; i < 4; ++i) config >> junk;	// Type of the pattern:
	config >> param;								// #type#
	Type type = getType(param);

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Initial big matrix:
	config >> param;								// #init#
	Matrix<size_t> result;
	if (param == "zero")
		result = Matrix<size_t>(N, N);				// zero matrix NxN
	else {
		result = Matrix<size_t>(N, N, param);
		initialized = true;
	}

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Map function approach:
	config >> param;								// #map#
	Map map = getMap(param);

	for (size_t i = 0; i < 4; ++i) config >> junk;	// Container for storing mappings:
	config >> param;								// #container#
	Map_container container = getContainer(param);

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Line ordering:
	config >> param;								// #order#
	Order order = getOrder(param);
	std::vector<size_t> custom_order;
	if (order == CUSTOM) {
		config >> param;							// file with custom order
		custom_order = readOrder(param);
	}


	for (size_t i = 0; i < 2; ++i) config >> junk;	// Program output:
	for (size_t i = 0; i < 3; ++i) config >> junk;	// Output matrix file:
	config >> param;								// #output#
	std::string output;
	if (param != "no" && param != "n")				// no output file
	{
		output = param;
		FILE* test = fopen(output.c_str(), "w");
		if (!test)
			std::cout << "Cannot open file \"" << output << "\". Please check all directories are created and accessible.\n";
		else
			fclose(test);
	}

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Write resulting matrix into console:
	config >> param;								// #cmatrix#
	bool console_matrix = getBool(param);

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Write the pattern into console:
	config >> param;								// #cpattern#
	bool console_pattern = getBool(param);

	for (size_t i = 0; i < 8; ++i) config >> junk;	// Write total time of the run into console:
	config >> param;								// #ctime#
	bool console_time = getBool(param);

	for (size_t i = 0; i < 4; ++i) config >> junk;	// Output performance statistics file:
	config >> param;								// #perf#
	std::string perf_file;
	if (param != "no" && param != "n")				// no output file
	{
		perf_file = param;
		FILE* test = fopen(perf_file.c_str(), "w");
		if (!test)
			std::cout << "Cannot open file \"" << perf_file << "\". Please check all directories are created and accessible.\n";
		else
			fclose(test);
	}

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Write performance statistics into console:
	config >> param;								// #cperf#
	bool performance = getBool(param);

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Output performance statistics csv file:
	config >> param;								// #csv#
	std::string csv_file;
	if (param != "no" && param != "n")				// no output file
	{
		csv_file = param;
		FILE* test = fopen(csv_file.c_str(), "w");
		if (!test)
			std::cout << "Cannot open file \"" << csv_file << "\". Please check all directories are created and accessible.\n";
		else
			fclose(test);
	}

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Matrix statistics:
	for (size_t i = 0; i < 3; ++i) config >> junk;	// Histogram refresh frequency:
	config >> freq;									// #freq#
	if (freq == 0)
		freq = (size_t)-1;

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Histogram refresh from:
	config >> from;									// #freq#

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Histogram refresh to:
	config >> to;									// #to#
	if (to == 0)
		to = iter;

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Histogram output file:
	config >> param;								// #histogram#
	std::string hist;
	if (param != "no" && param != "n")				// no output file
	{
		hist = param;
		FILE* test = fopen(hist.c_str(), "w");
		if (!test)
			std::cout << "Cannot open file \"" << hist << "\". Please check all directories are created and accessible.\n";
		else
			fclose(test);
	}

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Max ones matrix output file:
	config >> param;								// #max_ones#
	std::string max_ones;
	if (param != "no" && param != "n")				// no output file
	{
		max_ones = param;
		FILE* test = fopen(max_ones.c_str(), "w");
		if (!test)
			std::cout << "Cannot open file \"" << max_ones << "\". Please check all directories are created and accessible.\n";
		else
			fclose(test);
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//result = Matrix<size_t>("matrix/try.txt");
	//Walking_pattern wp(pattern, N);
	//std::vector<Counter> sizes;
	//bool ahoj = wp.avoid(result, sizes, 0, 0);


	Matrix_Statistics matrix_stats(from, to, N, freq);
	Performance_Statistics perf_stats(5, iter);
	std::vector<Counter> sizes;

	Patterns patterns;

	if (type == WALKING)
		patterns.add(new Walking_pattern(pattern, N));
	else if (type == SLOW)
		patterns.add(new Slow_pattern(pattern));
	else if (type == GENERAL && container == VECTOR)
	{
		if (order == AUTO)
		{
			if (initialized) 
			{
				General_pattern<std::vector<std::vector<size_t> > > gpSUM(pattern, SUM, map);
				t = clock();

				if (!gpSUM.avoid(result, sizes)) 
				{
					assert(!"Initial big matrix does not avoid the pattern");
					throw my_exception("Initial big matrix does not avoid the pattern");
				}

				t = clock() - t;
				auto sum_time = t;

				General_pattern<std::vector<std::vector<size_t> > > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(result, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::vector<std::vector<size_t> > > gpDESC(pattern, DESC, map);
				t = clock();
				gpDESC.avoid(result, sizes);
				t = clock() - t;
				auto desc_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
			else 
			{
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

				General_pattern<std::vector<std::vector<size_t> > > gpSUM(pattern, SUM, map);
				t = clock();
				gpSUM.avoid(rand, sizes);
				t = clock() - t;
				auto sum_time = t;

				General_pattern<std::vector<std::vector<size_t> > > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(rand, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::vector<std::vector<size_t> > > gpDESC(pattern, DESC, map);
				t = clock();
				gpDESC.avoid(rand, sizes);
				t = clock() - t;
				auto desc_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
		}

		patterns.add(new General_pattern<std::vector<std::vector<size_t> > >(pattern, order, map, std::move(custom_order)));
	}
	else if (type == GENERAL && container == SET)
	{
		if (order == AUTO)
		{
			if (initialized)
			{
				General_pattern<std::set<std::vector<size_t> > > gpDESC(pattern, DESC, map);
				t = clock();

				if (!gpDESC.avoid(result, sizes))
				{
					assert(!"Initial big matrix does not avoid the pattern");
					throw my_exception("Initial big matrix does not avoid the pattern");
				}

				t = clock() - t;
				auto desc_time = t;

				General_pattern<std::set<std::vector<size_t> > > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(result, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::set<std::vector<size_t> > > gpSUM(pattern, SUM, map);
				t = clock();
				gpSUM.avoid(result, sizes);
				t = clock() - t;
				auto sum_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
			else
			{
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

				General_pattern<std::set<std::vector<size_t> > > gpSUM(pattern, SUM, map);
				t = clock();
				gpSUM.avoid(rand, sizes);
				t = clock() - t;
				auto sum_time = t;

				General_pattern<std::set<std::vector<size_t> > > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(rand, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::set<std::vector<size_t> > > gpDESC(pattern, DESC, map);
				t = clock();
				gpDESC.avoid(rand, sizes);
				t = clock() - t;
				auto desc_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
		}

		patterns.add(new General_pattern<std::set<std::vector<size_t> > >(pattern, order, map, std::move(custom_order)));
	}
	else if (type == GENERAL && container == HASH)
	{
		if (order == AUTO)
		{
			if (initialized)
			{
				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpSUM(pattern, SUM, map);
				t = clock();

				if (!gpSUM.avoid(result, sizes))
				{
					assert(!"Initial big matrix does not avoid the pattern");
					throw my_exception("Initial big matrix does not avoid the pattern");
				}

				t = clock() - t;
				auto sum_time = t;

				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(result, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpDESC(pattern, DESC, map);
				t = clock();
				gpDESC.avoid(result, sizes);
				t = clock() - t;
				auto desc_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
			else
			{
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpSUM(pattern, SUM, map);
				t = clock();
				gpSUM.avoid(rand, sizes);
				t = clock() - t;
				auto sum_time = t;

				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpMAX(pattern, MAX, map);
				t = clock();
				gpMAX.avoid(rand, sizes);
				t = clock() - t;
				auto max_time = t;

				General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpDESC(pattern, DESC, map);
				t = clock();
				gpDESC.avoid(rand, sizes);
				t = clock() - t;
				auto desc_time = t;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
		}

		patterns.add(new General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> >(pattern, order, map, std::move(custom_order)));
	}
	else {
		assert(!"Something bad have happened.");
		throw my_exception("This shouldn't have happened. Please send a message to the developer team with your config file.");
	}

	//////////////////////////////////////////////////////
	t = clock();
	//parallelMCMCgenerator(iter, patterns, result, perf_stats, matrix_stats, 2);
	MCMCgenerator(iter, patterns, result, perf_stats, matrix_stats);
	t = clock() - t;
	//////////////////////////////////////////////////////

	// if output file is specified
	if (output != "") {
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

		matrix.WriteToFile(output.c_str());
	}

	// if performance stats file is specified
	if (perf_file != "") {
		std::ofstream opFile(perf_file);
		opFile << "Total running time: " << (double)t / CLOCKS_PER_SEC << " sec.\n\n";
		perf_stats.print_data(opFile);
		opFile.close();
	}

	// if performance stats csv file is specified
	if (csv_file != "") {
		std::ofstream opcFile(csv_file);
		opcFile << "Total running time: " << (double)t / CLOCKS_PER_SEC << " sec.\n\n";
		perf_stats.print_csv(opcFile);
		opcFile.close();
	}

	// if histogram file is specified
	if (hist != "")
		matrix_stats.print_histogram(hist.c_str());

	// if max ones file is specified
	if (max_ones != "")
		matrix_stats.print_max_ones(max_ones.c_str());

	if (console_matrix)
		std::cout << result.Print();
	if (console_pattern)
		std::cout << "\nAvoiding pattern:\n\n" << pattern.Print();
	if (console_time)
		std::cout << "\nRunning time: " << (double)t / CLOCKS_PER_SEC << " sec.\n";
	if (performance)
	{
		std::cout << "\nPerformance statistics:\n";
		perf_stats.print_data(std::cout);
	}

	getchar();
	return 0;
}

#endif
