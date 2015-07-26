#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <exception>
#include <time.h>

Type getType(const std::string& type)
{
	if (type == "general" || type == "GENERAL" || type == "G" || type == "g")
		return GENERAL;
	else if (type == "walking" || type == "WALKING" || type == "W" || type == "w")
		return WALKING;
	else {
		assert(!"Type not supported");
		throw new std::exception("Type not supported. Choose WALKING or GENERAL");
		return GENERAL;
	}
}

Map getMap(const std::string& map)
{
	if (map == "recursion" || map == "RECURSION" || map == "R" || map == "r" || map == "REC" || map == "Rec" || map =="rec")
		return RECURSION;
	else if (map == "compromise" || map == "COMPROMISE" || map == "C" || map == "c" || map == "COM" || map == "Com" || map == "com")
		return COMPROMISE;
	else if (map == "norecursion" || map == "NORECURSION" || map == "N" || map == "n" || map == "NOR" || map == "Nor" || map == "nor")
		return NORECURSION;
	else {
		assert(!"Mapping approach not supported");
		throw new std::exception("Mapping approach not supported. Choose RECURSION, NORECURSION or COMPROMISE");
		return RECURSION;
	}
}

Order getOrder(const std::string& order)
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
		assert(!"Line ordering not supported");
		throw new std::exception("Line ordering not supported. Choose DESC, MAX, SUM, AUTO or CUSTOM");
		return DESC;
	}
}

bool getBool(const std::string& write)
{
	if (write == "yes" || write == "YES" || write == "Y" || write == "y")
		return true;
	else if (write == "no" || write == "NO" || write == "N" || write == "n")
		return false;
	else {
		assert(!"Not sure what to write into the console.");
		throw new std::exception("Not sure what to write into the console. Choose yes or no in config.txt");
		return false;
	}
}

std::vector<size_t> readOrder(const std::string& custom_order)
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

int main(int argc, char* argv[])
{
	size_t N, iter;
	clock_t t;
	std::string junk, param;
	bool initialized = false;

#pragma region parsing config.txt
	std::ifstream config("config.txt");

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Program input:
	for (size_t i = 0; i < 5; ++i) config >> junk;	// Size of the generated matrix:
	config >> N;									// #size#

	for (size_t i = 0; i < 6; ++i) config >> junk;	// Number of iterations of the generator:
	config >> iter;									// #iterations#

	for (size_t i = 0; i < 2; ++i) config >> junk;	// Pattern file:
	config >> param;								// #pattern#
	matrix<size_t> pattern(param);

	for (size_t i = 0; i < 4; ++i) config >> junk;	// Type of the pattern:
	config >> param;								// #type#
	Type type = getType(param);

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Initial big matrix:
	config >> param;								// #init#
	matrix<size_t> result;
	if (param == "zero")
		result = matrix<size_t>(N, N);				// zero matrix NxN
	else {
		result = matrix<size_t>(N, N, param);
		initialized = true;
	}

	for (size_t i = 0; i < 3; ++i) config >> junk;	// Map function approach:
	config >> param;								// #map#
	Map map = getMap(param);

	for (size_t i = 0; i < 4; ++i) config >> junk;	// Container for storing mappings:
	config >> param;								// #container#
	std::string container = param;

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
		output = param;

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Write resulting matrix into console:
	config >> param;								// #cmatrix#
	bool console_matrix = getBool(param);

	for (size_t i = 0; i < 5; ++i) config >> junk;	// Write the pattern into console:
	config >> param;								// #cpattern#
	bool console_pattern = getBool(param);

	for (size_t i = 0; i < 8; ++i) config >> junk;	// Write total time of the run into console:
	config >> param;								// #ctime#
	bool console_time = getBool(param);
#pragma endregion


	if (type == WALKING) {
		walking_pattern wp(pattern, N);
		t = clock();
		MCMCgenerator(iter, wp, result);
		t = clock() - t;
	}
	else if (type == GENERAL) {
		if (order == AUTO) {
			if (initialized) {
				general_pattern gpDESC(pattern, DESC, map);
				size_t num_iter = 0;
				t = clock();

				while (clock() - t <= 10) {
					if (!gpDESC.avoid(result)) {
						assert(!"Initial big matrix does not avoid the pattern");
						throw new std::exception("Initial big matrix does not avoid the pattern");
					}
					++num_iter;
				}

				t = clock();
				for (size_t i = 0; i < num_iter; ++i)
					gpDESC.avoid(result);
				t = clock() - t;
				auto desc_time = (float)t / CLOCKS_PER_SEC;

				general_pattern gpMAX(pattern, MAX, map);
				t = clock();
				for (size_t i = 0; i < num_iter; ++i)
					gpDESC.avoid(result);
				t = clock() - t;
				auto max_time = (float)t / CLOCKS_PER_SEC;

				general_pattern gpSUM(pattern, SUM, map);
				t = clock();
				for (size_t i = 0; i < num_iter; ++i)
					gpDESC.avoid(result);
				t = clock() - t;
				auto sum_time = (float)t / CLOCKS_PER_SEC;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
			else {
				matrix<size_t> rand = matrix<size_t>::random_bin_matrix(N, N);

				general_pattern gpDESC(pattern, DESC, map);
				size_t num_iter = 0;
				t = clock();

				while (clock() - t <= 10) {
					gpDESC.avoid(result);
					++num_iter;
				}

				t = clock();
				gpDESC.avoid(rand);
				t = clock() - t;
				auto desc_time = (float)t / CLOCKS_PER_SEC;

				general_pattern gpMAX(pattern, MAX, map);
				t = clock();
				gpDESC.avoid(rand);
				t = clock() - t;
				auto max_time = (float)t / CLOCKS_PER_SEC;

				general_pattern gpSUM(pattern, SUM, map);
				t = clock();
				gpDESC.avoid(rand);
				t = clock() - t;
				auto sum_time = (float)t / CLOCKS_PER_SEC;

				if (desc_time <= max_time && desc_time <= sum_time)
					order = DESC;
				else if (max_time <= sum_time)
					order = MAX;
				else
					order = SUM;
			}
		}

		general_pattern gp(pattern, order, map, std::move(custom_order));
		t = clock();
		MCMCgenerator(iter, gp, result);
		t = clock() - t;
	}

	// if output file is specified
	if (output != "") {
		std::ofstream oFile(output);
		// print the resulting matrix into it
		oFile << result.Print();
		oFile.close();
	}

	if (console_matrix)
		std::cout << result.Print();
	if (console_pattern)
		std::cout << "\nAvoiding pattern:\n\n" << pattern.Print();
	if (console_time)
		std::cout << "\nRunning time: " << (float)t / CLOCKS_PER_SEC << " sec.\n";

	getchar();
	return 0;
}

#endif