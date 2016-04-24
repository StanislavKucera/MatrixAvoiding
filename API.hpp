#ifdef _MSC_VER 
// MS Visual Studio gives warnings when using fopen.
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef API_hpp_
#define API_hpp_

#include "HelpFunctionsAndStructures.hpp"
#include "Matrix.hpp"
#include "PatternHeaders.hpp"
#include "Statistics.hpp"
#include "MCMC.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>

struct Pattern_info
{
	Pattern_info() : pattern_file(""), init_matrix_file("default"), custom_order_file(""), type(GENERAL), map(5), map_container(HASH), order(AUTO) {}

	std::string pattern_file, init_matrix_file, custom_order_file;
	Type type;
	Map map;
	Map_container map_container;
	Order order;
};

inline Type get_type(const std::string& type)
{
	if (type == "general" ||  type == "g")
		return GENERAL;
	else if (type == "walking" || type == "w")
		return WALKING;
	else if (type == "slow" || type == "s")
		return SLOW;
	else {
		std::cerr << "Pattern type \"" << type << "\" not supported. Choose WALKING, GENERAL or SLOW." << std::endl;
		return GENERAL;
	}
}

inline Map get_map(const std::string& map)
{
	if (map == "5" || map == "super active" || map == "super_active" || map == "superactive")
		return Map(5);
	else if(map == "4" || map == "active")
		return Map(4);
	else if (map == "3" || map == "semi active" || map == "semi_active" || map == "semiactive")
		return Map(3);
	else if (map == "2" || map == "semi lazy" || map == "semi_lazy" || map == "semilazy")
		return Map(2);
	else if (map == "1" || map == "lazy")
		return Map(1);
	else if (map == "0" || map == "super lazy" || map == "super_lazy" || map == "superlazy")
		return Map(0);
	else {
		std::cerr << "Mapping approach \"" << map << "\" not supported. Choose SUPERACTIVE, ACTIVE, SEMIACTIVE, SEMILAZY, LAZY or SUPERLAZY." << std::endl;
		return Map(5);
	}
}

inline Map_container get_map_container(const std::string& container)
{
	if (container == "vector" || container == "v")
		return VECTOR;
	else if (container == "set" || container == "s")
		return SET;
	else if (container == "hash" || container == "h")
		return HASH;
	else {
		std::cerr << "Map container \"" << container << "\" not supported. Choose VECTOR, SET or HASH." << std::endl;
		return HASH;
	}
}

inline Order get_order(const std::string& order)
{
	if (order == "desc" || order == "d")
		return DESC;
	else if (order == "max" || order == "m")
		return MAX;
	else if (order == "sum" || order == "s")
		return SUM;
	else if (order == "auto" || order == "a" || "automatic")
		return AUTO;
	else
		return CUSTOM;
}

inline Parallel_mode get_parallel_mode(const std::string& mode)
{
	if (mode == "serial" || mode == "s")
		return SERIAL;
	else if (mode == "mcmc")
		return MCMC;
	else if (mode == "mcmc2")
		return MCMC2;
	else if (mode == "map")
		return MAP;
	else {
		std::cerr << "Parallel mode \"" << mode << "\" not supported. Choose SERIAL, MCMC or MAP." << std::endl;
		return SERIAL;
	}
}

inline bool get_bool(const std::string& write)
{
	if (write == "yes" || write == "YES" || write == "Y" || write == "y" || write == "1")
		return true;
	else if (write == "no" || write == "NO" || write == "N" || write == "n" || write == "0")
		return false;
	else
		throw my_exception("Not sure what to write into the console. Choose yes or no in config.txt.");
}

inline void set_file_or_console(const std::string& input, std::vector<std::string>& files, std::vector<Console_output>& console_outputs, Console_output output)
{
	if (input != "no" && input != "n" && input != "0")				// no output file
	{
		if (input == "console" || input == "c")
		{
			console_outputs.push_back(output);
			files.push_back("console");
			return;
		}

		FILE* test = fopen(input.c_str(), "w");

		if (!test)
			std::cerr << "Cannot open file \"" << input << "\". Please check that all directories are created and accessible.\n";
		else
		{
			fclose(test);
			files.push_back(input);
		}
	}
}

inline void set_matrix_file(const std::string& input, std::string& filename)
{
	if (input != "no" && input != "n" && input != "0")				// no output file
	{
		if (input == "zero")
		{
			filename = "zero";
			return;
		}

		filename = input;
		FILE* test = fopen(filename.c_str(), "r");

		if (!test)
			std::cerr << "Cannot open file \"" << filename << "\". Please check that all directories are created and accessible.\n";
		else
			fclose(test);
	}
	else
		filename = "default";
}

inline std::vector<int> read_order(const std::string& custom_order)
{
	std::vector<int> ret;
	std::ifstream iFile(custom_order);
	int val;

	while (!iFile.eof())
	{
		iFile >> val;
		ret.push_back(val);
	}

	return ret;
}

Pattern* create_new_pattern(Matrix<bool>&& pattern, const Type type, const Map map, const Map_container map_container, Order order, std::vector<int> custom_order,
	const int N, const bool initialized, const Matrix<bool>& result, const int threads_count)
{
	std::chrono::system_clock::time_point start, end;
	std::vector<Counter> sizes;
	Matrix<bool> init;

	if (type == WALKING)
		return new Walking_pattern(pattern, N);
	else if (type == SLOW)
		return new Slow_pattern(pattern);
	else if (type == GENERAL && map_container == VECTOR)
	{
		if (order == AUTO)
		{
			if (!initialized)
			{
				init = Matrix<bool>(N, N);
				Patterns patterns;
				patterns.add(new General_pattern<std::vector<std::vector<int> > >(pattern, threads_count - 1, MAX, map));
				Matrix_Statistics matrix_stats(-1, 0, N, -1, false);
				Performance_Statistics perf_stats(1, -1);
				MCMCgenerator(N*N, patterns, init, perf_stats, matrix_stats, 1, -1);
			}

			General_pattern<std::vector<std::vector<int> > > gpSUM(pattern, threads_count - 1, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::vector<std::vector<int> > > gpMAX(pattern, threads_count - 1, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::vector<std::vector<int> > > gpDESC(pattern, threads_count - 1, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto desc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			if (desc_time < max_time && desc_time < sum_time)
				order = DESC;
			else if (max_time < sum_time)
				order = MAX;
			else
				order = SUM;
		}

		return new General_pattern<std::vector<std::vector<int> > >(pattern, threads_count - 1, order, map, std::move(custom_order));
	}
	else if (type == GENERAL && map_container == SET)
	{
		if (order == AUTO)
		{
			if (!initialized)
			{
				init = Matrix<bool>(N, N);
				Patterns patterns;
				patterns.add(new General_pattern<std::set<std::vector<int> > >(pattern, threads_count - 1, MAX, map));
				Matrix_Statistics matrix_stats(-1, 0, N, -1, false);
				Performance_Statistics perf_stats(1, -1);
				MCMCgenerator(N*N, patterns, init, perf_stats, matrix_stats, 1, -1);
			}

			General_pattern<std::set<std::vector<int> > > gpSUM(pattern, threads_count - 1, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::set<std::vector<int> > > gpMAX(pattern, threads_count - 1, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::set<std::vector<int> > > gpDESC(pattern, threads_count - 1, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto desc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			if (desc_time < max_time && desc_time < sum_time)
				order = DESC;
			else if (max_time < sum_time)
				order = MAX;
			else
				order = SUM;
		}

		return new General_pattern<std::set<std::vector<int> > >(pattern, threads_count - 1, order, map, std::move(custom_order));
	}
	else if (type == GENERAL && map_container == HASH)
	{
		if (order == AUTO)
		{
			if (!initialized)
			{
				init = Matrix<bool>(N, N);
				Patterns patterns;
				patterns.add(new General_pattern<std::unordered_set<std::vector<int>, int_vector_hasher> >(pattern, threads_count - 1, MAX, map));
				Matrix_Statistics matrix_stats(-1, 0, N, -1, false);
				Performance_Statistics perf_stats(1, -1);
				MCMCgenerator(N*N, patterns, init, perf_stats, matrix_stats, 1, -1);
			}

			General_pattern<std::unordered_set<std::vector<int>, int_vector_hasher> > gpSUM(pattern, threads_count - 1, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::unordered_set<std::vector<int>, int_vector_hasher> > gpMAX(pattern, threads_count - 1, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::unordered_set<std::vector<int>, int_vector_hasher> > gpDESC(pattern, threads_count - 1, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid((initialized ? result : init), -1, -1, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto desc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			if (desc_time < max_time && desc_time < sum_time)
				order = DESC;
			else if (max_time < sum_time)
				order = MAX;
			else
				order = SUM;
		}

		return new General_pattern<std::unordered_set<std::vector<int>, int_vector_hasher> >(pattern, threads_count - 1, order, map, std::move(custom_order));
	}
	else
		throw my_exception("This shouldn't have happened. Please send a message to the developer team with your config file.");

	return new General_pattern<std::vector<std::vector<int> > >(pattern, threads_count - 1, order, map, std::move(custom_order));
}

std::vector<Pattern_info> parse_config(std::istream& config, int& N, int& iter, int& random_seed, int& hist_from, int& hist_to, int& hist_freq,
	std::vector<Console_output>& console_outputs, std::vector<std::string>& output_files, std::vector<std::string>& hist_files,
	std::vector<std::string>& max_ones_files, std::vector<std::string>& csv_files, std::vector<std::string>& perf_files, std::string& init_matrix,
	int& threads_count, Parallel_mode& parallel_mode)
{
	std::string input, var, value;
	bool new_pattern = false;
	Pattern_info pattern;
	std::vector<Pattern_info> patterns;

	while (std::getline(config, input))
	{
		if (input.empty())
			continue;

		// comment starts with ';'
		if (input.front() == ';')
			continue;

		// if config contains CR and getline doesn't expect it
		if (input.back() == '\r')
			input.resize(input.size() - 1);

		if (input.empty())
			continue;

		if (input.front() == '\r')
			input = input.substr(1);

		if (input.empty())
			continue;

		if (input == "[pattern]")
		{
			// add the previously defined pattern if there was any
			if (new_pattern)
				patterns.push_back(pattern);

			pattern = Pattern_info();
			new_pattern = false;
		}

		size_t pos = input.find('=');

		if (pos == std::string::npos || pos == 0 || pos == input.size() - 1)
			continue;

		int prewhitespaces = 0;
		int postwhitespaces = 0;

		for (; input[pos - prewhitespaces - 1] == ' ' || input[pos - prewhitespaces - 1] == '\t'; ++prewhitespaces) { if (prewhitespaces >= pos) break; }
		for (; input[pos + postwhitespaces + 1] == ' ' || input[pos + postwhitespaces + 1] == '\t'; ++postwhitespaces) { if (postwhitespaces >= input.size() - pos - 1) break; }

		var = input.substr(0, pos - prewhitespaces);
		std::transform(var.begin(), var.end(), var.begin(), ::tolower);			// never seen this ::sth before
		value = input.substr(pos + 1 + postwhitespaces);
		std::transform(value.begin(), value.end(), value.begin(), ::tolower);	// just seen this ::sth two lines above for the first time

		// input
		if (var == "size" || var == "n")
			N = my_stoi(value);
		else if (var == "iter" || var == "iterations")
			iter = my_stoi(value);
		else if (var == "random_seed" || var == "seed")
		{
			if (value == "random" || value == "rand" || value == "r")
				random_seed = -1;
			else
				random_seed = my_stoi(value);
		}
		else if (var == "threads" || var == "threads_count")
		{
			threads_count = my_stoi(value);

			if (threads_count == -1)
				threads_count = (int)std::thread::hardware_concurrency();
		}
		else if (var == "parallel_mode" || var == "parallelism")
			parallel_mode = get_parallel_mode(value);
		else if (var == "initial_matrix_file" || var == "init_matrix")
			set_matrix_file(value, init_matrix);
		// pattern
		else if (var == "pattern_file" || var == "pattern")
		{
			set_matrix_file(value, pattern.pattern_file);
			new_pattern = true;
		}
		else if (var == "matrix_file" || var == "matrix")
			set_matrix_file(value, pattern.init_matrix_file);
		else if (var == "pattern_type" || var == "type")
			pattern.type = get_type(value);
		else if (var == "map_approach" || var == "approach")
			pattern.map = get_map(value);
		else if (var == "one_entries" || var == "map_one_entries")
			pattern.map.enough_entries = get_bool(value);
		else if (var == "recursion" || var == "map_recursion")
			pattern.map.recursion = get_bool(value);
		else if (var == "orthogonal_bounds" || var == "map_orthogonal_bounds")
			pattern.map.orthogonal_bounds = get_bool(value);
		else if (var == "map_container" || var == "container")
			pattern.map_container = get_map_container(value);
		else if (var == "line_order" || var == "order")
		{
			pattern.order = get_order(value);

			if (pattern.order == CUSTOM)						// file with custom order
				pattern.custom_order_file = value;
		}
		// output
		else if (var == "matrix_output" || var == "output")
			set_file_or_console(value, output_files, console_outputs, MATRIX);
		else if (var == "pattern_to_console")
		{
			try
			{
				if (get_bool(value))
					console_outputs.push_back(PATTERNS);
			}
			catch (const my_exception& exc)
			{
				std::cout << "While processing \"" << var << "\", exception triggered: " << exc.what() << std::endl;
			}
		}
		else if (var == "time_to_console" || var == "time")
		{
			try
			{
				if (get_bool(value))
					console_outputs.push_back(TIME);
			}
			catch (const my_exception& exc)
			{
				std::cout << "While processing \"" << var << "\", exception triggered: " << exc.what() << std::endl;
			}
		}
		else if (var == "performance_stats" || var == "performance")
			set_file_or_console(value, perf_files, console_outputs, PERF);
		else if (var == "performance_csv_stats" || var == "performance_csv")
			set_file_or_console(value, csv_files, console_outputs, CSV);
		// matrix stats
		else if (var == "histogram_file" || var == "histogram")
			set_file_or_console(value, hist_files, console_outputs, HIST);
		else if (var == "max_ones_file" || var == "max_ones" || var == "max_ones_matrix_file")
			set_file_or_console(value, max_ones_files, console_outputs, MAX_ONES);
		else if (var == "histogram_frequency" || var == "frequency")
			hist_freq = my_stoi(value);
		else if (var == "histogram_initial" || var == "initial")
			hist_from = my_stoi(value);
		else if (var == "histogram_final" || var == "final")
			hist_to = my_stoi(value);
		else
			std::cerr << "Couldn't parse the line \"" << var << " = " << value << "\" from config file. Make sure there are no spelling errors." << std::endl;
	}

	// add the last defined pattern if it hasn't been added already
	if (new_pattern)
		patterns.push_back(pattern);

	return patterns;
}

void set_patterns(Patterns& patterns, const std::vector<Pattern_info>& pattern_info, const Matrix<bool>& init_matrix, const bool initialized, const int N, const int threads_count)
{
	for (const auto& pattern : pattern_info)
	{
		std::vector<int> custom_order;

		if (pattern.custom_order_file != "")
			custom_order = read_order(pattern.custom_order_file);

		Matrix<bool> pat;

		try
		{
			pat = Matrix<bool>(pattern.pattern_file);
		}
		catch (...)
		{
			std::cerr << "Cannot read the input pattern from \"" << pattern.pattern_file << "\"." << std::endl;
			continue;
		}

		if (pattern.init_matrix_file == "default" || (pattern.init_matrix_file == "zero" && !initialized))
			patterns.add(create_new_pattern(std::move(pat), pattern.type, pattern.map, pattern.map_container, pattern.order, custom_order, N, initialized, init_matrix, threads_count));
		else
		{
			Matrix<bool> init;

			if (pattern.init_matrix_file == "zero")
				Matrix<bool>(N, N);
			else
			{
				try
				{
					Matrix<bool>(pattern.init_matrix_file);
				}
				catch (...)
				{
					std::cerr << "Cannot read the initial matrix for a pattern from \"" << pattern.init_matrix_file << "\"." << std::endl;
					continue;
				}
			}

			patterns.add(create_new_pattern(std::move(pat), pattern.type, pattern.map, pattern.map_container, pattern.order, custom_order, N, pattern.init_matrix_file != "zero", init, threads_count));
		}
	}
}

#endif