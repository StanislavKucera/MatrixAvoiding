#ifdef _MSC_VER 
// MS Visual Studio gives warnings when using fopen.
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef API_hpp_
#define API_hpp_

#include "HelpFunctionsAndStructures.hpp"
#include "Matrix.hpp"
#include "PatternHeaders.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>

struct Pattern_info
{
	Pattern_info() : pattern_file(""), init_matrix_file("default"), custom_order_file(""), type(GENERAL), map(SUPERACTIVE), map_container(HASH), order(AUTO) {}

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
		return SUPERACTIVE;
	else if(map == "4" || map == "active")
		return ACTIVE;
	else if (map == "3" || map == "semi active" || map == "semi_active" || map == "semiactive")
		return SEMIACTIVE;
	else if (map == "2" || map == "semi lazy" || map == "semi_lazy" || map == "semilazy")
		return SEMILAZY;
	else if (map == "1" || map == "lazy")
		return ACTIVE;
	else if (map == "0" || map == "super lazy" || map == "super_lazy" || map == "superlazy")
		return SEMIACTIVE;
	else {
		std::cerr << "Mapping approach \"" << map << "\" not supported. Choose SUPERACTIVE, ACTIVE, SEMIACTIVE, SEMILAZY, LAZY or SUPERLAZY." << std::endl;
		return SUPERACTIVE;
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
	else if (mode == "map")
		return MAP;
	else {
		std::cerr << "Parallel mode \"" << mode << "\" not supported. Choose SERIAL, MCMC or MAP." << std::endl;
		return SERIAL;
	}
}

inline bool get_bool(const std::string& write)
{
	if (write == "yes" || write == "YES" || write == "Y" || write == "y")
		return true;
	else if (write == "no" || write == "NO" || write == "N" || write == "n")
		return false;
	else
		throw my_exception("Not sure what to write into the console. Choose yes or no in config.txt.");
}

inline void set_file_or_bool(const std::string& input, std::string& filename, bool& console)
{
	if (input != "no" && input != "n")				// no output file
	{
		if (input == "console" || input == "c")
		{
			console = true;
			return;
		}

		filename = input;
		FILE* test = fopen(filename.c_str(), "w");

		if (!test)
			std::cerr << "Cannot open file \"" << filename << "\". Please check that all directories are created and accessible.\n";
		else
			fclose(test);
	}
}

inline void set_matrix_file(const std::string& input, std::string& filename)
{
	if (input != "no" && input != "n")				// no output file
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

inline std::vector<size_t> read_order(const std::string& custom_order)
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

Pattern* create_new_pattern(Matrix<size_t>&& pattern, const Type type, const Map map, const Map_container map_container, Order order, std::vector<size_t> custom_order,
	const size_t N, const bool initialized, const Matrix<size_t>& result)
{
	std::chrono::system_clock::time_point start, end;
	std::vector<Counter> sizes;

	if (type == WALKING)
		return new Walking_pattern(pattern, N);
	else if (type == SLOW)
		return new Slow_pattern(pattern);
	else if (type == GENERAL && map_container == VECTOR)
	{
		if (order == AUTO)
		{
			if (!initialized)
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

			General_pattern<std::vector<std::vector<size_t> > > gpSUM(pattern, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::vector<std::vector<size_t> > > gpMAX(pattern, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::vector<std::vector<size_t> > > gpDESC(pattern, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid(result, sizes))
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

		return new General_pattern<std::vector<std::vector<size_t> > >(pattern, order, map, std::move(custom_order));
	}
	else if (type == GENERAL && map_container == SET)
	{
		if (order == AUTO)
		{
			if (!initialized)
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

			General_pattern<std::set<std::vector<size_t> > > gpSUM(pattern, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::set<std::vector<size_t> > > gpMAX(pattern, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::set<std::vector<size_t> > > gpDESC(pattern, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid(result, sizes))
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

		return new General_pattern<std::set<std::vector<size_t> > >(pattern, order, map, std::move(custom_order));
	}
	else if (type == GENERAL && map_container == HASH)
	{
		if (order == AUTO)
		{
			if (!initialized)
				Matrix<size_t> rand = Matrix<size_t>::random_bin_matrix(N, N);

			General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpSUM(pattern, SUM, map);
			start = std::chrono::system_clock::now();
			if (!gpSUM.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto sum_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpMAX(pattern, MAX, map);
			start = std::chrono::system_clock::now();
			if (!gpMAX.avoid(result, sizes))
				throw my_exception("Initial big matrix does not avoid the pattern");
			end = std::chrono::system_clock::now();
			auto max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

			General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> > gpDESC(pattern, DESC, map);
			start = std::chrono::system_clock::now();
			if (!gpDESC.avoid(result, sizes))
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

		return new General_pattern<std::unordered_set<std::vector<size_t>, size_t_vector_hasher> >(pattern, order, map, std::move(custom_order));
	}
	else
		throw my_exception("This shouldn't have happened. Please send a message to the developer team with your config file.");

	return new General_pattern<std::vector<std::vector<size_t> > >(pattern, order, map, std::move(custom_order));
}

std::vector<Pattern_info> parse_config(std::istream& config, size_t& N, size_t& iter, size_t& hist_from, size_t& hist_to, size_t& hist_freq,
	std::string& bmp_file, std::string& hist_file, std::string& max_ones_file, std::string& csv_file, std::string& perf_file, std::string& init_matrix,
	bool& console_time, bool& console_pattern, bool& console_matrix, bool& console_perf, bool& console_csv, bool& console_hist, bool& console_max_ones,
	size_t& threads_count, Parallel_mode& parallel_mode)
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

		var = input.substr(0, (input[pos - 1] == ' ' ? pos - 1 : pos));
		std::transform(var.begin(), var.end(), var.begin(), ::tolower);			// never seen this ::sth before
		value = input.substr((input[pos + 1] == ' ' ? pos + 2 : pos + 1));
		std::transform(value.begin(), value.end(), value.begin(), ::tolower);	// just seen this ::sth two lines above for the first time

		// input
		if (var == "size" || var == "n")
			N = you_must_be_kidding_me::stoul(value);
		else if (var == "iter" || var == "iterations")
			iter = you_must_be_kidding_me::stoul(value);
		else if (var == "threads" || var == "threads_count")
			threads_count = you_must_be_kidding_me::stoul(value);
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
			set_file_or_bool(value, bmp_file, console_matrix);
		else if (var == "pattern_to_console")
		{
			try
			{
				console_pattern = get_bool(value);
			}
			catch (const my_exception& exc)
			{
				std::cout << "While processing \"" << var << "\", exception triggered: " << exc.what() << std::endl;
				console_pattern = false;
			}
		}
		else if (var == "time_to_console" || var == "time")
		{
			try
			{
				console_time = get_bool(value);
			}
			catch (const my_exception& exc)
			{
				std::cout << "While processing \"" << var << "\", exception triggered: " << exc.what() << std::endl;
				console_time = false;
			}
		}
		else if (var == "performance_stats" || var == "performance")
			set_file_or_bool(value, perf_file, console_perf);
		else if (var == "performance_csv_stats" || var == "performance_csv")
			set_file_or_bool(value, csv_file, console_csv);
		// matrix stats
		else if (var == "histogram_file" || var == "histogram")
			set_file_or_bool(value, hist_file, console_hist);
		else if (var == "max_ones_file" || var == "max_ones" || var == "max_ones_matrix_file")
			set_file_or_bool(value, max_ones_file, console_max_ones);
		else if (var == "histogram_frequency" || var == "frequency")
			hist_freq = you_must_be_kidding_me::stoul(value);
		else if (var == "histogram_initial" || var == "initial")
			hist_from = you_must_be_kidding_me::stoul(value);
		else if (var == "histogram_final" || var == "final")
		{
			hist_to = you_must_be_kidding_me::stoul(value);

			if (hist_to == 0)
				hist_to = (size_t)-1;
		}
		else
			std::cerr << "Couldn't parse the line \"" << var << " = " << value << "\" from config file. Make sure there are no spelling errors." << std::endl;
	}

	// add the last defined pattern if it hasn't been added already
	if (new_pattern)
		patterns.push_back(pattern);

	return patterns;
}

void set_patterns(Patterns& patterns, const std::vector<Pattern_info>& pattern_info, const Matrix<size_t>& init_matrix, const bool initialized, const size_t N)
{
	for (const auto& pattern : pattern_info)
	{
		std::vector<size_t> custom_order;

		if (pattern.custom_order_file != "")
			custom_order = read_order(pattern.custom_order_file);

		Matrix<size_t> pat;

		try
		{
			pat = Matrix<size_t>(pattern.pattern_file);
		}
		catch (...)
		{
			std::cerr << "Cannot read the input pattern from \"" << pattern.pattern_file << "\"." << std::endl;
			continue;
		}

		if (pattern.init_matrix_file == "default" || (pattern.init_matrix_file == "zero" && !initialized))
			patterns.add(create_new_pattern(std::move(pat), pattern.type, pattern.map, pattern.map_container, pattern.order, custom_order, N, initialized, init_matrix));
		else
		{
			Matrix<size_t> init;

			if (pattern.init_matrix_file == "zero")
				Matrix<size_t>(N, N);
			else
			{
				try
				{
					Matrix<size_t>(N, N, pattern.init_matrix_file);
				}
				catch (...)
				{
					std::cerr << "Cannot read the initial matrix for a pattern from \"" << pattern.init_matrix_file << "\"." << std::endl;
					continue;
				}
			}

			patterns.add(create_new_pattern(std::move(pat), pattern.type, pattern.map, pattern.map_container, pattern.order, custom_order, N, pattern.init_matrix_file != "zero", init));
		}
	}
}

#endif