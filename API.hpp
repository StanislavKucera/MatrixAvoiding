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

inline Type get_type(const std::string& type)
{
	if (type == "general" ||  type == "g")
		return GENERAL;
	else if (type == "walking" || type == "w")
		return WALKING;
	else if (type == "slow" || type == "s")
		return SLOW;
	else {
		throw my_exception("Pattern type not supported. Choose WALKING, GENERAL or SLOW.");
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
		throw my_exception("Mapping approach not supported. Choose SUPERACTIVE, ACTIVE, SEMIACTIVE, SEMILAZY, LAZY or SUPERLAZY.");
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
		throw my_exception("Map container not supported. Choose VECTOR, SET or HASH.");
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
	else if (mode == "mcmc" || mode == "m")
		return MCMC;
	else if (mode == "map" || mode == "general")
		return MAP;
	else {
		throw my_exception("Pattern type not supported. Choose SERIAL, MCMC or MAP.");
		return SERIAL;
	}
}

inline bool get_bool(const std::string& write)
{
	if (write == "yes" || write == "YES" || write == "Y" || write == "y")
		return true;
	else if (write == "no" || write == "NO" || write == "N" || write == "n")
		return false;
	else {
		throw my_exception("Not sure what to write into the console. Choose yes or no in config.txt.");
		return false;
	}
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
			std::cout << "Cannot open file \"" << filename << "\". Please check that all directories are created and accessible.\n";
		else
			fclose(test);
	}
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

void parse_config(std::istream& config, size_t& N, size_t& iter, size_t& hist_from, size_t& hist_to, size_t& hist_freq, Patterns& patterns, Matrix<size_t>& result,
	std::string& bmp_file, std::string& hist_file, std::string& max_ones_file, std::string& csv_file, std::string& perf_file,
	bool& console_time, bool& console_pattern, bool& console_matrix, bool& console_perf, bool& console_csv, bool& console_hist, bool& console_max_ones, bool& initialized,
	size_t& threads_count, Parallel_mode& parallel_mode)
{
	std::string input, var, value;
	bool new_pattern = false;
	Type type = GENERAL;
	Map map = SUPERACTIVE;
	Map_container map_container = HASH;
	Order order = AUTO;
	std::vector<size_t> custom_order;
	Matrix<size_t> pattern;

	while (std::getline(config, input))
	{
		if (input.empty())
			continue;

		if (input == "[pattern]")
		{
			// add the previously defined pattern if there was any
			if (new_pattern)
				patterns.add(create_new_pattern(std::move(pattern), type, map, map_container, order, custom_order, N, initialized, result));

			new_pattern = false;
			// set default values
			type = GENERAL;
			map = SUPERACTIVE;
			map_container = HASH;
			order = AUTO;
			custom_order.clear();
		}

		size_t pos = input.find('=');

		if (pos == std::string::npos)
			continue;

		var = input.substr(0, (input[pos - 1] == ' ' ? pos - 1 : pos));
		std::transform(var.begin(), var.end(), var.begin(), ::tolower);			// never seen this ::sth before
		value = input.substr((input[pos + 1] == ' ' ? pos + 2 : pos + 1));
		std::transform(value.begin(), value.end(), value.begin(), ::tolower);	// only seen this ::sth two lines above for the first time

		// input
		if (var == "size" || var == "n")
			N = std::stoul(value);
		else if (var == "iter" || var == "iterations")
			iter = std::stoul(value);
		else if (var == "threads" || var == "threads_count")
			threads_count = std::stoul(value);
		else if (var == "parallel_mode" || var == "parallelism")
			parallel_mode = get_parallel_mode(value);
		else if (var == "matrix_file" || var == "matrix")
		{
			if (value == "zero")
				result = Matrix<size_t>(N, N);			// zero matrix NxN
			else {
				result = Matrix<size_t>(N, N, value);
				initialized = true;
			}
		}
		// pattern
		else if (var == "pattern_file" || var == "pattern")
		{
			pattern = Matrix<size_t>(value);
			new_pattern = true;
		}
		else if (var == "pattern_type" || var == "type")
			type = get_type(value);
		else if (var == "map_approach" || var == "approach")
			map = get_map(value);
		else if (var == "map_container" || var == "container")
			map_container = get_map_container(value);
		else if (var == "line_order" || var == "order")
		{
			order = get_order(value);

			if (order == CUSTOM)						// file with custom order
				custom_order = read_order(value);
		}
		// output
		else if (var == "matrix_output" || var == "output")
			set_file_or_bool(value, bmp_file, console_matrix);
		//else if (var == "pattern_to_console")
		//	console_pattern = get_bool(value);
		else if (var == "time_to_console" || var == "time")
			console_time = get_bool(value);
		else if (var == "performance_stats" || var == "performance")
			set_file_or_bool(value, perf_file, console_perf);
		else if (var == "performace_csv_stats" || var == "performace_csv")
			set_file_or_bool(value, csv_file, console_csv);
		// matrix stats
		else if (var == "histogram_file" || var == "histogram")
			set_file_or_bool(value, hist_file, console_hist);
		else if (var == "max_ones_file" || var == "max_ones")
			set_file_or_bool(value, max_ones_file, console_max_ones);
		else if (var == "histogram_frequency" || var == "frequency")
			hist_freq = std::stoul(value);
		else if (var == "histogram_initial" || var == "initial")
			hist_from = std::stoul(value);
		else if (var == "histogram_final" || var == "final")
		{
			hist_to = std::stoul(value);

			if (hist_to == 0)
				hist_to = (size_t)-1;
		}
	}

	// add the last defined pattern if it hasn't been added already
	if (new_pattern)
		patterns.add(create_new_pattern(std::move(pattern), type, map, map_container, order, custom_order, N, initialized, result));
}

#endif