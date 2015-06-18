#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>
#include <time.h>

int main(int argc, char* argv[])
{
	size_t N = 15, iter = 500;
	enum Type { general, walking };
	Type type = general;
	std::string input = "input.txt", output = "";
	clock_t t = 0;

	// the size of resulting matrix is given
	if (argc > 1)
	{
		std::istringstream convert(argv[1]);
		convert >> N;
	}
	// the number of iterations is given
	if (argc > 2)
	{
		std::istringstream convert(argv[2]);
		convert >> iter;
	}
	// the type of the pattern is given
	if (argc > 3)
	{
		std::istringstream convert(argv[3]);
		short typ;
		convert >> typ;
		type = (Type)typ;
	}
	// the input pattern file is given
	if (argc > 4)
		input = argv[4];
	// the output result file is given
	if (argc > 5)
		output = argv[5];

	// create a pattern from the file - it has a constructor from initializer lists as well
	matrix<size_t> pattern(input);
	matrix<size_t> result;
	
	if (type == general)
	{
		general_pattern gp(pattern);
		t = clock();
		result = MCMCgenerator(N, iter, gp);
		t = clock() - t;
	}
	else if (type == walking)
	{
		// there is a possibility, that pattern is not a walking pattern
		try { 
			walking_pattern wp(pattern, N);
			t = clock();
			result = MCMCgenerator(N, iter, wp);
			t = clock() - t;
		}
		catch (std::invalid_argument& e) 
		{
			std::cout << e.what() << std::endl;
		}
	}

	// if output file is specified
	if (output != "")
	{
		std::ofstream oFile(output);
		// print the resulting matrix into it
		oFile << result.Print();
		oFile.close();
	}

	std::cout << result.Print();
	std::cout << "\nAvoiding pattern:\n\n";
	std::cout << pattern.Print();
	std::cout << "\nRunning time: ";
	std::cout << (((float)t) / CLOCKS_PER_SEC);
	std::cout << " sec.\n";
	getchar();
	return 0;
}

#endif