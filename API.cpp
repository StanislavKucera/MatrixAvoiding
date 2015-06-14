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

	if (argc > 1)
	{
		std::istringstream convert(argv[1]);
		convert >> N;
	}
	if (argc > 2)
	{
		std::istringstream convert(argv[2]);
		convert >> iter;
	}
	if (argc > 3)
	{
		std::istringstream convert(argv[3]);
		short typ;
		convert >> typ;
		type = (Type)typ;
	}
	if (argc > 4)
		input = argv[4];
	if (argc > 5)
		output = argv[5];

	matrix<size_t> pattern(input);
	//matrix<size_t> pattern(2, 2, 1);
	//for (size_t i = 0; i < 1; i++)	{ pattern.at(1, i) = 0; }
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

	if (output != "")
	{
		std::ofstream oFile(output);
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