#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>

void main(int argc, char** argv)
{
	size_t N = 15, iter = 50000, type = 1;
	matrix<int> pattern("input.txt");
	//matrix<int> pattern(2, 2, 1);
	//for (size_t i = 0; i < 1; i++)	{ pattern.at(1, i) = 0; }
	
	walking_pattern pat(pattern, N);
	matrix<int> result(MCMCgenerator(N, iter, pat));
	std::cout << result.Print();
	std::cout << "\nAvoids:\n\n";
	std::cout << pattern.Print();
	getchar();
}

#endif