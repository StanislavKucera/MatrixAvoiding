#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	size_t N = 15, iter = 5000, type = 1;
	if (argc > 1)
		N = (size_t)std::stoull(argv[1]);
	if (argc > 2)
		iter = (size_t)std::stoull(argv[2]);
	if (argc > 3)
		type = (size_t)std::stoull(argv[3]);

	matrix<size_t> pattern("input.txt");
	//matrix<size_t> pattern(2, 2, 1);
	//for (size_t i = 0; i < 1; i++)	{ pattern.at(1, i) = 0; }
	type += 0;

	general_pattern pat3(pattern, N);	// debugging purposes
	walking_pattern pat1(pattern, N);
	grandfather_pattern& rpat = pat1;

	matrix<size_t> result(MCMCgenerator(N, iter, rpat));
	std::cout << result.Print();
	std::cout << "\nAvoids:\n\n";
	std::cout << pattern.Print();
	getchar();
	return 0;
}

#endif