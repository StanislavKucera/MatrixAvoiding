#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	size_t N = 15, iter = 50000, type = 1;
	if (argc > 1)
		N = std::stoull(argv[1]);
	if (argc > 2)
		iter = std::stoull(argv[2]);
	if (argc > 3)
		type = std::stoull(argv[3]);

	matrix<int> pattern("input.txt");
	//matrix<int> pattern(2, 2, 1);
	//for (size_t i = 0; i < 1; i++)	{ pattern.at(1, i) = 0; }
	type += 0;

	grandfather_pattern* ppat;
	walking_pattern pat1(pattern, N);
	ppat = &pat1;

	matrix<int> result(MCMCgenerator(N, iter, ppat));
	std::cout << result.Print();
	std::cout << "\nAvoids:\n\n";
	std::cout << pattern.Print();
	getchar();
	return 0;
}

#endif