#ifndef API_cpp_
#define API_cpp_

#include "AvoidanceTests.hpp"
#include "Matrix.hpp"
#include "MCMC.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <time.h>

int main(int argc, char* argv[])
{
	size_t N = 15, iter = 100, type = 1;
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

	general_pattern pat3(pattern, N);
	try { walking_pattern pat1(pattern, N); } 
	catch (std::exception& e) { std::cout << e.what() << std::endl; }
	grandfather_pattern& rpat = pat3;

	clock_t t = clock();
	matrix<size_t> result(MCMCgenerator(N, iter, rpat));
	t = clock() - t;

	/*	If you want the result matrix to be in a text file
	std::ofstream oFile("output.txt");
	oFile << result.Print();
	oFile.close();
	*/

	std::cout << result.Print();
	std::cout << "\nAvoids:\n\n";
	std::cout << pattern.Print();
	std::cout << "\nRunning time: ";
	std::cout << (((float)t) / CLOCKS_PER_SEC);
	std::cout << " sec.\n";
	getchar();
	return 0;
}

#endif