#include <bits/stdc++.h>

int main()
{
	std::vector<size_t> matrix(50*50, 0);
	std::ifstream iFile("changes.txt");
	
	size_t i, j;
	
	while(iFile >> i >> j)
		matrix[i * 50 + j] = matrix[i * 50 + j] ? 0 : 1;
	
	std::cout << 50 << " " << 50 << std::endl;
	for (size_t k = 0; k < 50 * 50; ++k)
		std::cout << matrix[k] << " ";
	
	return 0;
}