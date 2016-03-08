#include <bits/stdc++.h>

using namespace std;

int main(int argc, char** argv)
{
	bool odd = false;
	const int size = 10;
	
	cout << size << " " << size << endl;
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			cout << ((odd) ? 1 : 0) << " ";
			odd = !odd;
		}
		if (size % 2 == 0)
		   odd = !odd;
		cout << endl;
	}
}
