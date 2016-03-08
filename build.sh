#!/bin/sh
g++ -std=c++11 *.cpp EasyBMP/*.cpp -o matrix.exe -Wall -Wextra -pedantic -O3 -pthread -fopenmp