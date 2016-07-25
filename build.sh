#!/bin/sh

NAMEEXE="matrix.exe"
g++ -std=c++11 source/*.cpp source/EasyBMP/*.cpp -o $NAMEEXE -Wall -Wextra -pedantic -O3 -pthread
