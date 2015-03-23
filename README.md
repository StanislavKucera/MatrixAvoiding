# MatrixAvoiding
Generator of matrices avoiding a given pattern

Giving a binary matrix pattern, N, number of iteration and pattern type,
program returns a random binary matrix N x N, which avoids the pattern.

Randomness of the pattern is provided by Markov chain Monte Carlo method.
Although noone knows how fast the method converges to random matrix,
you can specify the number of iterations as a input paramater.

At the time being, only patterns for which this project can find random matrix
are those, which contain a walk from left-upper corner to right-bottom corner,
which covers all the one entries of the pattern.
