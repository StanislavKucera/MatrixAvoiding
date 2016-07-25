# MatrixAvoiding
A generator of matrices avoiding a given pattern as a submatrix

Given a *binary matrix* pattern, number N, number of iteration and pattern type,
program returns a random binary matrix N x N that avoids the pattern.

Randomness of the result is achieved by *Markov chain Monte Carlo* method.

There are two types of patterns program can use: <br />
  1. **General pattern:** Any binary matrix is a general pattern. For given pattern the program finds the best possible order of line mapping and then tries to map lines of pattern to the lines of generated matrix in that order.
 
  2. **Walking pattern:** A binary matrix is a walking pattern if there is a walk from upper-left corner to bottom-right corner of the matrix, which contains all one-entries of the matrix (but can also contain some zero-entries). For this type of pattern deciding whether generated matrix avoids the pattern is much faster, because it can be done by a simple algorithm of dynamic programming.

To speed up your computations, you can also use a parallel version of the program.
  
Before you use the program, you should read the user documentation, which is in Chapter 6 of Bc/thesis/thesis.pdf. As an example of a correct configuration file, use config.txt.

To find out more about the theory of Markov chains used in the program or about the program itself, read Bc/thesis/thesis.pdf.