# MatrixAvoiding
Generator of matrices avoiding a given pattern

Giving a *binary matrix* pattern, N, number of iteration and pattern type,
program returns a random binary matrix N x N, which avoids the pattern.

Randomness of the result is provided by *Markov chain Monte Carlo* method.
Although noone knows how fast the method converges to random matrix,
you can specify the number of iterations as a input paramater.

*There are two types of patterns program can use:* <br />
  1. **General pattern:** Any square binary matrix is a general pattern. For given pattern the program finds the best possible order of line mapping and then tries to map lines of pattern to the lines of generated matrix in that order.
 
  2. **Walking pattern:** Square binary matrix is a walking pattern if there is a walk from upper-left corner to bottom-right corner of the matrix, which contains all one-entries of the matrix (but can also contain some zero-entries). For this type of pattern deciding whether generated matrix avoids the pattern is much faster, because it can be done by a simple algorithm of dynamic programming.
 
**Program input:** <br />
  The program expects the pattern to be in a text file "input" having format of two natural numbers "rows" and "columns" and then "rows" lines of "columns" numbers (0/1) separated by spaces. <br />
  It can be given up to five arguments when calling:
  1. "N" - size of the result matrix (N x N)
  2. "iterations" - number of iterations of Markov chain Monte Carlo method.
  3. "type" - type of a pattern as explained above:
    - 0 for a general pattern
    - 1 for a walking pattern
  4. "input" - name of the file from which pattern will be readed.
  5. "output" - name of the file which result will be written to.
  
**Default behaviour:** <br />
  MatrixAvoiding.exe 15 500 0 "input.txt" ""
   
**Program output:** <br />
  The program prints generated N x N matrix to the console as well as the pattern which the matrix is supposed to avoid. As a last thing it writes down the time the generating took. If "output" is not empty string it will also print generated matrix to the file.