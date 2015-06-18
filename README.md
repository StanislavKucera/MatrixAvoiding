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
  MatrixAvoiding.exe 20 1000 0 "input.txt" ""
   
**Program output:** <br />
  The program prints generated N x N matrix to the console as well as the pattern which the matrix is supposed to avoid. As a last thing it writes down the time the generating took. If "output" is not empty string it will also print generated matrix to the file.
  
**Algorithms:** <br />
  1. *General pattern* <br />
  The program takes one line of the pattern after another and tries to map them to every possible
  line of the resulting matrix. It is, as it sounds, a brute force method. To make it more efficient,
  the mappings, which have "important lines" mapped to the same lines of big matrix are shrinked into
  one mapping.
	- Line order finding:
	  1. Descending order - lines of the pattern are ordered according to the number of its one-entries
	  descendingly. This algorithm is trying to avoid as many mappings as possible, because it maps the
	  hardest lines to map at the beginning.
	  2. "DAG" order - lines of the pattern are ordered in a way that as many as possible lines can be
	  "forgotten" during the entire algorithm. This algorithm is trying to put together as many mappings 
	  as possible.
	- Finding what to remember:
	  Above, I used words "important lines" and "forget lines". Important line is a line of the pattern,
	  which bounds another not mapped line. This means either that those lines are parallel and
	  neighbouring each other or that they are perpendicular and intersect each other in a one-entry.
	  Algorithm itself is pretty easy, for given line order it just extends previous important lines by
	  the one added in current step and tries if any line can be forgotten.
	- Mapping a line of the pattern to a line of the big matrix:
	  First it needs to find bounds for the line according to lines, that are already mapped. 
	  For each line in between those bounds, it checks whether every one-entry of the line can be placed
	  to a one-entry on the line of the big matrix. This means that for one-entries on the intersection
	  with already mapped lines it checks the intersection of mapped lines and for those one-entries which
	  are on the intersection with a line that haven't been mapped yet, it finds its bounds, if there is
	  enough one-entries and if placing the line makes sense for the rest of already mapped lines.
	- Matrix avoiding:
	  Finds line order, what to remember when mapping each line and finds every possible mapping of the
	  subset of lines. If it success to find a mapping for the set of all lines, matrix doesn't avoid the
	  pattern. Otherwise it does.
  2. *Walking pattern:* <br />	
  The walk in the pattern is indexed from 1 to 2k. For each entry of a big matrix, there are two numbers c_h (horizontal) and c_v (vertical), which have the highest index of the walk that the part of the walk can be mapped into a submatrix from the left-upper corner to the entry. If there is 2k in any entry of the matrix then the whole pattern can be mapped to the matrix and the matrix doesn't avoid it.
    - c_v and c_h:
      For each entry of the pattern on the walk, the entry with index greater by one is either right to
	  the left of it (horizontal) or right underneath it (vertical). When computing c_v and c_h for an
	  entry of the matrix, it is sufficient to look at the closest element to the left and up and compute
	  c_v and c_h according to their values.	