# MatrixAvoiding
Generator of matrices avoiding a given pattern as a submatrix

Giving a *binary matrix* pattern, N, number of iteration and pattern type,
program returns a random binary matrix N x N, which avoids the pattern.

Randomness of the result is provided by *Markov chain Monte Carlo* method.
Although noone knows how fast the method converges to random matrix,
you can specify the number of iterations as a input paramater.

There are two types of patterns program can use: <br />
  1. **General pattern:** Any binary matrix is a general pattern. For given pattern the program finds the best possible order of line mapping and then tries to map lines of pattern to the lines of generated matrix in that order.
 
  2. **Walking pattern:** A binary matrix is a walking pattern if there is a walk from upper-left corner to bottom-right corner of the matrix, which contains all one-entries of the matrix (but can also contain some zero-entries). For this type of pattern deciding whether generated matrix avoids the pattern is much faster, because it can be done by a simple algorithm of dynamic programming.
 
**Program input:** <br />
  The program expects the pattern to be in a text file "input" having format of two natural numbers "rows" and "columns" and then "rows" lines of "columns" numbers (0/1) separated by spaces. <br />
  All the parameters are read from file config.txt having this format:
  1. Size of the generated matrix:	#size# (integer)
  2. Number of iterations of the generator:  #iterations# (integer)
  3. Pattern file: #pattern filename# (string)
  4. Type of the pattern:  #type# (string)
    - general
    - walking
  5. Initial big matrix:  #initial# (string)
    - "zero" for no initial matrix, algorithm will then start with a zero matrix
	- #initial matrix filename# from which algorithm will read the matrix; since the size is given in the argument before, file only contains the matrix itself
  6. Map function approach:  #map# (string)
    - recursion: map is called recursively for lines that intersect currently mapped line in a one-entry
	- compromise: there is no recursive call, but for those lines it atleast checks mandatory conditions
	- norecursion: it does not even check those mandatory conditions
  7. Container for storing mappings:  #container# (string)
    - vector: found mappings are stored in the std::vector
	- set: found mappings are stored in the std::set
  8. Line ordering:  #order# (string) [#order filename# (string)]
    - desc: lines are ordered according to the number of one-entries
	- max: lines are ordered in such a way that the maximum of remembered lines in one iteration is the smallest
	- sum: lines are ordered in such a way that the sum of remembered lines in all iteartions is the smallest
	- auto: all previous orderings are tested and the fastest one is taken
	- custom: order of lines is read from the file given in the next argument
  
**Program output:** <br />
  All parameters are read from file config.txt having this format:
  9. Output matrix file:  #output# (string)
    - "no" for no output file
	- #output filename# to which generated matrix will be written
  10. Write resulting matrix into console:  #console matrix# (string)
    - "yes": generated matrix will be written into the console
	- "no": generated matrix won't be written into the console
  11. Write the pattern into console:  #console pattern# (string)
    - "yes": the pattern will be written into the console
	- "no": the pattern won't be written into the console
  12. Write total time of the run into console:  #console time# (string)
    - "yes": the time spent in MCMC generator will be written into the console
	- "no": the time spent in MCMC generator won't be written into the console

**Default behaviour:** <br />
<table>
 <tr>
  <td>Program input:</td> <td></td>
 </tr>
 <tr>
  <td>Size of the generated matrix:</td> <td>40</td>
 </tr>
 <tr>
  <td>Number of iterations of the generator:</td> <td>10000</td>
 </tr>
 <tr>
  <td>Pattern file:</td> <td>input.txt</td>
 </tr>
 <tr>
  <td>Type of the pattern:</td> <td>general</td>
 </tr>
 <tr>
  <td>Initial big matrix:</td> <td>zero</td>
 </tr>
 <tr>
  <td>Map function approach:</td> <td>recursion</td>
 </tr>
 <tr>
  <td>Container for storing mappings:</td> <td>vector</td>
 </tr>
 <tr>
  <td>Line ordering:</td> <td>desc</td>
 </tr>
 <tr>
  <td></td> <td></td>
 </tr>
 <tr>
  <td>Program output:</td> <td></td>
 </tr>
 <tr>
  <td>Output matrix file:</td> <td>no</td>
 </tr>
 <tr>
  <td>Write resulting matrix into console:</td> <td>yes</td>
 </tr>
 <tr>
  <td>Write the pattern into console:</td> <td>yes</td>
 </tr>
 <tr>
  <td>Write total time of the run into console:</td> <td>yes</td>
 </tr>
</table>
     
**Algorithms:** <br />
<ol>
  <li>
  *General pattern* <br />
  The program takes one line of the pattern after another and tries to map them to every possible line of the resulting matrix. It is, as it sounds, a brute force method. To make it more efficient, the mappings, which have "important lines" mapped to the same lines of big matrix are shrinked into one mapping. <br />
	<ul>
	  <li>
	  Line order finding: <br />
	  <ol>
	    <li>
		Descending order - lines of the pattern are ordered according to the number of its
		one-entries descendingly. This algorithm is trying to avoid as many mappings as possible,
		because it maps the hardest lines to map at the beginning.
		</li>
		<li>
		"DAG" order - lines of the pattern are ordered in a way that as many as possible lines
		can be "forgotten" during the entire algorithm. This algorithm is trying to put together
		as many mappings as possible.
		</li>
	  </ol>
	  </li>
	  <li>
	  Finding what to remember: <br />
	  Above, I used words "important lines" and "forget lines". Important line is a line of the
	  pattern, which bounds another not mapped line. This means either that those lines are
	  parallel and neighbouring each other or that they are perpendicular and intersect each
	  other in a one-entry. Algorithm itself is pretty easy, for given line order it just
	  extends previous important lines by the one added in current step and tries if any line
	  can be forgotten.
	  </li>
	  <li>
	  Mapping a line of the pattern to a line of the big matrix: <br />
	  First it needs to find bounds for the line according to lines, that are already mapped.
	  For each line in between those bounds, it checks whether every one-entry of the line can
	  be placed to a one-entry on the line of the big matrix. This means that for one-entries
	  on the intersection with already mapped lines it checks the intersection of mapped lines
	  and for those one-entries which are on the intersection with a line that haven't been
	  mapped yet, it finds its bounds, if there is enough one-entries and if placing the line
	  makes sense for the rest of already mapped lines.
	  </li>
	  <li>
	  Matrix avoiding: <br />
	  Finds line order, what to remember when mapping each line and finds every possible
	  mapping of the subset of lines. If it success to find a mapping for the set of all
	  lines, matrix doesn't avoid the pattern. Otherwise it does.
	  </li>
	</ul>
  </li>
  <li>
  *Walking pattern:* <br />	
  The walk in the pattern is indexed from 1 to 2k. For each entry of a big matrix, there are
  two numbers c_h (horizontal) and c_v (vertical), which have the highest index of the walk
  that the part of the walk can be mapped into a submatrix from the left-upper corner to the
  entry. If there is 2k in any entry of the matrix then the whole pattern can be mapped to the
  matrix and the matrix doesn't avoid it. <br />
  <ul>
	<li>
    c_v and c_h: <br />
      For each entry of the pattern on the walk, the entry with index greater by one is either
	  right to the left of it (horizontal) or right underneath it (vertical). When computing
	  c_v and c_h for an entry of the matrix, it is sufficient to look at the closest element
	  to the left and up and compute c_v and c_h according to their values.	
	</li>
  </ul>
  </li>
</ol>