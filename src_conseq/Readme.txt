pdsat_conseq
// formula p*n^3 + i*n^2 + j*n + z (p - number of square, i - num of row, j - num of column, z - value )
[1.05] 02.12.2013 
	[+] Using several known rows
[1.04] 08.11.2013
	[+] Using own variant of Cartesian product - via next_cartesian()
[1.03] 31.10.2013 
	[*] PermutationWithRepetition() changed to MakePermutations()
[1.02] 03.10.2012
	[+] using known values of 2nd row in 1st square
[1.01] 31.08.2012 
	[+] checking of all possible variants for unknown cells in used rows. collecting of interrupted is not needed now
	[*] Work with order 10 added
	[-] dminisat and minisat2
	[*] row_index for every PermutationWithRepetition() was 1 (some correct vectors were skipped)
[1.00] 25.05.2012