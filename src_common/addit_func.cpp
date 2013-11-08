#include "addit_func.h"

using namespace Addit_func;

bool IfValueInArray( int value, int *arr, int arr_len )
{
	for ( int i = 0; i < arr_len; i++ )
		if ( value == arr[i] ) 
			return true;
	return false;
}

int Addit_func :: conseq_multipl( int low_bound, int high_bound )
{
// multiplication low_bound * (low_bound + 1) * ... * high_bound
	int final_val = 1;
	for ( int i = low_bound; i <= high_bound; i++ )
		final_val *= i;
	return final_val;
}

void Addit_func :: MakeCombinations( int n, int k, vector< vector<int> > &combinations )
{
// Generation of set of all k-combinations of a set n
	int val;
	vector<int> index_arr;
	index_arr.resize(k);
	unsigned comb_index = 0;
	combinations.resize( conseq_multipl(n-k+1,n) / conseq_multipl(1, k) );
	for ( unsigned i = 0; i < combinations.size(); i++ )
		combinations[i].resize(k);
	for ( unsigned i = 0; i < index_arr.size(); i++ )
		index_arr[i] = i; // start indexes
	
	for ( ;; ) {
		combinations[comb_index++] = index_arr;
		if ( index_arr[k-1] != n-1 ) { // increase last value if we can
			index_arr[k-1]++;
			continue;
		}
		val = k-1;
		while ( ( val >= 0 ) && ( index_arr[val] == (n-k) + val ) )
			val--; // val - index of rightest value that is not on final position

		if ( val >= 0 ) { // if some values on final positions but not all
			index_arr[val]++;
			for ( int i = val+1; i < k; i++ )
				index_arr[i] = index_arr[i-1] + 1; // set initial state to all values right to changed one
		} // if val < 0 then final state - all '1' in tail
		else
			break; // all values are on final positions
	}
	index_arr.resize(0);
}

void Addit_func :: MakePermutations( int n, int k, vector< vector<int> > &permutations )
{
	vector< vector<int> > combinations;
	vector<int> cur_permutation;
	permutations.reserve( conseq_multipl(n-k+1,n) ); // reserve count of permutations
	MakeCombinations(n,k,combinations);
	for( unsigned i=0; i<combinations.size(); ++i ){
		cur_permutation = combinations[i];
		do permutations.push_back( cur_permutation );
		while ( next_permutation( cur_permutation.begin(), cur_permutation.end() ) );
	}
}