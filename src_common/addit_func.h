#ifndef addit_func_h
#define addit_func_h

#include <iostream>
#include <vector>
#include <algorithm>
using std::vector;

namespace Addit_func {

void MakeCombinations( int n, int k, vector< vector<int> > &combinations );
extern void MakePermutations( int n, int k, vector< vector<int> > &permutations );
extern int combinations_count( int n, int k );
extern int conseq_multipl( int low_bound, int high_bound );
extern bool IfValueInArray( int value, int *arr, int arr_len);

template< typename T > 
bool next_cartesian( vector<T> &vii, vector<int> &index_arr, T &cur_vi )
{
	if( index_arr.size() == 0 ) { // init
		index_arr.resize( vii.size() );
		for( auto &x : index_arr )
			x = 0;
	}
	if( index_arr[0] == -1 )
		return false;
	// get current value
	cur_vi.resize(vii.size());
	for( unsigned i = 0; i < index_arr.size(); ++i )
		cur_vi[i] = vii[i][index_arr[i]];
		// check if last iteration
	bool IsLastValue = true; 
	for( unsigned i = 0; i < index_arr.size(); ++i ) {
		if ( index_arr[i] != vii[i].size() - 1 ) {
			IsLastValue = false;
			break;
		}
	}
	if( IsLastValue )	
		index_arr[0] = -1; // condition of stopping
	else {
		// find last changable row to increase its value
		unsigned last_changable = index_arr.size()-1;
		while( last_changable != -1 ){
			if( index_arr[last_changable] < vii[last_changable].size() - 1 )
				break;
			--last_changable;
		}
		index_arr[last_changable]++;
		for( unsigned i = last_changable+1; i < index_arr.size(); ++i )
			index_arr[i] = 0;
	}
	
	return true;
}

}

#endif

