#ifndef addit_func_h
#define addit_func_h

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/dynamic_bitset.hpp>

#ifdef _WIN32
#include "./win_headers/dirent.h"
#else
#include <dirent.h>
#endif

namespace Addit_func {

double cpu_time(void);
extern int strtoint( std::string str );
//extern string inttostr( int num );
//extern string doubletostr( double num );
extern bool isNumber( char num );
extern bool isNumberOrMinus( char num );
extern void cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, int &real_seconds );
extern void MakeCombinations(int n, int k, std::vector< std::vector<int> > &combinations, bool isOnlySizeCalc = false );
extern void MakePermutations( int n, int k, std::vector< std::vector<int> > &permutations );
extern int ConseqMultip( int low_bound, int high_bound );
extern int BitCount( unsigned u );
extern bool getdir( std::string dir, std::vector<std::string> &files );
//extern boost::dynamic_bitset<> IntVecToBitset( unsigned bitset_len, vector<int> &vec_int );
//extern vector<int> BitsetToIntVec( boost::dynamic_bitset<> &bs );
extern void shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count );
extern unsigned long long BitsetToUllong( boost::dynamic_bitset<> cur_bitset );
extern void UllongToBitset( unsigned long long ull, boost::dynamic_bitset<> &bs );
extern unsigned uint_rand( boost::random::mt19937 &gen );
extern bool bool_rand( boost::random::mt19937 &gen );
extern std::string exec( std::string cmd_str );

template< typename T > 
bool next_cartesian( std::vector<T> &vii, std::vector<int> &index_arr, T &cur_vi )
{
	if( index_arr.size() == 0 ) { // init
		index_arr.resize( vii.size() );
		//for( auto &x : index_arr )
		//	x = 0;
		for( std::vector<int> :: iterator it = index_arr.begin(); it != index_arr.end(); ++it )
			*it = 0;
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
		unsigned last_changable = (unsigned)index_arr.size()-1;
		while( last_changable != -1 ){
			if( index_arr[last_changable] < (int)vii[last_changable].size() - 1 )
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

