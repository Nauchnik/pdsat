#include "addit_func.h"

using namespace Addit_func;

double Addit_func :: cpuTime(void) { 
	return (double)clock() / CLOCKS_PER_SEC;
}

int Addit_func :: strtoint( string str )
{
	stringstream sstream;
	sstream << str;
	int val;
	sstream >> val;
	return val;
}

/*string Addit_func :: inttostr( int num )
{
	stringstream sstream;
	sstream << num;
	string s = sstream.str();
	return s;
}

string Addit_func :: doubletostr( double num )
{
	stringstream sstream;
	sstream << num;
	string s = sstream.str();
	return s;
}*/

bool Addit_func :: isNumber(char num)
{ return num>='0' && num<='9'; }

bool Addit_func :: isNumberOrMinus(char num)
{ return ( ( num>='0' && num<='9' ) || ( num == '-' ) ); }

void Addit_func :: cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, int &real_seconds ) 
{
// Time of work in hours, minutes, seconds
	int full_minutes = ( int )full_seconds / 60;
	real_seconds = ( int )full_seconds % 60;
	real_hours   = full_minutes / 60;
	real_minutes = full_minutes % 60;
}

int Addit_func :: BitCount( unsigned u )
{
	unsigned uCount = u - ((u >> 1) & 033333333333) - ((u >> 2) & 011111111111);
    return ((uCount + (uCount >> 3)) & 030707070707) % 63;
}

int Addit_func :: ConseqMultip( int low_bound, int high_bound )
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
	combinations.resize( ConseqMultip(n-k+1,n) / ConseqMultip(1, k) );
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
	permutations.reserve( ConseqMultip(n-k+1,n) ); // reserve count of permutations
	MakeCombinations(n,k,combinations);
	for( unsigned i=0; i<combinations.size(); ++i ){
		cur_permutation = combinations[i];
		do permutations.push_back( cur_permutation );
		while ( next_permutation( cur_permutation.begin(), cur_permutation.end() ) );
	}
}

bool Addit_func :: getdir( string dir, vector<string> &files )
{
    DIR *dp;
	string cur_name;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        std::cout << std::endl << "Error in opening " << dir << std::endl;
        return false;
    }
    while ((dirp = readdir(dp)) != NULL) { 
		cur_name = string(dirp->d_name);
		if ( cur_name[0] != '.' ) files.push_back(cur_name); 
	}
    closedir(dp);
    return true;
}

/*
boost::dynamic_bitset<> Addit_func :: IntVecToBitset( unsigned bitset_len, vector<int> &int_vec )
{
	boost::dynamic_bitset<> bs( bitset_len );
	for ( unsigned i=0; i<int_vec.size(); i++ )
		bs.set( int_vec[i]-1 ); // set element to 1
	return bs;
}

vector<int> Addit_func :: BitsetToIntVec( boost::dynamic_bitset<> &bs )
{
	vector<int> vec_int;
	for ( unsigned i=0; i<bs.size(); i++ )
		if ( bs[i] ) vec_int.push_back( (int)(i+1) );
	return vec_int;
}*/

void Addit_func :: shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count )
{
	unsigned int val1, val2;
	if ( ( bit_count > 30 ) && ( bit_count < 61 ) ) {
		val1 = 30;
		val2 = bit_count - val1;
		val_for_left_shift =  ( unsigned long long int )( 1 << val1 );
		val_for_left_shift *= ( unsigned long long int )( 1 << val2 );
	}
	else if ( bit_count < 31 )
		val_for_left_shift =  ( unsigned long long int )( 1 << bit_count );
	else
		std::cout << "\n bit_count " <<  bit_count << " is too large ";
}

unsigned long long Addit_func :: BitsetToUllong( boost::dynamic_bitset<> cur_bitset )
{
	if ( cur_bitset.size() > 60 ) {
		std::cerr << "cur_bitset.size() > 60" << std::endl;
		exit(1);
	}
	unsigned long long ull = 0, val;
	for( unsigned i=0; i < cur_bitset.size(); ++i ) {
		if (cur_bitset[i]) { 
			val = 1;
			shl64(val,i);
			ull += val;
		}
	}
	return ull;
}

extern void Addit_func :: UllongToBitset( unsigned long long ull, boost::dynamic_bitset<> &bs )
{
	if ( bs.size() == 0 )
		bs.resize( 64 );
	bs.reset();
	unsigned index = 0;
	do {
		if ( index > bs.size() - 1 )
			bs.resize( index + 1 );
		if ( ull & 1 )
			bs[index] = 1;
		ull /= 2;
		index++;
	} while ( ull > 0 );
}

// TODO move to Addit_func namespace
extern unsigned Addit_func :: uint_rand( boost::random::mt19937 &gen ) {
	static boost::random::uniform_int_distribution<uint32_t> dist;
	return dist(gen);
}

extern bool Addit_func :: bool_rand( boost::random::mt19937 &gen ) {
	static boost::random::uniform_int_distribution<uint32_t> dist(0,1);
	return ( dist(gen) == 0 );
}

// ececute command via system process
extern std::string Addit_func :: exec( std::string cmd_str ) {
	char* cmd = new char[cmd_str.size() + 1];
	strcpy( cmd, cmd_str.c_str() );
	cmd[cmd_str.size()] = '\0';
#ifdef _WIN32
    FILE* pipe = _popen(cmd, "r");
#else
	FILE* pipe = popen(cmd, "r");
#endif
	delete[] cmd;
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
#ifdef _WIN32
    _pclose(pipe);
#else
	pclose(pipe);
#endif
    return result;
}