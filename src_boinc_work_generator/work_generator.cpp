// work generator for SAT@home
// creating files for further creating workunit

#ifndef _WIN32
#include <mysql.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include "../src_common/common.h"
#include "../src_common/addit_func.h"

using namespace std;
using namespace Addit_func;

const long long MIN_WUS_FOR_CREATION = 100;

// Number of results we have received so far
static long long all_processed_wus;
static long long processed_wus;
static long long unsent_wus;
static long long running_wus;

char *pass_file_name = NULL;
bool IsTasksFile;
string prev_path;

unsigned long long assumptions_count = 0;
bool isRangeMode = false;

// Command line options

struct config_params_crypto {
	string problem_type;
	string settings_file;
	string data_file;
	unsigned long long cnf_variables;
	unsigned long long cnf_clauses;
	unsigned long long problems_in_wu;
	unsigned long long unsent_needed_wus;
	unsigned long long total_wus;
	unsigned long long created_wus;
};

static void print_help(const char *prog);
bool do_work( string master_config_file_name, vector<long long> &wu_id_vec );
void ParseConfigFile( string &cnf_head, string master_config_file_name, stringstream &config_sstream );
static void create_wus( string master_config_file_name, stringstream &config_sstream, config_params_crypto &config_p, 
					    string cnf_head, long long wus_for_creation_count, vector<long long> &wu_id_vec, bool &IsLastGenerating );
#ifndef _WIN32
void GetCountOfUnsentWUs( long long &unsent_count );
bool ProcessQuery( MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec );
#endif
//bool find_sat( int cnf_index );
//double cpuTime( void ) { return ( double )clock( ) / CLOCKS_PER_SEC; }

int main( int argc, char *argv[] )
{
	//double start_time = cpuTime();
	vector<long long> wu_id_vec;
	string str;
	IsTasksFile= false;
	
	// find full path to file
	string master_config_file_name = argv[1];
	std::cout << "master_config_file_name " << master_config_file_name << std::endl;
	int pos = -1, last_pos = 0;
	for(;;){
		pos = master_config_file_name.find("/", pos+1);
		if ( pos != string::npos )
			last_pos = pos;
		else
			break;
	}
	prev_path = master_config_file_name.substr(0, last_pos+1);
	cout << "prev_path " << prev_path << endl;
	cout << "master_config_file_name " << master_config_file_name << endl;
	
	do_work( master_config_file_name, wu_id_vec );
	//double total_time = cpuTime() - start_time;
	//cout << "total time " << total_time << endl;
	
	return 0;
}

void ParseConfigFile( config_params_crypto &config_p, string master_config_file_name, string &cnf_head, stringstream &config_sstream )
{
	fstream master_config_file;
	string input_str, str1, str2, str3;
	stringstream sstream;
;
	master_config_file.open( master_config_file_name.c_str() );
	
	cout << "In ParseConfigFile() master_config_file_name " << master_config_file_name << endl;
	if ( !master_config_file.is_open() ) {
		cerr << "file " << master_config_file_name << " doesn't exist" << endl;
		exit(1);
	}
	cout << endl << "input file opened" << endl;
	
	// problem_type
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.problem_type;
	sstream.str(""); sstream.clear();
	// settings_file
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.settings_file;
	sstream.str(""); sstream.clear();
	// data_file
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.data_file;
	sstream.str(""); sstream.clear();
	// cnf_variables
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_variables;
	sstream.str(""); sstream.clear();
	// cnf_clauses
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_clauses;
	sstream.str(""); sstream.clear();
	// problems_in_wu
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.problems_in_wu;
	sstream.str(""); sstream.clear();
	// unsent_needed_wus
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.unsent_needed_wus;
	sstream.str(""); sstream.clear();
	// total_wus
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.total_wus;
	sstream.str(""); sstream.clear();
	// created_wus
	getline( master_config_file, input_str );
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.created_wus;
	sstream.str(""); sstream.clear();
	
	// make cnf_head
	sstream << "p cnf " << config_p.cnf_variables << " " << config_p.cnf_clauses;
	cnf_head = sstream.str();
	sstream.str(""); sstream.clear();

	config_p.settings_file = prev_path + config_p.settings_file;
	config_p.data_file     = prev_path + config_p.data_file;
	
	cout << "problem_type "      << config_p.problem_type << endl;
	cout << "settings_file "     << config_p.settings_file << endl;
	cout << "data_file "         << config_p.data_file << endl;
	cout << "cnf_variables "     << config_p.cnf_variables << endl;
	cout << "cnf_clauses "       << config_p.cnf_clauses << endl;
	cout << "problems_in_wu "    << config_p.problems_in_wu << endl;
	cout << "unsent_needed_wus " << config_p.unsent_needed_wus << endl;
	cout << "total_wus "         << config_p.total_wus << endl;
	cout << "created_wus "       << config_p.created_wus << endl;
	cout << "*** cnf_head "      << cnf_head << endl;
	cout << endl;
	
	master_config_file.close();
}

bool do_work( string master_config_file_name, vector<long long> &wu_id_vec )
{
	//double start_time = cpuTime();

	processed_wus     = 0;
	all_processed_wus = 0;
	long long unsent_count = 0;
	config_params_crypto config_p;
	stringstream config_sstream;
	string cnf_head;
	long long wus_for_creation_count = 0;
	
	ParseConfigFile( config_p, master_config_file_name, cnf_head, config_sstream );
	ifstream ifile;
	ifile.open( config_p.data_file.c_str(), ios_base :: in | ios_base :: binary );
	if ( !ifile.is_open() ) {
		isRangeMode = true;
		cout << "isRangeMode " << isRangeMode << endl;
	}
	else
		ifile.close();
	bool IsLastGenerating = false;
	
	/*if ( IsTasksFile ) // get problems_in_wu assumptions for every task
		create_wus( ls, config_sstream, config_p, cnf_head, wus_for_creation_count, wu_id_vec, IsLastGenerating );
	else
	{*/
		// create wus - supply needed level of unsent wus
	wus_for_creation_count = 0;
	long long old_wus_for_creation_count = 0;
	for (;;) {
		if ( IsLastGenerating ) {
			if ( config_p.created_wus == config_p.total_wus ) {
				cout << "config_p.created_wus == config_p.total_wus" << endl;
				cout << config_p.created_wus << " == " << config_p.total_wus << endl;
			}
			cout << "IsLastGenerating " << IsLastGenerating << endl;
			cout << "generation done" << endl;
			break;
		}
#ifndef _WIN32
		GetCountOfUnsentWUs( unsent_count );
#endif
		cout << "unsent_count " << unsent_count << endl;
			
		if ( unsent_count < 0 ) {
			cout << "SQL error. unsent_count < 0. Waiting 60 seconds and try again" << endl;
#ifndef _WIN32
			sleep( 60 );
#endif
			continue; // try to execute SQL again
		}
				
		wus_for_creation_count = config_p.unsent_needed_wus - unsent_count;

		if ( wus_for_creation_count + config_p.created_wus >= config_p.total_wus ) {
			wus_for_creation_count = config_p.total_wus - config_p.created_wus; // create last several task
			IsLastGenerating = true;
			cout << "IsLastGenerating " << IsLastGenerating << endl;
		}

		cout << "wus_for_creation_count " << wus_for_creation_count << endl;

		if ( ( wus_for_creation_count >= MIN_WUS_FOR_CREATION ) || ( IsLastGenerating ) ) {
			// ls can be used many times - each launch vectore will be resized and filled again
			// ls.skip_valus is updated too
			create_wus( master_config_file_name, config_sstream, config_p, cnf_head, wus_for_creation_count, 
				        wu_id_vec, IsLastGenerating );
		}
		else {
			cout << "wus_for_creation_count < MIN_WUS_FOR_CREATION" << endl;
			cout << wus_for_creation_count << " < " << MIN_WUS_FOR_CREATION << endl;
			if ( old_wus_for_creation_count != wus_for_creation_count )
				cout << "wus_for_creation_count " << wus_for_creation_count << endl;
			old_wus_for_creation_count = wus_for_creation_count;
		}
			
		if ( !IsLastGenerating ) {
#ifndef _WIN32
			sleep( 1800 ); // wait
#endif
		}
	}
	//}
	
	cout << "wus_for_creation_count " << wus_for_creation_count << endl;
	
	//double total_time = cpuTime() - start_time;
	//cout << "total time " << total_time << endl;

	return 0;
}

void create_wus( string master_config_file_name, stringstream &config_sstream, config_params_crypto &config_p, 
				 string cnf_head, long long wus_for_creation_count, vector<long long> &wu_id_vec, bool &IsLastGenerating )
{
	ofstream output;
	string wu_tag_str;
	stringstream sstream, header_sstream;
	
	cout << "Start create_wus()" << endl;
	cout << "cnf_head " << cnf_head << endl;
	cout << "wus_for_creation_count " << wus_for_creation_count << endl;
	long long wu_index = 0;
	bool IsAddingWUneeded;
	bool IsFastExit = false;
	long long new_created_wus = 0;
	string str, word1;
	ifstream ifile;
	vector<int> var_choose_order;
	long long skip_byte;
	unsigned long long values_index;
	
	// read header data once - it's common for every wu
	ifile.open( config_p.settings_file.c_str() ); // write common head to every wu
	if ( !ifile.is_open() ) {
		cerr << "!ifile.is_open() " << config_p.settings_file << endl;
		exit(1);
	}
	cout << "file " << config_p.settings_file << " opened" << endl;
	//cout << "header:" << endl;
	unsigned header_str_count = 0;
	while ( getline( ifile, str ) ) {
		sstream << str;
		sstream >> word1;
		if ( word1 == "var_set" ) {
			int val;
			while ( sstream >> val )
				var_choose_order.push_back( val );
		}
		sstream.clear(); sstream.str("");
		header_sstream << str << endl;
		header_str_count++;
		//cout << str << endl;
	}
	ifile.close();
	cout << "header_str_count " << header_str_count << endl;
	short int si;
	unsigned long long ul = 0;
	
	if ( isRangeMode ) {
		cout << "isRangeMode" << endl;
		shl64( assumptions_count, var_choose_order.size() );
		cout << "var_choose_order.size() " << var_choose_order.size() << endl;
		cout << "assumptions_count " << assumptions_count << endl;
	}
	else {
		// count blocks of data in file
		if ( !assumptions_count ) {
			ifile.open( config_p.data_file.c_str(), ios_base :: in | ios_base :: binary );
			if ( !ifile.is_open() ) {
				cerr << "!ifile.is_open() " << config_p.data_file << endl;
				exit(1);
			}
			cout << "file " << config_p.data_file << " opened" << endl;
			/*ifile.read( (char*)&si, sizeof(si) ); // read header
			while ( ifile.read( (char*)&ul, sizeof(ul) ) ) {
				assumptions_count++;
				if ( assumptions_count % 10000000 == 0 )
					cout << assumptions_count << " time " << cpuTime() - assumption_counting_start_time << " s" << endl;
			}*/
			ifile.seekg( 0, ifile.end );
			long long total_byte_length = ifile.tellg();
			ifile.close();
			assumptions_count = (total_byte_length - 2) / sizeof(ul); // skip 2 byte of prefix 
			cout << "assumptions_count:" << endl;
			cout << assumptions_count << " time " << endl;
			ifile.close();
		}
		
		values_index = config_p.created_wus * config_p.problems_in_wu;
		skip_byte = 2 + sizeof(ul)*values_index;
		cout << "skip_byte " << skip_byte << endl;
		
		ifile.open( config_p.data_file.c_str(), ios_base :: in | ios_base :: binary );
		ifile.read( (char*)&si, sizeof(si) );
		// skip already sended values
		if ( skip_byte > 0 ) {
			ifile.clear();
			ifile.seekg( skip_byte, ifile.beg );
		}
	}
	
	cout << "created_wus "         << config_p.created_wus << endl;
	unsigned long long total_wu_data_count = ceil( double(assumptions_count) / double(config_p.problems_in_wu) );
	cout << "total_wu_data_count " << total_wu_data_count  << endl;
	if ( total_wu_data_count > config_p.total_wus )
		total_wu_data_count = config_p.total_wus;
	cout << "total_wu_data_count changed to " << total_wu_data_count << endl;
	
	cout << "before creating wus" << endl;
	unsigned long long now_created = 0;
	unsigned long long range1, range2;
	for( unsigned long long wu_index = config_p.created_wus; wu_index < config_p.created_wus + wus_for_creation_count; wu_index++ ) {
		if ( IsFastExit )
			break;
		output.open( "wu-input.txt", ios_base :: out );
		if ( !output.is_open() ) {
			std::cerr << "Failed to create wu-input.txt" << std::endl;
			exit(1);
		}
		output << header_sstream.str();
		
		if ( isRangeMode ) {
			if ( header_sstream.str().find( "before_range" ) == string::npos )
				output << "before_range" << endl; // add if forgot
			range1 = wu_index*config_p.problems_in_wu;
			IsAddingWUneeded = ( range1 < assumptions_count ) ? true : false;
			range2 = (wu_index+1)*config_p.problems_in_wu - 1;
			if ( range2 >= assumptions_count ) {
				range2 = assumptions_count - 1;
				cout << "range2 changed to " << range2 << endl;
				IsFastExit = true; // add last values to WU and exit
				IsLastGenerating = true; // tell to high-level function about ending of generation
			}
			output << range1 << " " << range2;
		}
		else {
			if ( header_sstream.str().find( "before_assignments" ) == string::npos )
				output << "before_assignments" << endl; // add if forgot
			output.close();
			output.open( "wu-input.txt", ios_base::out | ios_base::app | ios_base::binary );
			output.write( (char*)&si, sizeof(si) ); // write first 2 symbols
			IsAddingWUneeded = false; // if no values will be added then WU not needed
			for ( unsigned long long i = 0; i < config_p.problems_in_wu; i++ ) {
				if ( values_index >= assumptions_count ) {
					cout << "in create_wus() last data was added to WU" << endl;
					cout << "values_index " << values_index << endl;
					IsFastExit = true; // add last values to WU and exit
					IsLastGenerating = true; // tell to high-level function about ending of generation
					break;
				}
				ifile.read( (char*)&ul, sizeof(ul) );
				output.write( (char*)&ul, sizeof(ul) );
				values_index++;
				IsAddingWUneeded = true;
			}
		}
		output.close();
		
		if ( !IsAddingWUneeded ) {
			cout << "break cause of IsAddingWUneeded " << IsAddingWUneeded << endl;
			break; // don't create new WU
		}
		
		sstream.clear(); sstream.str("");
		sstream << config_p.problem_type;
		sstream << "--" << wu_index + 1; // save info about CNF name
		wu_tag_str = sstream.str();
		sstream.str( "" ); sstream.clear();
		if ( !now_created ) {
			cout << "isRangeMode " << isRangeMode << endl;
			cout << "first wu_tag_str " << wu_tag_str << endl;
		}
		now_created++;
		/*wu = DC_createWU( "pdsat_crypto", NULL, 0, wu_tag_str.c_str() );
		if ( !wu ) {
			DC_log( LOG_ERR, "Work unit creation has failed" );
			exit(1);
		}*/
		
		new_created_wus++;
	}
	ifile.close();
	
	cout << "new_created_wus " << new_created_wus << endl;
	config_p.created_wus += new_created_wus;
	
	if ( !IsTasksFile ) { // don't update if additional wus from file ewre created
		ofstream master_config_file;
		master_config_file.open( master_config_file_name.c_str() );
		master_config_file << config_sstream.str();
		master_config_file << "created_wus = " << config_p.created_wus << endl; // update master config file
		master_config_file.close();
	}

	cout << "created_wus " << config_p.created_wus << endl;
	cout << "---***---" << endl;
}

void add_result_to_file( string output_filename, char *tag, char *id )
{
	ofstream output_file;
	string final_output_name;
	stringstream final_sstream;
	stringstream sstream;
	
	sstream << "final_output.txt";
	final_output_name = sstream.str( );
	sstream.str( "" );
	sstream.clear();
	output_file.open( final_output_name.c_str(), ios_base :: out | ios_base :: app );
	if ( !output_file.is_open() ) {
		std::cerr << "Cannot open final_output.txt for writing" << std::cerr;
		exit(1);
	}
	
	string input_str;
	ifstream result_file;
	result_file.open( output_filename.c_str(), ios_base :: in );
	if ( !result_file.is_open( ) ) {
		std::cerr <<  "Cannot open result_file" << std::endl;
		output_file.close();
		exit(1);
	}
	
	final_sstream << "result_" << tag << " " << "id " << id << " ";
	while ( getline( result_file, input_str ) )
		final_sstream << input_str;
	final_sstream << endl;
	
	output_file << final_sstream.rdbuf();

	output_file.close();
	result_file.close();
}

#ifndef _WIN32
void GetCountOfUnsentWUs( long long &unsent_count )
{
	char *host = "localhost";
    char *db;
	char *user;
    char *pass;
	MYSQL *conn;
	
	ifstream pass_file;
	pass_file.open( pass_file_name );
	if ( !pass_file.is_open() ) {
		cerr << "psswd_file not open" << endl;
		exit(1);
	}
	string str;
	getline( pass_file, str );
	db = new char[str.length() + 1];
	strcpy( db, str.c_str() );
	getline( pass_file, str );
	user = new char[str.length() + 1];
	strcpy( user, str.c_str() );
	getline( pass_file, str );
	pass = new char[str.length() + 1];
	strcpy( pass, str.c_str() );
	cout << "db "   << db   << endl;
	cout << "user " << user << endl;
	//cout << "pass " << pass << endl;
	
	conn = mysql_init(NULL);
	if(conn == NULL)
		cerr << "Error: can't create MySQL-descriptor\n";

	if(!mysql_real_connect(conn, host, user, pass, db, 0, NULL, 0)) {
		cerr << "Error: can't connect to MySQL server" << endl;
		exit(1);
	}
	delete[] db;
	delete[] user;
	delete[] pass;

	vector< vector<stringstream *> > result_vec;
	str = "SELECT COUNT(*) FROM workunit WHERE id IN(SELECT workunitid FROM result WHERE server_state = 2)";
	cout << str << endl;

	if ( ProcessQuery( conn, str, result_vec ) ) {
		*result_vec[0][0] >> unsent_count;
		result_vec.clear();
		mysql_close(conn);
	}
	else
		unsent_count = -1;
}

bool ProcessQuery( MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec )
{
	MYSQL_RES *res;
	MYSQL_ROW row;
	int num_fields;
	
	if ( mysql_query(conn, str.c_str()) != 0 ) {
		cerr << "Error: can't execute SQL-query\n";
		return false;
	}
	
	res = mysql_store_result( conn );

	if( res == NULL ) 
		cerr << "Error: can't get the result description\n";

	num_fields = mysql_num_fields(res);
	stringstream *sstream_p;
	vector<stringstream *> result_data;

	if ( mysql_num_rows( res ) > 0 ) {
		while((row = mysql_fetch_row(res)) != NULL) {
			for( int i = 0; i < num_fields; ++i ) {
				sstream_p = new stringstream();
				*sstream_p << row[i]; // get value
				result_data.push_back( sstream_p );
			}
			result_vec.push_back( result_data );
			result_data.clear();
		}
	}

	mysql_free_result(res);

	return true;
}
#endif