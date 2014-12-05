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

#include "../src_common/addit_func.h"

using namespace Addit_func;

const long long MIN_WUS_FOR_CREATION = 100;

// Number of results we have received so far
static long long all_processed_wus;
static long long processed_wus;
static long long unsent_wus;
static long long running_wus;

std::string pass_file_name;
bool IsTasksFile;
std::string prev_path;

unsigned long long assumptions_count = 0;
bool isRangeMode = false;

// Command line options

struct config_params_crypto {
	std::string problem_type;
	std::string settings_file;
	std::string data_file;
	unsigned long long cnf_variables;
	unsigned long long cnf_clauses;
	unsigned long long problems_in_wu;
	unsigned long long unsent_needed_wus;
	unsigned long long total_wus;
	unsigned long long created_wus;
};

static void print_help(const char *prog);
bool do_work( std::string master_config_file_name, std::vector<long long> &wu_id_vec );
void ParseConfigFile( std::string &cnf_head, std::string master_config_file_name, std::stringstream &config_sstream );
static void create_wus( std::string master_config_file_name, std::stringstream &config_sstream, config_params_crypto &config_p, 
					    std::string cnf_head, long long wus_for_creation_count, std::vector<long long> &wu_id_vec, bool &IsLastGenerating );
#ifndef _WIN32
void GetCountOfUnsentWUs( long long &unsent_count );
bool ProcessQuery( MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec );
#endif
//bool find_sat( int cnf_index );
//double cpuTime( void ) { return ( double )clock( ) / CLOCKS_PER_SEC; }

int main( int argc, char *argv[] )
{
	//double start_time = cpuTime();
	std::vector<long long> wu_id_vec;
	std::string str;
	IsTasksFile= false;
	
	// find full path to file
	std::string master_config_file_name = argv[1];
	std::cout << "master_config_file_name " << master_config_file_name << std::endl;
	int pos = -1, last_pos = 0;
	for(;;){
		pos = master_config_file_name.find("/", pos+1);
		if ( pos != std::string::npos )
			last_pos = pos;
		else
			break;
	}
	prev_path = master_config_file_name.substr(0, last_pos+1);
	std::cout << "prev_path " << prev_path << std::endl;
	std::cout << "master_config_file_name " << master_config_file_name << std::endl;
	
	do_work( master_config_file_name, wu_id_vec );
	//double total_time = cpuTime() - start_time;
	//cout << "total time " << total_time << endl;
	
	return 0;
}

void ParseConfigFile( config_params_crypto &config_p, std::string master_config_file_name, std::string &cnf_head, 
					  std::stringstream &config_sstream )
{
	std::fstream master_config_file;
	std::string input_str, str1, str2, str3;
	std::stringstream sstream;
	master_config_file.open( master_config_file_name.c_str() );
	
	std::cout << "In ParseConfigFile() master_config_file_name " << master_config_file_name << std::endl;
	if ( !master_config_file.is_open() ) {
		std::cerr << "file " << master_config_file_name << " doesn't exist" << std::endl;
		exit(1);
	}
	std::cout << std::endl << "input file opened" << std::endl;
	
	// problem_type
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.problem_type;
	sstream.str(""); sstream.clear();
	// settings_file
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.settings_file;
	sstream.str(""); sstream.clear();
	// data_file
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.data_file;
	sstream.str(""); sstream.clear();
	// cnf_variables
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_variables;
	sstream.str(""); sstream.clear();
	// cnf_clauses
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_clauses;
	sstream.str(""); sstream.clear();
	// problems_in_wu
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.problems_in_wu;
	sstream.str(""); sstream.clear();
	// unsent_needed_wus
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.unsent_needed_wus;
	sstream.str(""); sstream.clear();
	// total_wus
	getline( master_config_file, input_str );
	config_sstream << input_str << std::endl;
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
	
	std::cout << "problem_type "      << config_p.problem_type << std::endl;
	std::cout << "settings_file "     << config_p.settings_file << std::endl;
	std::cout << "data_file "         << config_p.data_file << std::endl;
	std::cout << "cnf_variables "     << config_p.cnf_variables << std::endl;
	std::cout << "cnf_clauses "       << config_p.cnf_clauses << std::endl;
	std::cout << "problems_in_wu "    << config_p.problems_in_wu << std::endl;
	std::cout << "unsent_needed_wus " << config_p.unsent_needed_wus << std::endl;
	std::cout << "total_wus "         << config_p.total_wus << std::endl;
	std::cout << "created_wus "       << config_p.created_wus << std::endl;
	std::cout << "*** cnf_head "      << cnf_head << std::endl;
	std::cout << std::endl;
	
	master_config_file.close();
}

bool do_work( std::string master_config_file_name, std::vector<long long> &wu_id_vec )
{
	//double start_time = cpuTime();

	processed_wus     = 0;
	all_processed_wus = 0;
	long long unsent_count = 0;
	config_params_crypto config_p;
	std::stringstream config_sstream;
	std::string cnf_head;
	long long wus_for_creation_count = 0;
	
	ParseConfigFile( config_p, master_config_file_name, cnf_head, config_sstream );
	std::ifstream ifile;
	ifile.open( config_p.data_file.c_str(), std::ios_base :: in | std::ios_base :: binary );
	if ( !ifile.is_open() ) {
		isRangeMode = true;
		std::cout << "isRangeMode " << isRangeMode << std::endl;
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
				std::cout << "config_p.created_wus == config_p.total_wus" << std::endl;
				std::cout << config_p.created_wus << " == " << config_p.total_wus << std::endl;
			}
			std::cout << "IsLastGenerating " << IsLastGenerating << std::endl;
			std::cout << "generation done" << std::endl;
			break;
		}
#ifndef _WIN32
		GetCountOfUnsentWUs( unsent_count );
#endif
		std::cout << "unsent_count " << unsent_count << std::endl;
			
		if ( unsent_count < 0 ) {
			std::cout << "SQL error. unsent_count < 0. Waiting 60 seconds and try again" << std::endl;
#ifndef _WIN32
			sleep( 60 );
#endif
			continue; // try to execute SQL again
		}
				
		wus_for_creation_count = config_p.unsent_needed_wus - unsent_count;

		if ( wus_for_creation_count + config_p.created_wus >= config_p.total_wus ) {
			wus_for_creation_count = config_p.total_wus - config_p.created_wus; // create last several task
			IsLastGenerating = true;
			std::cout << "IsLastGenerating " << IsLastGenerating << std::endl;
		}

		std::cout << "wus_for_creation_count " << wus_for_creation_count << std::endl;

		if ( ( wus_for_creation_count >= MIN_WUS_FOR_CREATION ) || ( IsLastGenerating ) ) {
			// ls can be used many times - each launch vectore will be resized and filled again
			// ls.skip_valus is updated too
			create_wus( master_config_file_name, config_sstream, config_p, cnf_head, wus_for_creation_count, 
				        wu_id_vec, IsLastGenerating );
		}
		else {
			std::cout << "wus_for_creation_count < MIN_WUS_FOR_CREATION" << std::endl;
			std::cout << wus_for_creation_count << " < " << MIN_WUS_FOR_CREATION << std::endl;
			if ( old_wus_for_creation_count != wus_for_creation_count )
				std::cout << "wus_for_creation_count " << wus_for_creation_count << std::endl;
			old_wus_for_creation_count = wus_for_creation_count;
		}
			
		if ( !IsLastGenerating ) {
#ifndef _WIN32
			sleep( 1800 ); // wait
#endif
		}
	}
	//}
	
	std::cout << "wus_for_creation_count " << wus_for_creation_count << std::endl;
	
	//double total_time = cpuTime() - start_time;
	//cout << "total time " << total_time << endl;

	return 0;
}

void create_wus( std::string master_config_file_name, std::stringstream &config_sstream, config_params_crypto &config_p, 
				 std::string cnf_head, long long wus_for_creation_count, std::vector<long long> &wu_id_vec, bool &IsLastGenerating )
{
	std::ofstream output;
	std::string wu_tag_str;
	std::stringstream sstream, header_sstream;
	
	std::cout << "Start create_wus()" << std::endl;
	std::cout << "cnf_head " << cnf_head << std::endl;
	std::cout << "wus_for_creation_count " << wus_for_creation_count << std::endl;
	long long wu_index = 0;
	bool IsAddingWUneeded;
	bool IsFastExit = false;
	long long new_created_wus = 0;
	std::string str, word1;
	std::ifstream ifile;
	std::vector<int> var_choose_order;
	long long skip_byte;
	unsigned long long values_index;
	
	// read header data once - it's common for every wu
	ifile.open( config_p.settings_file.c_str() ); // write common head to every wu
	if ( !ifile.is_open() ) {
		std::cerr << "!ifile.is_open() " << config_p.settings_file << std::endl;
		exit(1);
	}
	std::cout << "file " << config_p.settings_file << " opened" << std::endl;
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
		header_sstream << str << std::endl;
		header_str_count++;
		//cout << str << endl;
	}
	ifile.close();
	std::cout << "header_str_count " << header_str_count << std::endl;
	short int si;
	unsigned long long ul = 0;
	
	if ( isRangeMode ) {
		std::cout << "isRangeMode" << std::endl;
		shl64( assumptions_count, var_choose_order.size() );
		std::cout << "var_choose_order.size() " << var_choose_order.size() << std::endl;
		std::cout << "assumptions_count " << assumptions_count << std::endl;
	}
	else {
		// count blocks of data in file
		if ( !assumptions_count ) {
			ifile.open( config_p.data_file.c_str(), std::ios_base :: in | std::ios_base :: binary );
			if ( !ifile.is_open() ) {
				std::cerr << "!ifile.is_open() " << config_p.data_file << std::endl;
				exit(1);
			}
			std::cout << "file " << config_p.data_file << " opened" << std::endl;
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
			std::cout << "assumptions_count:" << std::endl;
			std::cout << assumptions_count << " time " << std::endl;
			ifile.close();
		}
		
		values_index = config_p.created_wus * config_p.problems_in_wu;
		skip_byte = 2 + sizeof(ul)*values_index;
		std::cout << "skip_byte " << skip_byte << std::endl;
		
		ifile.open( config_p.data_file.c_str(), std::ios_base :: in | std::ios_base :: binary );
		ifile.read( (char*)&si, sizeof(si) );
		// skip already sended values
		if ( skip_byte > 0 ) {
			ifile.clear();
			ifile.seekg( skip_byte, ifile.beg );
		}
	}
	
	std::cout << "created_wus "         << config_p.created_wus << std::endl;
	unsigned long long total_wu_data_count = (unsigned long long)ceil( double(assumptions_count) / double(config_p.problems_in_wu) );
	std::cout << "total_wu_data_count " << total_wu_data_count  << std::endl;
	if ( total_wu_data_count > config_p.total_wus )
		total_wu_data_count = config_p.total_wus;
	std::cout << "total_wu_data_count changed to " << total_wu_data_count << std::endl;
	
	std::cout << "before creating wus" << std::endl;
	unsigned long long now_created = 0;
	unsigned long long range1, range2;
	for( unsigned long long wu_index = config_p.created_wus; wu_index < config_p.created_wus + wus_for_creation_count; wu_index++ ) {
		if ( IsFastExit )
			break;
		output.open( "wu-input.txt", std::ios_base :: out );
		if ( !output.is_open() ) {
			std::cerr << "Failed to create wu-input.txt" << std::endl;
			exit(1);
		}
		output << header_sstream.str();
		
		if ( isRangeMode ) {
			if ( header_sstream.str().find( "before_range" ) == std::string::npos )
				output << "before_range" << std::endl; // add if forgot
			range1 = wu_index*config_p.problems_in_wu;
			IsAddingWUneeded = ( range1 < assumptions_count ) ? true : false;
			range2 = (wu_index+1)*config_p.problems_in_wu - 1;
			if ( range2 >= assumptions_count ) {
				range2 = assumptions_count - 1;
				std::cout << "range2 changed to " << range2 << std::endl;
				IsFastExit = true; // add last values to WU and exit
				IsLastGenerating = true; // tell to high-level function about ending of generation
			}
			output << range1 << " " << range2;
		}
		else {
			if ( header_sstream.str().find( "before_assignments" ) == std::string::npos )
				output << "before_assignments" << std::endl; // add if forgot
			output.close();
			output.open( "wu-input.txt", std::ios_base::out | std::ios_base::app | std::ios_base::binary );
			output.write( (char*)&si, sizeof(si) ); // write first 2 symbols
			IsAddingWUneeded = false; // if no values will be added then WU not needed
			for ( unsigned long long i = 0; i < config_p.problems_in_wu; i++ ) {
				if ( values_index >= assumptions_count ) {
					std::cout << "in create_wus() last data was added to WU" << std::endl;
					std::cout << "values_index " << values_index << std::endl;
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
			std::cout << "break cause of IsAddingWUneeded " << IsAddingWUneeded << std::endl;
			break; // don't create new WU
		}
		
		sstream.clear(); sstream.str("");
		sstream << config_p.problem_type;
		sstream << "--" << wu_index + 1; // save info about CNF name
		wu_tag_str = sstream.str();
		sstream.str( "" ); sstream.clear();
		if ( !now_created ) {
			std::cout << "isRangeMode " << isRangeMode << std::endl;
			std::cout << "first wu_tag_str " << wu_tag_str << std::endl;
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
	
	std::cout << "new_created_wus " << new_created_wus << std::endl;
	config_p.created_wus += new_created_wus;
	
	if ( !IsTasksFile ) { // don't update if additional wus from file ewre created
		std::ofstream master_config_file;
		master_config_file.open( master_config_file_name.c_str() );
		master_config_file << config_sstream.str();
		master_config_file << "created_wus = " << config_p.created_wus << std::endl; // update master config file
		master_config_file.close();
	}

	std::cout << "created_wus " << config_p.created_wus << std::endl;
	std::cout << "---***---" << std::endl;
}

void add_result_to_file( std::string output_filename, char *tag, char *id )
{
	std::ofstream output_file;
	std::string final_output_name;
	std::stringstream final_sstream;
	std::stringstream sstream;
	
	sstream << "final_output.txt";
	final_output_name = sstream.str( );
	sstream.str( "" );
	sstream.clear();
	output_file.open( final_output_name.c_str(), std::ios_base :: out | std::ios_base :: app );
	if ( !output_file.is_open() ) {
		std::cerr << "Cannot open final_output.txt for writing" << std::cerr;
		exit(1);
	}
	
	std::string input_str;
	std::ifstream result_file;
	result_file.open( output_filename.c_str(), std::ios_base :: in );
	if ( !result_file.is_open( ) ) {
		std::cerr <<  "Cannot open result_file" << std::endl;
		output_file.close();
		exit(1);
	}
	
	final_sstream << "result_" << tag << " " << "id " << id << " ";
	while ( getline( result_file, input_str ) )
		final_sstream << input_str;
	final_sstream << std::endl;
	
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
	pass_file.open( pass_file_name.c_str() );
	if ( !pass_file.is_open() ) {
		std::cerr << "psswd_file not open" << std::endl;
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
	std::cout << "db "   << db   << std::endl;
	std::cout << "user " << user << std::endl;
	//cout << "pass " << pass << endl;
	
	conn = mysql_init(NULL);
	if(conn == NULL)
		std::cerr << "Error: can't create MySQL-descriptor\n" << std::endl;

	if(!mysql_real_connect(conn, host, user, pass, db, 0, NULL, 0)) {
		std::cerr << "Error: can't connect to MySQL server" << std::endl;
		exit(1);
	}
	delete[] db;
	delete[] user;
	delete[] pass;

	vector< vector<stringstream *> > result_vec;
	str = "SELECT COUNT(*) FROM workunit WHERE id IN(SELECT workunitid FROM result WHERE server_state = 2)";
	std::cout << str << std::endl;

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