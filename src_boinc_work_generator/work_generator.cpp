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
#include <ctime>

const unsigned long long MIN_WUS_FOR_CREATION = 2;
const unsigned long long MAX_WUS_FOR_CREATION = 1000;

// Number of results we have received so far
static long long all_processed_wus;
static long long processed_wus;
static long long unsent_wus;
static long long running_wus;

std::string pass_file_name;
std::string master_config_file_name;
bool IsTasksFile;
std::string prev_path;

unsigned long long assumptions_count = 0;
bool isRangeMode = false;

double cpu_time(void) { return (double)clock() / CLOCKS_PER_SEC; }

void shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count )
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
	unsigned seconds_between_launches;
	unsigned credits_per_wu;
};

static void print_help(const char *prog);
bool do_work();
void parse_config_file( std::string &cnf_head, std::stringstream &config_without_created_wus_sstream );
static void create_wus( std::stringstream &cconfig_without_created_wus_sstream, config_params_crypto &config_p, std::string cnf_head, 
					    long long wus_for_creation_count, bool &IsLastGenerating );
#ifndef _WIN32
void GetCountOfUnsentWUs( unsigned long long &unsent_count );
bool ProcessQuery( MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec );
#endif

int main( int argc, char *argv[] )
{
	double start_time = cpu_time();
	std::string str;
	IsTasksFile= false;
	if ( argc < 3 ) {
		std::cerr << "Usage : program master_config_file_name pass_file_name" << std::endl;
		return 1;
	}
	// find full path to file
	master_config_file_name = argv[1];
	std::cout << "master_config_file_name " << master_config_file_name << std::endl;
	// read password for database
	pass_file_name = argv[2];
	std::cout << "pass_file_name " << pass_file_name << std::endl;
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
	std::cout << "new master_config_file_name " << master_config_file_name << std::endl;
	
	do_work();
	std::cout << "total time " << cpu_time() - start_time << std::endl;
	
	return 0;
}

void parse_config_file( config_params_crypto &config_p, std::string &cnf_head, std::stringstream &config_without_created_wus_sstream )
{
	std::fstream master_config_file;
	std::string input_str, str1, str2, str3;
	std::stringstream sstream;
	master_config_file.open( master_config_file_name.c_str() );
	
	std::cout << "In parse_config_file() master_config_file_name " << master_config_file_name << std::endl;
	if ( !master_config_file.is_open() ) {
		std::cerr << "file " << master_config_file_name << " doesn't exist" << std::endl;
		exit(1);
	}
	std::cout << std::endl << "input file opened" << std::endl;
	
	// problem_type
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.problem_type;
	sstream.str(""); sstream.clear();
	// settings_file
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.settings_file;
	sstream.str(""); sstream.clear();
	// data_file
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.data_file;
	sstream.str(""); sstream.clear();
	// cnf_variables
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_variables;
	sstream.str(""); sstream.clear();
	// cnf_clauses
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_clauses;
	sstream.str(""); sstream.clear();
	// problems_in_wu
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	if ( input_str.find( "problems_in_wu = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'problems_in_wu = '" << std::endl; 
		exit(1);
	}
	sstream >> str1 >> str2 >> config_p.problems_in_wu;
	sstream.str(""); sstream.clear();
	// unsent_needed_wus
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	if ( input_str.find( "unsent_needed_wus = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'unsent_needed_wus = '" << std::endl; 
		exit(1);
	}
	sstream >> str1 >> str2 >> config_p.unsent_needed_wus;
	sstream.str(""); sstream.clear();
	// read total_wus
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	if ( input_str.find( "total_wus = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'total_wus = '" << std::endl; 
		exit(1);
	}
	sstream >> str1 >> str2 >> config_p.total_wus;
	sstream.str(""); sstream.clear();
	// read seconds_between_launches
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;     
	if ( input_str.find( "seconds_between_launches = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'seconds_between_launches = '" << std::endl; 
		exit(1);
	}
	sstream >> str1 >> str2 >> config_p.seconds_between_launches;
	sstream.str(""); sstream.clear();
	// read credits_per_wu
	getline( master_config_file, input_str );
	config_without_created_wus_sstream << input_str << std::endl;
	sstream << input_str;
	if ( input_str.find( "credits_per_wu = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'credits_per_wu = '" << std::endl; 
		exit(1);
	}
	sstream >> str1 >> str2 >> config_p.credits_per_wu;
	sstream.str(""); sstream.clear();
	// created_wus
	getline( master_config_file, input_str );
	sstream << input_str;
	if ( input_str.find( "created_wus = " ) == std::string::npos ) {
		std::cerr << "string " << input_str << " doesn't include 'created_wus = '" << std::endl; 
		exit(1);
	}
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
	std::cout << "seconds_between_launches " << config_p.seconds_between_launches << std::endl;
	std::cout << "credits_per_wu "    << config_p.credits_per_wu << std::endl;
	std::cout << "created_wus "       << config_p.created_wus << std::endl;
	std::cout << "*** cnf_head "      << cnf_head << std::endl;
	std::cout << std::endl;
	
	master_config_file.close();
}

bool do_work()
{
	//double start_time = cpuTime();

	processed_wus     = 0;
	all_processed_wus = 0;
	unsigned long long unsent_count = 0;
	config_params_crypto config_p;
	std::stringstream config_without_created_wus_sstream;
	std::string cnf_head;
	unsigned long long wus_for_creation_count = 0;
	
	parse_config_file( config_p, cnf_head, config_without_created_wus_sstream );
	std::ifstream ifile;
	ifile.open( config_p.data_file.c_str(), std::ios_base :: in | std::ios_base :: binary );
	if ( !ifile.is_open() ) {
		isRangeMode = true;
		std::cout << "isRangeMode " << isRangeMode << std::endl;
	}
	else
		ifile.close();
	bool IsLastGenerating = false;
	
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

		if ( wus_for_creation_count > MAX_WUS_FOR_CREATION ) {
			std::cout << "wus_for_creation_count > MAX_WUS_FOR_CREATION" << std::endl;
			wus_for_creation_count = MAX_WUS_FOR_CREATION;
			std::cout << "changed to " << MAX_WUS_FOR_CREATION << std::endl;
		}

		if ( ( wus_for_creation_count >= MIN_WUS_FOR_CREATION ) || ( IsLastGenerating ) ) {
			// ls can be used many times - each launch vectore will be resized and filled again
			// ls.skip_valus is updated too
			create_wus( config_without_created_wus_sstream, config_p, cnf_head, wus_for_creation_count, IsLastGenerating );
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
			std::cout << "Waiting " << config_p.seconds_between_launches << " seconds" << std::endl;
			sleep( config_p.seconds_between_launches ); // wait
#endif
		}
	}
	
	std::cout << "wus_for_creation_count " << wus_for_creation_count << std::endl;
	
	//double total_time = cpuTime() - start_time;
	//cout << "total time " << total_time << endl;

	return 0;
}

void create_wus( std::stringstream &config_without_created_wus_sstream, config_params_crypto &config_p, 
				 std::string cnf_head, long long wus_for_creation_count, bool &IsLastGenerating )
{
	std::ofstream temp_wu_file_name;
	std::string cur_wu_input_file_name;
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
	unsigned long long ul = 0;
	
	if ( isRangeMode ) {
		std::cout << "isRangeMode" << std::endl;
		shl64( assumptions_count, var_choose_order.size() );
		std::cout << "var_choose_order.size() " << var_choose_order.size() << std::endl;
		std::cout << "assumptions_count " << assumptions_count << std::endl;
	}
	else {
		std::cerr << "isRangeMode " << isRangeMode << std::endl;
		exit(1);
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
	std::string system_str, wu_name;
	for( unsigned long long wu_index = config_p.created_wus; wu_index < config_p.created_wus + wus_for_creation_count; wu_index++ ) {
		if ( IsFastExit )
			break;
		
		range1 = wu_index*config_p.problems_in_wu;
		IsAddingWUneeded = ( range1 < assumptions_count ) ? true : false;
		range2 = (wu_index+1)*config_p.problems_in_wu - 1;
		if ( range2 >= assumptions_count ) {
			range2 = assumptions_count - 1;
			std::cout << "range2 changed to " << range2 << std::endl;
			IsFastExit = true; // add last values to WU and exit
			IsLastGenerating = true; // tell to high-level function about ending of generation
		}
		
		if ( !IsAddingWUneeded ) {
			std::cout << "break cause of IsAddingWUneeded " << IsAddingWUneeded << std::endl;
			break; // don't create new WU
		}
		
		sstream.clear(); sstream.str("");
		sstream << config_p.problem_type;
		sstream << "--" << wu_index + 1; // save info about CNF name

		wu_name = sstream.str();
		cur_wu_input_file_name = "input_" + wu_name;
		sstream.str( "" ); sstream.clear();
		
		temp_wu_file_name.open( "tmp_wu_file", std::ios_base :: out );
		if ( !temp_wu_file_name.is_open() ) {
			std::cerr << "Failed to create wu-input.txt" << std::endl;
			exit(1);
		}
		// write input data to WU file
		temp_wu_file_name << header_sstream.str();
		if ( header_sstream.str().find( "before_range" ) == std::string::npos )
			temp_wu_file_name << "before_range" << std::endl; // add if missed in header
		temp_wu_file_name << range1 << " " << range2;
		temp_wu_file_name.close(); 
		temp_wu_file_name.clear();
		
		system_str = "cp tmp_wu_file `dir_hier_path " + cur_wu_input_file_name + "`";
		std::cout << "before system command : " << system_str << std::endl; 
		system( system_str.c_str() );
		std::cout << "after system command" << std::endl;
		system_str = "create_work -appname pdsat_crypto -wu_name " + wu_name +
			         " -wu_template ./templates/workunit_template_bivium9.xml" + 
					 " -result_template ./templates/result_template_bivium9.xml " + cur_wu_input_file_name;
		std::cout << "before system command : " << system_str << std::endl;
		system( system_str.c_str() );
		std::cout << "after system command" << std::endl;
		
		if ( !now_created ) {
			std::cout << "isRangeMode " << isRangeMode << std::endl;
			std::cout << "first cur_wu_input_file_name " << cur_wu_input_file_name << std::endl;
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
	
	if ( !IsTasksFile ) { // don't update if additional wus from file were created
		std::ofstream master_config_file;
		master_config_file.open( master_config_file_name.c_str() );
		master_config_file << config_without_created_wus_sstream.str();
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
void GetCountOfUnsentWUs( unsigned long long &unsent_count )
{
	char *host = "localhost";
    char *db;
	char *user;
    char *pass;
	MYSQL *conn;
	
	std::ifstream pass_file;
	pass_file.open( pass_file_name.c_str() );
	if ( !pass_file.is_open() ) {
		std::cerr << "psswd_file not open" << std::endl;
		exit(1);
	}
	std::string str;
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

	std::vector< std::vector<std::stringstream *> > result_vec;
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

bool ProcessQuery( MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec )
{
	MYSQL_RES *res;
	MYSQL_ROW row;
	int num_fields;
	
	if ( mysql_query(conn, str.c_str()) != 0 ) {
		std::cerr << "Error: can't execute SQL-query\n";
		return false;
	}
	
	res = mysql_store_result( conn );

	if( res == NULL ) 
		std::cerr << "Error: can't get the result description\n";

	num_fields = mysql_num_fields(res);
	std::stringstream *sstream_p;
	std::vector<std::stringstream *> result_data;

	if ( mysql_num_rows( res ) > 0 ) {
		while((row = mysql_fetch_row(res)) != NULL) {
			for( int i = 0; i < num_fields; ++i ) {
				sstream_p = new std::stringstream();
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