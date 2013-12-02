// master DC-API PD-SAT
// 04.10.2013 - added limitaion in using of positive_literals; if values ended then stop creating new WUs
// 03.10.2013 - added using of known 2nd row for 1s square
//            - in WU file fixed c ls_10_3_inc70
// 30.08.2013 - added generating of WUs with particular numbers (file with flag -w)
// 19.07.2013 - added generating of triple inc60
// 01.07.2013 - added writing problem_type (diag, inc0, inc70) into ls.problem_count
// 30.06.2013 - added generating of triple inc70
// 23.06.2013 - added new experiment: latin_triple
// 01.03.2013 - interval of SQL launch was increased from 5 min to 1 hour
// if SQL fails with error, skip it and try again

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mysql.h>

#include <dc.h>
#include <dc_client.h>

#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "../src_common/latin_squares.h"
#include "../src_common/common.h"

using namespace std;

const int PROBLEMS_IN_WU = 20;
const int MIN_WUS_FOR_CREATION = 200;

// Number of results we have received so far
static int all_processed_wus;
static int processed_wus;
static int unsent_wus;
static int running_wus;

char *dcapi_config_file_name = NULL;
char *master_config_file_name = NULL;
char *wu_id_file_name = NULL; // id of wu for additional creation
bool IsTasksFile;

// Command line options

struct config_params {
	string cnf_name;
	int cnf_variables;
	int cnf_clauses;
	int N;
	int rows_count;
	int K;
	int unsent_needed_wus;
	int total_wus;
	int created_wus;
	int skip_values;
	vector< vector<int>> known_values_vec;
};

static void print_help(const char *prog);
bool do_work( vector<int> &wu_id_vec );
void ParseConfigFile( string &cnf_name, string &cnf_head, stringstream &config_sstream );
static void create_wus( latin_square &ls, stringstream &config_sstream, config_params &config_p, 
	                    string cnf_head, int wus_for_creation_count, vector<int> &wu_id_vec, bool &IsLastGenerating );
bool find_sat( int cnf_index );
void add_result_to_file( string output_filename, char *tag, char *id );
void GetCountOfUnsentWUs( int &unsent_count );
bool ProcessQuery( MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec );

static const struct option longopts[] =
{
	{ "dcapi_config",	 required_argument, NULL, 'c' },
	{ "master_config",   required_argument, NULL, 'm' },
	{ "master_config",   optional_argument, NULL, 'w' },
	{ "help",	         no_argument,		NULL, 'h' },
	{ NULL }
};

int main( int argc, char *argv[] )
{
	int c;
	double start_time = cpuTime();
	vector<int> wu_id_vec;
	string str;
	IsTasksFile= false;

	while ( ( c = getopt_long( argc, argv, "c:m:w:h", longopts, NULL ) ) != -1 ) {
		switch ( c ) {
			case 'c':
				dcapi_config_file_name = optarg;
				break;
			case 'm':
				master_config_file_name = optarg;
				break;
			case 'w':
				wu_id_file_name = optarg;
				break;
			case 'h':
				print_help( argv[0] );
				break;
			default:
				exit( 1 );
		}
	}

	if (optind != argc) {
		fprintf(stderr, "Extra arguments on the command line\n");
		exit( 1 );
	}

	// Specifying the config file is mandatory, otherwise the master can't
	// run as a BOINC daemon
	if ( !dcapi_config_file_name ) {
		fprintf( stderr, "You must specify the dcapi config file\n" );
		exit( 1 );
	}

	if ( !master_config_file_name ) {
		fprintf( stderr, "You must specify the master config file file\n" );
		exit( 1 );
	}

	if ( wu_id_file_name ) {
		IsTasksFile = true;
		cout << "wu_id_file " << wu_id_file_name << endl; 
		cout << "IsTasksFile " << IsTasksFile << endl;
		ifstream infile( wu_id_file_name );
		while ( getline( infile, str ) )
			wu_id_vec.push_back( atoi(str.c_str()) );
		cout << "wu_id_vec" << endl;
		for ( unsigned i=0; i < wu_id_vec.size(); ++i )
			cout << wu_id_vec[i] << endl;
		infile.close();
	}
	
	cout << "config_file_name "        << dcapi_config_file_name  << endl;
	cout << "master_config_file_name " << master_config_file_name << endl;
	
	// Initialize the DC-API
	if ( DC_initMaster( dcapi_config_file_name ) ) {
		fprintf(stderr, "Master: DC_initMaster failed, exiting.\n");
		exit(1);
	}

	// We need the result callback function only
	//DC_setMasterCb( process_result, NULL, NULL );

	do_work( wu_id_vec );
	double total_time = cpuTime() - start_time;
	cout << "total time " << total_time << endl;
		
	return 0;
}

static void print_help(const char *prog)
{
	const char *p;

	// Strip the path component if present
	p = strrchr(prog, '/');
	if (p)
		prog = p;

	printf("Usage: %s {-c|--config} <dcapi_config file> (-m|--master_config) <master_config_file>\n", prog);
	printf("Available options:\n");
	printf("\t-c <file>, --config <file>\tUse the specified dcapi config. file\n");
	printf("\t-h, --help\t\t\tThis help text\n");
	exit(0);
}

void ParseConfigFile( config_params &config_p, string &cnf_head, stringstream &config_sstream )
{
	fstream master_config_file;
	string input_str, str1, str2, str3;
	stringstream sstream;

	string master_config_file_name_str = master_config_file_name;
	master_config_file.open( master_config_file_name_str.c_str() );
	
	cout << "In ParseConfigFile() master_config_file_name " << master_config_file_name_str << endl;
	if ( !master_config_file.is_open() ) {
		cerr << "file " << master_config_file_name_str << " doesn't exist" << endl;
		exit(1);
	}
	cout << endl << "input file opened" << endl;
	
	// cnf_name
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.cnf_name;
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
	// N
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.N;
	sstream.str(""); sstream.clear();
	// rows_count
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.rows_count;
	sstream.str(""); sstream.clear();
	// K
	getline( master_config_file, input_str );
	config_sstream << input_str << endl;
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.K;
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
	// skip_values
	getline( master_config_file, input_str );
	sstream << input_str;
	sstream >> str1 >> str2 >> config_p.skip_values;
	sstream.str(""); sstream.clear();
	// known_values
	getline( master_config_file, input_str );
	sstream << input_str;
	int int_val;
	vector<int> vec_int;
	if ( input_str.find("known_values") != string::npos ) {
		cout << "known_values_vec" << endl;
		while( getline( master_config_file, input_str ) ) {
			sstream.str(""); sstream.clear();
			sstream << input_str;
			while ( sstream >> int_val )
				vec_int.push_back( int_val );
			for ( unsigned i=0; i < vec_int.size(); ++i )
				cout << vec_int[i] << " ";
			cout << endl;
			config_p.known_values_vec.push_back( vec_int );
			vec_int.clear();
			sstream.str(""); sstream.clear();
		}
	}
	
	sstream.str(""); sstream.clear();

	// make cnf_head
	sstream << "p cnf " << config_p.cnf_variables << " " << config_p.cnf_clauses;
	cnf_head = sstream.str();
	sstream.str(""); sstream.clear();

	cout << "cnf_name "          << config_p.cnf_name << endl;
	cout << "cnf_variables "     << config_p.cnf_variables << endl;
	cout << "cnf_clauses "       << config_p.cnf_clauses << endl;
	cout << "N "                 << config_p.N << endl;
	cout << "rows_count "        << config_p.rows_count << endl;
	cout << "K "                 << config_p.K << endl;
	cout << "unsent_needed_wus " << config_p.unsent_needed_wus << endl;
	cout << "total_wus "         << config_p.total_wus << endl;
	cout << "created_wus "       << config_p.created_wus << endl;
	cout << "skip_values "       << config_p.skip_values << endl;
	cout << "*** cnf_head "      << cnf_head << endl;
	cout << endl;

	master_config_file.close();
}

bool do_work( vector<int> &wu_id_vec )
{
	double start_time = cpuTime();

	DC_log( LOG_NOTICE, "Master: Creating work units" );

	processed_wus     = 0;
	all_processed_wus = 0;
	int unsent_count = 0;
	config_params config_p;
	stringstream config_sstream;
	string cnf_head;
	int wus_for_creation_count = 0;
	
	ParseConfigFile( config_p, cnf_head, config_sstream );
	latin_square ls;
	ls.N = config_p.N;
	ls.K = config_p.K;
	ls.rows_count = config_p.rows_count;
	ls.skip_values = config_p.skip_values;
	if ( config_p.cnf_name.find( "diag" ) != string::npos )
		ls.problem_type = "diag";
	if ( config_p.cnf_name.find( "inc0" ) != string::npos )
		ls.problem_type = "inc0";
	else if ( config_p.cnf_name.find( "inc70" ) != string::npos )
		ls.problem_type = "inc70";
	cout << "ls.N "            << ls.N << endl;
	cout << "ls.rows_count "   << ls.rows_count << endl;
	cout << "ls.K "            << ls.K << endl;
	cout << "ls.skip_values "  << ls.skip_values << endl;
	cout << "ls.problem_type " << ls.problem_type << endl;
	bool IsLastGenerating = false;

	if ( IsTasksFile ) // get PROBLEMS_IN_WU assumptions for every task
		create_wus( ls, config_sstream, config_p, cnf_head, wus_for_creation_count, wu_id_vec, IsLastGenerating );
	else
	{
		// create wus - supply needed level of unsent wus
		wus_for_creation_count = 0;
		int old_wus_for_creation_count = 0;
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

			GetCountOfUnsentWUs( unsent_count );
			cout << "unsent_count " << unsent_count << endl;
		
			if ( unsent_count < 0 ) {
				cout << "SQL error. Skipping creating wus" << endl;
				sleep( 60 );
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
				create_wus( ls, config_sstream, config_p, cnf_head, wus_for_creation_count, wu_id_vec, IsLastGenerating );
			}
			else {
				if ( old_wus_for_creation_count != wus_for_creation_count )
					cout << "wus_for_creation_count " << wus_for_creation_count << endl;
				old_wus_for_creation_count = wus_for_creation_count;
			}
			
			if ( !IsLastGenerating )
				sleep( 7200 ); // wait 2 hours
		}
	}

	cout << "wus_for_creation_count " << wus_for_creation_count << endl;
	DC_log( LOG_NOTICE, "\n Work finished" );

	double total_time = cpuTime() - start_time;
	cout << "total time " << total_time << endl;

	return 0;
}

void create_wus( latin_square &ls, stringstream &config_sstream, config_params &config_p, 
	             string cnf_head, int wus_for_creation_count, vector<int> &wu_id_vec, bool &IsLastGenerating )
{
	int index = 0;
	DC_Workunit *wu;
	ofstream output;
	string wu_tag_str;
	stringstream sstream;

	cout << "Start create_wus()" << endl;
	cout << "cnf_head " << cnf_head << endl;
	cout << "cnf_name " << config_p.cnf_name << endl;

	if ( IsTasksFile ) {
		cout << "IsTasksFile mode" << endl;
		ls.max_values_len = wu_id_vec[ wu_id_vec.size()-1 ] * PROBLEMS_IN_WU;
		wus_for_creation_count = wu_id_vec.size();
		config_p.created_wus = 0;
	}
	else
		ls.max_values_len = wus_for_creation_count * PROBLEMS_IN_WU;
	
	cout << "wus_for_creation_count " << wus_for_creation_count << endl;
	cout << "max_values_len " << ls.max_values_len << endl;
	ls.verbosity = 0;

	if ( config_p.known_values_vec.size() > 0 )
		ls.known_values_vec = config_p.known_values_vec;

	ls.Make_Permutations_n_k( );
	cout << "Make_Permutations_n_k() done" << endl;

	if ( IsTasksFile ) {
		int assumption_index;
		vector< vector<int> > tmp_positive_literals;
		for ( int i = 0; i < wu_id_vec.size(); i++ )
			for ( int j = 0; j < PROBLEMS_IN_WU; j++ ) {
				assumption_index = (wu_id_vec[i]-1) * PROBLEMS_IN_WU + j;
				tmp_positive_literals.push_back( ls.positive_literals[assumption_index] );
			}
		ls.positive_literals = tmp_positive_literals;
		cout << "new positive_literals.size() " << ls.positive_literals.size() << endl;
	}
	
	ofstream tmp_out;
	ls.WriteCurrentState( tmp_out, 0 );
	cout << "Got latin values" << endl;
	vector<int> :: iterator vec_it;
	int values_index = 0;
	int wu_index = 0;
	bool IsAddingWUneeded;
	bool IsFastExit = false;
	unsigned new_created_wus = 0;
	for( index = config_p.created_wus; index < config_p.created_wus + wus_for_creation_count; index++ ) {
		if ( IsFastExit )
			break;
		IsAddingWUneeded = false; // if no values will be added then WU not needed
		//sstream << "c " << config_p.cnf_name << endl; // TODO return it after changing client app
		sstream << "c ls_10_3_inc70" << endl;
		for ( int i = 0; i < PROBLEMS_IN_WU; i++ ) {
			sstream << "c ";
			if ( values_index >= ls.positive_literals.size() ) {
				cout << "in create_wus() last ls.positive_literals was added to WU" << endl;
				cout << "values_index " << values_index << endl;
				IsFastExit = true; // add last values to WU and exit
				IsLastGenerating = true; // tell to high-level function about ending of generation
				break;
			}
			for ( vec_it = ls.positive_literals[values_index].begin(); vec_it != ls.positive_literals[values_index].end() - 1; vec_it++ )
				sstream << *vec_it << " ";
			sstream << *vec_it << endl;
			values_index++;
			IsAddingWUneeded = true;
		}

		if ( !IsAddingWUneeded ) {
			cout << "IsAddingWUneeded true" << endl;
			break; // don't create new WU
		}
		sstream << cnf_head;
		output.open( "wu-input.txt", ios_base :: out );
		if ( !output.is_open() ) {
			DC_log(LOG_ERR, "Failed to create wu-input.txt: %s",
			strerror(errno));
			exit(1);
		}
		output << sstream.rdbuf();
		sstream.str( "" ); sstream.clear();
		output.close( );
		
		sstream << config_p.cnf_name;
		if ( IsTasksFile ) // mark name with "additional" label
			sstream << "-add";
		sstream << "--" << index + 1; // save info about CNF name
		wu_tag_str = sstream.str();
		sstream.str( "" ); sstream.clear();
		wu = DC_createWU( "pdsat", NULL, 0, wu_tag_str.c_str() );
		if ( !wu ) {
			DC_log( LOG_ERR, "Work unit creation has failed" );
			exit(1);
		}

		if (DC_addWUInput( wu, INPUT_LABEL, "wu-input.txt",
				           DC_FILE_VOLATILE ) ) {
			DC_log( LOG_ERR, "Failed to register WU input file" );
			exit(1);
		}
		if (DC_addWUOutput( wu, OUTPUT_LABEL) ) {
			DC_log(LOG_ERR, "Failed to register WU output file");
			exit(1);
		}
		if (DC_submitWU(wu)) {
			DC_log(LOG_ERR, "Failed to submit WU");
			exit(1);
		}
		new_created_wus++;
	}
	
	cout << "new_created_wus " << new_created_wus << endl;
	config_p.created_wus += new_created_wus;
	// update ls value
	config_p.skip_values = ls.skip_values = ls.values_checked;
	// update config file
	
	if ( !IsTasksFile ) { // don't update if additional wus from file ewre created
		ofstream master_config_file;
		string master_config_file_name_str = master_config_file_name;
		master_config_file.open( master_config_file_name_str.c_str() );
		master_config_file << config_sstream.str();
		master_config_file << "created_wus = " << config_p.created_wus << endl; // update master config file
		master_config_file << "skip_values = " << config_p.skip_values << endl; 
		if ( config_p.known_values_vec.size() > 0 ) {
			master_config_file << "known_values" << endl;
			for ( unsigned i=0; i < config_p.known_values_vec.size(); ++i ) {
				master_config_file << config_p.known_values_vec[i];
				if ( i != config_p.known_values_vec.size() - 1 )
					master_config_file << " ";
			}
		}
		master_config_file.close();
	}

	cout << "created_wus " << config_p.created_wus << endl;
	cout << "skip_values " << config_p.skip_values << endl;
	cout << "---***---" << endl;
}

// Concatenate all results in their original order to form the final output
// Try to find SAT in results
bool find_sat( int cnf_index )
{
	bool flag = false;
	ofstream output_file;
	stringstream final_sstream;
	stringstream sstream;
	string final_output_name; 
	
	sstream << "final_output" << cnf_index << ".txt";
	final_output_name = sstream.str( );
	sstream.str( "" );
	sstream.clear();
	output_file.open( final_output_name.c_str() );
	if ( !output_file.is_open() ) {
		DC_log(LOG_ERR, "Cannot open final_output.txt for writing: %s",
			strerror(errno));
		exit(1);
	}

	string result_filename, input_str;
	ifstream result_file;
	/*for( int i = 1; i < created_wus + 1; i++ ) { // try all files
		sstream << "result_" << i << ".txt";
		result_filename = sstream.str();
		sstream.str( "" );
		sstream.clear();
		result_file.open( result_filename.c_str(), ios_base :: in );
		
		if ( result_file.is_open( ) ) {     
			//DC_log(LOG_NOTICE, "The result for WU %d exists ", i);
			getline( result_file, input_str );
			final_sstream << input_str;
			final_sstream << " result_" << i << endl;
			int sat_pos = input_str.find( "SAT" );
			if ( !sat_pos ) { // if SAT then
				flag = true;
				while ( getline( result_file, input_str ) )
					final_sstream << input_str << endl;
			}
			result_file.close();
			result_file.clear();
		}
	}*/
    
	output_file << final_sstream.rdbuf();
	output_file.close();
	
	return flag;
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
		DC_log(LOG_ERR, "Cannot open final_output.txt for writing: %s",
			strerror(errno));
		exit(1);
	}

	string input_str;
	ifstream result_file;
	result_file.open( output_filename.c_str(), ios_base :: in );
	if ( !result_file.is_open( ) ) {
		DC_log(LOG_ERR, "Cannot open result_file: %s",
			strerror(errno));
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

void GetCountOfUnsentWUs( int &unsent_count )
{
	char *host = "localhost";
    char *db = "boinc_pdsat";
	char *user;
    char *pass;
	MYSQL *conn;
	
	ifstream pass_file("pass");
	if ( !pass_file.is_open() ) {
		cerr << "psswd_file not open" << endl;
		exit(1)
	}
	string str;
	getline( pass_file, str );
	user = str.c_str();
	getline( pass_file, str );
	pass = str.c_str();
	cout << "user " << user << endl;
	cout << "pass " << pass << endl;
	
	conn = mysql_init(NULL);
	if(conn == NULL)
		cerr << "Error: can't create MySQL-descriptor\n";

	if(!mysql_real_connect(conn, host, user, pass, db, 0, NULL, 0))
		cerr << "Error: can't connect to MySQL server\n";

	vector< vector<stringstream *> > result_vec;
	string str = "SELECT COUNT(*) FROM workunit WHERE id IN(SELECT workunitid FROM result WHERE server_state = 2)";
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