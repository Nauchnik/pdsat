
// version 11/01/2012
// [+] parameter --results+to_skip= added for non creating new wus for some forst results
// [*] adding new WU for every processed result, not only from current launch
// [+] parameter --skip_wus= added
// saving of WU_id added
// FIRST_WU_CREATION changed to 65536
// added procesing all CNFs in folder 
// for every CNF - own final_output file
// 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dc.h>

#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <string>
#include <fstream>
#include <sstream>
#include <dirent.h>

#include "mpi_base.h"
#include "common.h"

using namespace std;

#define UNSENT_WU_COUNT 2048
#define FIRST_WU_CREATION 4096

// Number of WUs we have created
static int created_wus;
static int new_created_wus;
// Number of results we have received so far
static int all_processed_wus;
static int processed_wus;
static int unsent_wus;
static int running_wus;
static int skip_wus;
static int results_to_skip;
unsigned long long max_wus;

DC_Workunit **wu_arr;
DC_Workunit **unsent_wu_arr;
int *wu_status_arr; // 0 - unsent, 1 - running, 2 - finished

// Command line options

static inline double cpuTime( void ) 
{
    return ( double )clock( ) / CLOCKS_PER_SEC; 
}

static const struct option longopts[] =
{
	{ "config",	         required_argument, NULL, 'c' },
	{ "help",	         no_argument,		 NULL, 'h' },
	{ "skip_wus",        optional_argument, NULL, 's' },
	{ "results_to_skip", optional_argument, NULL, 'r' },
	{ NULL }
};

int getdir( string dir, vector<string> &files )
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << endl << "Error in opening " << dir;
        return 1;
    }
    while ((dirp = readdir(dp)) != NULL) 
	{ files.push_back(string(dirp->d_name)); }
    closedir(dp);
    return 0;
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
	if ( !output_file.is_open() )
	{
		DC_log(LOG_ERR, "Cannot open final_output.txt for writing: %s",
			strerror(errno));
		exit(1);
	}

	string input_str;
	ifstream result_file;
	result_file.open( output_filename.c_str(), ios_base :: in );
	if ( !result_file.is_open( ) )
	{
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

// Copies the output files into the master working directory
static void process_result(DC_Workunit *wu, DC_Result *result)
{
	char *output_filename, *tag, *id, cmd[256];

	// Extract our own private identifier
	if (!result)
	{
		DC_log(LOG_WARNING, "Work unit %s has failed", tag);
		// In a real application we would create a new work unit
		// instead of the failed one and submit it again. Here we
		// just let it fail.
		return;
	}

	tag = DC_getWUTag(wu);
	id  = DC_getWUId(wu);

	output_filename = DC_getResultOutput( result, OUTPUT_LABEL );
	if ( !output_filename )
	{
		DC_log(LOG_WARNING, "Work unit %s contains no output file", 
			tag);
		free(tag);
		free(id);
		return;
	}

	// add it anyway, for different CNFs tags different
	add_result_to_file( output_filename, tag, id );
	
	DC_log(LOG_NOTICE, "Work unit %s with id %s has completed", tag, id);
		
	// We no longer need the work unit
	DC_destroyWU(wu);
	free(output_filename);
	free(tag);
	free(id);
}

static void print_help(const char *prog) __attribute__((noreturn));
static void print_help(const char *prog)
{
	const char *p;

	// Strip the path component if present
	p = strrchr(prog, '/');
	if (p)
		prog = p;

	printf("Usage: %s {-c|--config} <config file>\n", prog);
	printf("Available options:\n");
	printf("\t-c <file>, --config <file>\tUse the specified config. file\n");
	printf("\t-h, --help\t\t\tThis help text\n");
	printf("\t--skip_wus\t\t\tCount of wus to skip\n");
	printf("\t--results_to_skip\t\t\tCount of results with skiping creating new wus\n");
	exit(0);
}

int main( int argc, char *argv[] )
{
	char *config_file = NULL;
	int c;

	double start_time = cpuTime();

	while ( ( c = getopt_long( argc, argv, "c:h", longopts, NULL ) ) != -1 )
	{
	    switch ( c )
	    {
		case 'c':
		config_file = optarg;
		break;
	    case 'h':
		print_help( argv[0] );
		break;
	    default:
	    exit( 1 );
	}
																		    }
	// Specifying the config file is mandatory, otherwise the master can't
	// run as a BOINC daemon
	if ( !config_file ) {
		fprintf( stderr, "You must specify the config file\n" );
		exit( 1 );
	}

	// Initialize the DC-API
	if ( DC_initMaster( config_file ) ) {
		fprintf(stderr, "Master: DC_initMaster failed, exiting.\n");
		exit(1);
	}

	// We need the result callback function only
	DC_setMasterCb( process_result, NULL, NULL );
	
	string cmd;
	cmd = "rm -r *final*"; // clear old  final files
	system( cmd.c_str() );
	cmd = "rm -r *result*";
	cout << endl << "Deleting old results and final files in folder";
	system( cmd.c_str( ) );

	bool flag = true;
	while ( flag )
	    DC_processMasterEvents( 300 ); // wait new results 5 min
	
	double total_time = cpuTime() - start_time;
	printf("Total time:%-12.2f s \n", total_time);
	//cout << "total time " << total_time;
		
	return 0;
}
