#include <dc_client.h>
/* Don't forget to put the DC-API libraries 
/* (libdc-client-boinc.lib and 
/* libdc-client-boinc-debug.lib) in the lib dir.
/* These files are included in the DC-API client 
/* libraries for windows package */

#ifndef _WIN32
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#else
#include <dc_win32.h>
#include <windows.h>
#endif
//#include <unistd.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <ctype.h>
#include <fstream>
#include <string>
#include <sstream>

#include "mpi_base.h"
#include "common.h"

using namespace std;

#pragma warning( disable : 4996 )

/*#define INPUT_FILENAME	"in.txt"
#define OUTPUT_FILENAME	"out.txt"*/

ifstream infile;
fstream outfile;
ifstream chptfile;
string infile_name;
string outfile_name;
string ckptfile_name;
unsigned long long int last_iteration_done;

/************************************/
/* Application Specific Prototypes */
/**********************************/

/**************************/
/* Additional prototypes */
/************************/

static void prepare_toexit(void);
static void init_files(void);
static void do_work(void);
/**********************/

/************************************/
/* Application Specific Functions  */
/**********************************/

void prepare_toexit (void )
{
	string first_string = "";
	stringstream sstream;

	cout << "Preparing to exit..." << endl;

	if ( infile.is_open() ) infile.close();
	if ( chptfile.is_open() ) chptfile.close();

	/* If output file is empty, we have to print something into it,
	* becasue BOINC server does not create 0 length file on the data server
	*/
	if ( outfile.is_open() ) 
	{
		outfile.seekg(0, ios::end); // put the "cursor" at the end of the file
		int length = outfile.tellg(); // find the position of the cursor
		if ( !length )
		{
			sstream << "Program terminated without output \n";
			sstream << "APP: Output is empty, placing msg in out.txt";
			outfile << sstream.rdbuf();
		}
		outfile.close();
	}
}

void init_files(void)
{
	char *resolved_name;
  
	/* open input file */ 
	resolved_name = DC_resolveFileName(DC_FILE_IN, INPUT_LABEL);
    if (!resolved_name) {
        fprintf(stderr, "APP: Could not resolve input file\n");
        prepare_toexit();
        DC_finishClient(1);
	}
	else
	{
	    infile_name = resolved_name;
	    free( resolved_name );
	}
	infile.open( infile_name.c_str( ), ios_base :: in );
	
	if ( !infile.is_open( ) )
	{
		cerr << "APP: Input file " << infile_name << " does not exists";
		prepare_toexit();
		DC_finishClient(1);
	}
	
	// open output file
	resolved_name = DC_resolveFileName(DC_FILE_OUT, OUTPUT_LABEL);
	if (!resolved_name) 
	{
		cerr << "APP: Could not resolve output file" << endl;
		prepare_toexit();
		DC_finishClient(1);
	}
	else
	{
	    outfile_name = resolved_name;
	    free( resolved_name );
	}
	outfile.open( outfile_name.c_str( ), ios_base :: out );
	if ( !outfile.is_open( ) ) {
		cerr << "APP: Failed to open output file" << endl;
		prepare_toexit();
		DC_finishClient(1);
	}

	// Check the checkpoint file
	last_iteration_done = 0;
	resolved_name = DC_resolveFileName(DC_FILE_IN, DC_CHECKPOINT_FILE);
	if ( resolved_name )
	{
		ckptfile_name = resolved_name;
		free( resolved_name );
		chptfile.open( ckptfile_name.c_str( ), ios_base :: in );
		string line_str;
		stringstream sstream;
		// ckpt file exists: read and set everything according to it
		getline( chptfile, line_str );
		sstream << line_str;
		sstream >> last_iteration_done;
		//sstream << " last_iteration_done " << endl;
		//outfile << sstream.rdbuf();
	}
}

void do_work(void)
{
	stringstream sstream;
	string input_str;

	char *input_cnf_name;
	int IsSATFinded;
	int current_obj_val = 0;
	double cnf_time_from_node = 0;
	MPI_Base mpi_base;
	unsigned int value[FULL_MASK_LEN];
	for ( int i = 0; i < FULL_MASK_LEN; i++ )
		mpi_base.full_mask[i] = 0;
	for ( int i = 0; i < FULL_MASK_LEN; i++ )
		mpi_base.part_mask[i] = 0;
	for ( int i = 0; i < FULL_MASK_LEN; i++ )
		value[i] = 0;

	string str1, str2;
	getline( infile, input_str );
	sstream << input_str;
	sstream >> str1 >> str2;
	if ( str2 == "full_mask" )
		sstream >> mpi_base.full_mask[0] >> mpi_base.full_mask[1] >> mpi_base.full_mask[2] >> mpi_base.full_mask[3];

	sstream.str( "" ), sstream.clear();
	getline( infile, input_str );
	sstream << input_str;
	sstream >> str1 >> str2;
	if ( str2 == "part_mask" )
		sstream >> mpi_base.part_mask[0] >> mpi_base.part_mask[1] >> mpi_base.part_mask[2] >> mpi_base.part_mask[3];

	sstream.str( "" ), sstream.clear();
	getline( infile, input_str );
	sstream << input_str;
	sstream >> str1 >> str2;
	if ( str2 == "value" )
		sstream >> value[0] >> value[1] >> value[2] >> value[3];

	IsSATFinded = 0;

	input_cnf_name = new char[strlen( infile_name.c_str() ) + 1];
	strcpy( input_cnf_name, infile_name.c_str( ) );

	mpi_base.input_cnf_name = input_cnf_name;
	if ( !mpi_base.ReadIntCNF( ) ) // Read original CNF
		{ cout << "\n Error in ReadIntCNF" << endl; exit(1); }
	//
	int IsHardProblem = 0;
	int sort_type = 0;
	int max_attempt = 0;
	
	cout << "mpi_base.core_len " << mpi_base.core_len << endl;
	if ( !dminisat_solve( input_cnf_name, mpi_base.full_mask, mpi_base.part_mask, value,
						  mpi_base.solver_type, mpi_base.core_len, mpi_base.corevars_activ_type, 
						  &IsSATFinded, &mpi_base.b_SAT_set_array, sort_type,
					      &cnf_time_from_node, 0, IsHardProblem, last_iteration_done,
					      &max_attempt ) )
	{ printf( "\n Error in dminisat_solve" ); }

	cout << endl << "Answer is " << IsSATFinded << endl;
	sstream.str( "" ), sstream.clear();
	if ( IsSATFinded )
	{
		sstream << "SAT" << endl;
		sstream << "full_mask ";
		for ( int i = 0; i < FULL_MASK_LEN; i++ )
			sstream << mpi_base.full_mask[i] << " ";
		sstream << endl;
		for ( int i = 0; i < FULL_MASK_LEN; i++ )
			sstream << mpi_base.part_mask[i] << " ";
		sstream << endl;
		for ( int i = 0; i < FULL_MASK_LEN; i++ )
			sstream << value[i] << " ";
		sstream << endl;
		sstream << "core_len " << mpi_base.core_len << endl;
		for ( int i = 0; i < mpi_base.core_len; i++ )
		{
			if ( mpi_base.b_SAT_set_array[i] ) sstream << i + 1;
				else sstream << -( i + 1 );
			sstream << " ";
		}
	}
	else
		sstream << "UNSAT";
	if ( max_attempt > 0 )
	    sstream << " max_attempt " << max_attempt;
	    
	outfile << sstream.rdbuf();
	cout << sstream.str() << endl;
	cout << "Must be Ok" << endl;
	outfile.close();
	if ( outfile.is_open() ) {
		perror("APP: Closing the output file has failed");
		DC_finishClient(1);
	}
	delete[] input_cnf_name;
}

int main(int ac, char* av [])
{
  int i;
  //DC_ClientEvent *event;
  
  i = DC_initClient();
  if (i) {
      fprintf(stderr, "APP: Failed to initialize DC-API. Ret: %d", i);
      prepare_toexit();
      DC_finishClient(1);
  }

  init_files();

  fprintf(stdout, "APP: Input processed...notifying core client\n");
  do_work();
  prepare_toexit();
  DC_finishClient(0);
  return 0;
}

#ifdef _WIN32

extern int parse_command_line( char *, char ** );

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode)
{
    LPSTR command_line;
    char* argv[100];
    int argc;
    
    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif
