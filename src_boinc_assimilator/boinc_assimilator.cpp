#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#ifdef _WIN32
#include "dirent.h"
#else
#include <dirent.h>
#endif

int getdir( std::string dir, std::vector<std::string> &files )
{
    DIR *dp;
	std::string cur_name;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        std::cout << std::endl << "Error in opening " << dir;
        return 1;
    }
    while ((dirp = readdir(dp)) != NULL) { 
		cur_name = std::string(dirp->d_name);
		if ( cur_name[0] != '.' ) files.push_back(cur_name); 
	}
    closedir(dp);
    return 0;
}

int main( int argc, char *argv[] )
{
#ifdef _DEBUG
	argc = 1;
	argv[0] = "./assimilator";
#endif
	std::vector<std::string> file_names = std::vector<std::string>();
	getdir( ".", file_names );
	std::cout << "file_names.size() " << file_names.size() << std::endl;
	std::string system_str;
	std::ifstream ifile;
	std::string str, program_name = argv[0];
	program_name.erase( std::remove(program_name.begin(), program_name.end(), '.'), program_name.end() );
	program_name.erase( std::remove(program_name.begin(), program_name.end(), '/'), program_name.end() );
	std::cout << "program_name " << program_name << std::endl;
	std::stringstream sstream, final_sstream, sat_sstream;
	bool isCorrectFile;
	
	for ( unsigned i=0; i < file_names.size(); i++ ) {
		if ( ( file_names[i].find( "assimilator" ) != std::string::npos ) || 
			 ( file_names[i].find( "output" ) != std::string::npos ) ||
			 ( file_names[i].find( "out" ) != std::string::npos ) 
			 )
			continue;
		ifile.open(file_names[i].c_str());
		isCorrectFile = false;
		while ( getline( ifile, str ) ) {
			// check if file contain result SAT@home info
			if ( ( str.find( "SAT" ) != std::string::npos ) && ( str.find( "INTERRUPTED" ) != std::string::npos ) )
				isCorrectFile = true;
			if ( str.find( " SAT" ) != std::string::npos ) {
				std::cout << "SAT found" << std::endl;
				std::ofstream sat_out("sat_output", std::ios_base::app );
				sat_out << file_names[i] << " " << str << std::endl;
				sat_out.close(); sat_out.clear();
			}
			sstream << str;
		}
		ifile.close(); ifile.clear();
		if ( isCorrectFile ) {
			final_sstream << file_names[i] << " " << sstream.str() << std::endl;
			system_str = "rm ";
			system_str += file_names[i];
			std::cout << "system_str " << system_str << std::endl;
			system( system_str.c_str() );
		}
		sstream.clear(); sstream.str("");
	}
	
	std::ofstream ofile( "final_output", std::ios_base::app );
	ofile << final_sstream.rdbuf();
	ofile.close();

	std::cout << "*** done" << std::endl;
}