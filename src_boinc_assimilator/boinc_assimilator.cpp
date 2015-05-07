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

#ifdef _WIN32
#include <my_global.h> // Include this file first to avoid problems
#endif
#include <mysql.h>

void ProcessQuery(MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec);
void MakeHTMLfromWU(MYSQL *conn, string wu_id_str);

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

	// read sat from file
	std::ifstream sat_file("sat_output");
	std::string wu_name;
	while (getline(sat_file, str)) {
		sstream << str;
		sstream >> wu_name;
		sstream << "SELECT userid, teamid, mod_time FROM result WHERE validate_state = 1 AND id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		cout << str << endl;
		ProcessQuery(conn, str, result_vec);
	}
	sat_file.close();
}

void MakeHTMLfromWU(MYSQL *conn, string wu_id_str)
{
	cout << "wu_id_str " << wu_id_str << endl;
	vector< vector<stringstream *> > result_vec;
	string str = "SELECT id FROM result WHERE workunitid=" + wu_id_str;
	cout << str << endl;

	ProcessQuery(conn, str, result_vec);

	if (result_vec.size() == 0) {
		cout << "result_vec.size() == 0" << endl;
		return;
	}

	stringstream sstream;
	vector<int> resultid_vec;
	vector<int> userid_vec;
	vector<string> username_vec;
	vector<int> teamid_vec;
	vector<string> teamname_vec;
	vector<string> mod_time_vec;
	unsigned u_val;
	for (unsigned i = 0; i < result_vec.size(); i++)
		for (unsigned j = 0; j < result_vec[i].size(); j++) {
			*result_vec[i][j] >> u_val;
			resultid_vec.push_back(u_val);
			delete result_vec[i][j];
		}
	result_vec.clear();

	for (vector<int>::iterator it = resultid_vec.begin(); it != resultid_vec.end(); it++) {
		sstream << "SELECT userid, teamid, mod_time FROM result WHERE validate_state = 1 AND id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		cout << str << endl;
		ProcessQuery(conn, str, result_vec);
		cout << "workunitid " << wu_id_str << endl;
		cout << "resultid " << *it << endl;
		for (unsigned i = 0; i < result_vec.size(); i++) {
			*result_vec[i][0] >> u_val; // get userid
			userid_vec.push_back(u_val);
			cout << "userid " << u_val << endl;
			*result_vec[i][1] >> u_val; // get teamid
			teamid_vec.push_back(u_val);
			cout << "teamid " << u_val << endl;
			str = (*result_vec[i][2]).str(); // get mod_time
			cout << str << endl;
			cout << "mod_time " << str << endl;
			cout << endl;
			mod_time_vec.push_back(str);
			for (unsigned j = 0; j < result_vec[i].size(); j++)
				delete result_vec[i][j];
		}
		result_vec.clear();
	}

	cout << "teamid_vec:" << endl;
	for (unsigned i = 0; i < teamid_vec.size(); i++)
		cout << teamid_vec[i] << endl;

	// get names of users
	for (vector<int>::iterator it = userid_vec.begin(); it != userid_vec.end(); it++) {
		sstream << "SELECT name FROM user WHERE id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		cout << str << endl;
		ProcessQuery(conn, str, result_vec);

		for (unsigned i = 0; i < result_vec.size(); i++) {
			str = (*result_vec[i][0]).str();
			cout << str << endl;
			username_vec.push_back(str);
			delete result_vec[i][0];
		}
		result_vec.clear();
	}

	// get names of teams
	for (vector<int>::iterator it = teamid_vec.begin(); it != teamid_vec.end(); it++) {
		if (*it == 0)
			teamname_vec.push_back("");

		sstream << "SELECT name FROM team WHERE id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		cout << str << endl;
		ProcessQuery(conn, str, result_vec);

		for (unsigned i = 0; i < result_vec.size(); i++) {
			str = (*result_vec[i][0]).str();
			cout << str << endl;
			teamname_vec.push_back(str);
			delete result_vec[i][0];
		}
		result_vec.clear();
	}

	sstream << "<tr>" << endl << "<td> <b>" << mod_time_vec[0] << " UTC </b> </td>" << endl;
	sstream << "<td> <a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[0] <<
		"'>" << username_vec[0] << "</a>";
	if (teamname_vec[0] != "")
		sstream << " from " << teamname_vec[0];
	sstream << " /" << endl;
	sstream << "<a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[1] <<
		"'>" << username_vec[1] << "</a>";
	if (teamname_vec[1] != "")
		sstream << " from " << teamname_vec[1];
	sstream << " </td>" << endl;
	sstream << "<td>Bivium</td>" << endl << "<td> </td>" << endl << "</tr>" << endl;

	cout << sstream.rdbuf();
}

void ProcessQuery(MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec)
{
	// Дескриптор результирующей таблицы
	MYSQL_RES *res;
	// Дескриптор строки
	MYSQL_ROW row;
	int num_fields;

	if (mysql_query(conn, str.c_str()) != 0)
		cerr << "Error: can't execute SQL-query\n";

	// Получаем дескриптор результирующей таблицы
	res = mysql_store_result(conn);

	if (res == NULL)
		cerr << "Error: can't get the result description\n";

	num_fields = mysql_num_fields(res);
	stringstream *sstream_p;
	vector<stringstream *> result_data;

	// Если имеется хотя бы одна запись - выводим список каталогов
	if (mysql_num_rows(res) > 0) {
		// В цикле перебираем все записи результирующей таблицы
		while ((row = mysql_fetch_row(res)) != NULL) {
			for (int i = 0; i < num_fields; i++) {
				//cout << row[i] << endl;
				sstream_p = new stringstream();
				*sstream_p << row[i]; // get value
				result_data.push_back(sstream_p);
			}
			result_vec.push_back(result_data);
			result_data.clear();
		}
	}

	// Освобождаем память, занятую результирующей таблицей
	mysql_free_result(res);
}