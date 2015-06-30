#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Mols.h"

#ifdef _WIN32
#include "dirent.h"
#else
#include <dirent.h>
#endif

#ifdef _WIN32
#include <my_global.h> // Include this file first to avoid problems
#endif

#include <mysql.h>
				
void getMOLS(std::vector<MOLS> &pairs_vec);
void touchBoincResultFiles();
void ProcessQuery(MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec);
void MakeHTMLfromWU(MYSQL *conn, std::string wu_id_str, MOLS pair_MOLS);
int getdir(std::string dir, std::vector<std::string> &files);

int main( int argc, char *argv[] )
{
#ifdef _DEBUG
	argc = 1;
	argv[0] = "./assimilator";
#endif
	if (argc < 2) {
		std::cerr << "program [DB password]" << std::endl;
		return 1;
	}
	
	//connection params
	char *host = "localhost";
	char *db = "boinc_pdsat";
	char *user = "boinc_pdsat";
	char *pass = argv[1];

	touchBoincResultFiles(); // touch files from errors files
									
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
		while (getline(ifile, str)) {
			// check if file contain result SAT@home info
			if ((str.find("UNSAT") != std::string::npos) || (str.find("INTERRUPTED") != std::string::npos))
				isCorrectFile = true;
			if (str.find(" SAT") != std::string::npos) {
				std::cout << "SAT found" << std::endl;
				std::ofstream sat_out("sat_output", std::ios_base::app);
				sat_out << file_names[i] << " " << str << std::endl;
				sat_out.close(); sat_out.clear();
			}
			sstream << str;
		}
		ifile.close(); ifile.clear();
		if (isCorrectFile) {
			final_sstream << file_names[i] << " " << sstream.str() << std::endl;
			system_str = "rm ";
			system_str += file_names[i];
			std::cout << "system_str " << system_str << std::endl;
			system(system_str.c_str());
		}
		sstream.clear(); sstream.str("");
	}
	
	std::ofstream ofile( "final_output", std::ios_base::app );
	ofile << final_sstream.rdbuf();
	ofile.close();
	
	// read sat from file
	std::ifstream sat_file("sat_output");
	MYSQL *conn;
	conn = mysql_init(NULL);
	if (conn == NULL)
		std::cerr << "Error: can't create MySQL-descriptor" << std::endl;
	
	// Устанавливаем соединение с базой данных
	if (!mysql_real_connect(conn, host, user, pass, db, 0, NULL, 0))
		std::cerr << "Error: can't connect to MySQL server" << std::endl;
	
	// Устанавливаем кодировку соединения, чтобы предотвратить
	// искажения русского текста
	//if(mysql_query(conn, "SET NAMES 'utf8'") != 0)
	//   cerr << "Error: can't set character set\n";

	std::vector<MOLS> pairs_vec;
	getMOLS(pairs_vec);
	// clean file
	std::ofstream MOLS_out_file("MOLS_out", std::ios_base::out);
	MOLS_out_file.close();
	
	std::string wu_part_name;
	std::vector< std::vector<std::stringstream *> > result_vec;
	unsigned MOLS_pair_index = 0;
	while (getline(sat_file, str)) {
		sstream << str;
		sstream >> wu_part_name;
		sstream.str(""); sstream.clear();
		if (MOLS_pair_index > pairs_vec.size()) {
			std::cerr << "MOLS_pair_index > pairs_vec.size()" << std::endl;
			std::cerr << MOLS_pair_index << " " << pairs_vec.size() << std::endl;
			exit(1);
		}
		MakeHTMLfromWU(conn, wu_part_name, pairs_vec[MOLS_pair_index++]);
	}
	sat_file.close();
	
	std::cout << "*** done" << std::endl;
}

void getMOLS(std::vector<MOLS> &pairs_vec)
{
	ifstream ifile("sat_output");
	string str;
	int n = 10, r = 2;

	unsigned str_count = 0;
	while (getline(ifile, str)) {
		str_count++;
		cout << str << endl;
		stringstream sstream;
		sstream << str;
		string sat_assign_str = "";
		while ((sat_assign_str.length() < 100) && (!sstream.eof()))
			sstream >> sat_assign_str;

		unsigned ones_count = 0;
		for (unsigned i = 0; i < sat_assign_str.size(); i++) {
			if (sat_assign_str[i] == '1')
				ones_count++;
		}

		if (ones_count != sat_assign_str.size() / 10) {
			std::cerr << "ones_count != sat_assign_str.size() / 2" << std::endl;
			std::cerr << "ones_count " << ones_count << std::endl;
			std::cerr << "sat_assign_str.size() / 2 " << sat_assign_str.size() / 2 << std::endl;
			continue;
		}

		if (sat_assign_str.size() != r*pow(n, 3)) {
			std::cerr << "sat_assign_str.size() != r*pow(n, 3)" << std::endl;
			std::cerr << "sat_assign_str.size() " << sat_assign_str.size() << std::endl;
			std::cerr << "r*pow(n, 3) " << r*pow(n, 3) << std::endl;
			continue;
		}

		std::cout << "new pair" << std::endl;
		MOLS mols(sat_assign_str, 10, 2, false);
		//mols.printToCout();
		if (mols.Squares[0].check(true) && mols.Squares[1].check(true) && mols.ortogonalitycheck())
			pairs_vec.push_back(mols);
	}
	ifile.close();
	std::cout << "str_count " << str_count << std::endl;
	std::cout << "pairs_vec.size() " << pairs_vec.size() << std::endl;

	std::ofstream mols_file("MOLS");
	for (unsigned i = 0; i < pairs_vec.size(); i++) {
		pairs_vec[i].Squares[0].normalizefirstrow();
		pairs_vec[i].Squares[1].normalizefirstrow();
		mols_file << pairs_vec[i].HtmlstringView() << std::endl;
	}
	mols_file.close();
}

void touchBoincResultFiles()
{
	std::string error_file_name = "errors";
	std::ifstream error_file(error_file_name.c_str());
	if (!error_file.is_open()) {
		std::cerr << "can't open file with name errors";
		exit(1);
	}
	
	std::vector< std::string > files_names_vec;
	std::string str, prefix_str, launch_str;
	prefix_str = "couldn't copy file ";
	while (getline(error_file, str)) {
		if (str.find(prefix_str) != std::string::npos) {
			str = str.substr(prefix_str.size(), str.size() - prefix_str.size());
			if (std::find(files_names_vec.begin(), files_names_vec.end(), str) == files_names_vec.end()) {
				files_names_vec.push_back(str);
				launch_str = "touch " + str;
				system(launch_str.c_str());
				std::cout << "system command" << std::endl;
				std::cout << launch_str << std::endl;
			}
		}
	}
	error_file.close();
}

void MakeHTMLfromWU(MYSQL *conn, std::string wu_name_part, MOLS pair_MOLS )
{
	std::cout << "wu_name_part " << wu_name_part << std::endl;
	std::vector< std::vector<std::stringstream *> > result_vec;
	std::string str = "SELECT id FROM result WHERE workunitid IN(SELECT id FROM workunit WHERE name LIKE '%" + wu_name_part + "%')";
	std::cout << str << std::endl;
	
	ProcessQuery(conn, str, result_vec);

	if (result_vec.size() == 0) {
		std::cout << "result_vec.size() == 0" << std::endl;
		return;
	}

	std::stringstream sstream;
	std::vector<int> resultid_vec;
	std::vector<int> userid_vec;
	std::vector<std::string> username_vec;
	std::vector<int> teamid_vec;
	std::vector<std::string> teamname_vec;
	std::vector<std::string> mod_time_vec;
	unsigned u_val;
	for (unsigned i = 0; i < result_vec.size(); i++)
		for (unsigned j = 0; j < result_vec[i].size(); j++) {
			*result_vec[i][j] >> u_val;
			resultid_vec.push_back(u_val);
			delete result_vec[i][j];
		}
	result_vec.clear();

	for (std::vector<int>::iterator it = resultid_vec.begin(); it != resultid_vec.end(); it++) {
		sstream << "SELECT userid, teamid, mod_time FROM result WHERE validate_state = 1 AND id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		std::cout << str << std::endl;
		ProcessQuery(conn, str, result_vec);
		//std::cout << "workunitid " << wu_id_str << std::endl;
		std::cout << "resultid " << *it << std::endl;
		for (unsigned i = 0; i < result_vec.size(); i++) {
			*result_vec[i][0] >> u_val; // get userid
			userid_vec.push_back(u_val);
			std::cout << "userid " << u_val << std::endl;
			*result_vec[i][1] >> u_val; // get teamid
			teamid_vec.push_back(u_val);
			std::cout << "teamid " << u_val << std::endl;
			str = (*result_vec[i][2]).str(); // get mod_time
			std::cout << str << std::endl;
			std::cout << "mod_time " << str << std::endl;
			std::cout << std::endl;
			mod_time_vec.push_back(str);
			for (unsigned j = 0; j < result_vec[i].size(); j++)
				delete result_vec[i][j];
		}
		result_vec.clear();
	}

	std::cout << "teamid_vec:" << std::endl;
	for (unsigned i = 0; i < teamid_vec.size(); i++)
		std::cout << teamid_vec[i] << std::endl;

	// get names of users
	for (std::vector<int>::iterator it = userid_vec.begin(); it != userid_vec.end(); it++) {
		sstream << "SELECT name FROM user WHERE id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		std::cout << str << std::endl;
		ProcessQuery(conn, str, result_vec);

		for (unsigned i = 0; i < result_vec.size(); i++) {
			str = (*result_vec[i][0]).str();
			std::cout << str << std::endl;
			username_vec.push_back(str);
			delete result_vec[i][0];
		}
		result_vec.clear();
	}

	// get names of teams
	for (std::vector<int>::iterator it = teamid_vec.begin(); it != teamid_vec.end(); it++) {
		if (*it == 0)
			teamname_vec.push_back("");

		sstream << "SELECT name FROM team WHERE id=" << *it;
		str = sstream.str();
		sstream.clear(); sstream.str("");
		std::cout << str << std::endl;
		ProcessQuery(conn, str, result_vec);

		for (unsigned i = 0; i < result_vec.size(); i++) {
			str = (*result_vec[i][0]).str();
			std::cout << str << std::endl;
			teamname_vec.push_back(str);
			delete result_vec[i][0];
		}
		result_vec.clear();
	}

	std::ofstream MOLS_out_file("MOLS_out", std::ios_base::app);
	std::stringstream MOLS_out_sstream;

	MOLS_out_sstream << "<tr>" << std::endl << "<td> <b>" << mod_time_vec[0] << " UTC </b> </td>" << std::endl;
	MOLS_out_sstream << "<td> <a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[0] <<
		"'>" << username_vec[0] << "</a>";
	if (teamname_vec[0] != "")
		MOLS_out_sstream << " from " << teamname_vec[0];
	MOLS_out_sstream << " /" << std::endl;
	MOLS_out_sstream << "<a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[1] <<
		"'>" << username_vec[1] << "</a>";
	if (teamname_vec[1] != "")
		MOLS_out_sstream << " from " << teamname_vec[1];
	MOLS_out_sstream << " </td>" << std::endl;
	MOLS_out_sstream << "<td> diag10_2 </td>" << std::endl << "<td> \n" << "<FONT SIZE = -2>\n" << pair_MOLS.HtmlstringView() << "</FONT>\n </td>\n </tr>" << std::endl;
	MOLS_out_file << MOLS_out_sstream.str();
	MOLS_out_file.close();
	
	sstream << MOLS_out_sstream.str();
	std::cout << sstream.rdbuf();
}

void ProcessQuery(MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec)
{
	// Дескриптор результирующей таблицы
	MYSQL_RES *res;
	// Дескриптор строки
	MYSQL_ROW row;
	int num_fields;

	if (mysql_query(conn, str.c_str()) != 0)
		std::cerr << "Error: can't execute SQL-query\n";

	// Получаем дескриптор результирующей таблицы
	res = mysql_store_result(conn);

	if (res == NULL)
		std::cerr << "Error: can't get the result description\n";

	num_fields = mysql_num_fields(res);
	std::stringstream *sstream_p;
	std::vector<std::stringstream *> result_data;

	// Если имеется хотя бы одна запись - выводим список каталогов
	if (mysql_num_rows(res) > 0) {
		// В цикле перебираем все записи результирующей таблицы
		while ((row = mysql_fetch_row(res)) != NULL) {
			for (int i = 0; i < num_fields; i++) {
				//cout << row[i] << endl;
				sstream_p = new std::stringstream();
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

int getdir(std::string dir, std::vector<std::string> &files)
{
	DIR *dp;
	std::string cur_name;
	struct dirent *dirp;
	if ((dp = opendir(dir.c_str())) == NULL) {
		std::cout << std::endl << "Error in opening " << dir;
		return 1;
	}
	while ((dirp = readdir(dp)) != NULL) {
		cur_name = std::string(dirp->d_name);
		if (cur_name[0] != '.') files.push_back(cur_name);
	}
	closedir(dp);
	return 0;
}