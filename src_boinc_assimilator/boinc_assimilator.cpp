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

void ProcessQuery(MYSQL *conn, std::string str, std::vector< std::vector<std::stringstream *> > &result_vec);
void MakeHTMLfromWU(MYSQL *conn, std::string wu_id_str);
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
	std::string wu_part_name;
	std::vector< std::vector<std::stringstream *> > result_vec;
	while (getline(sat_file, str)) {
		sstream << str;
		sstream >> wu_part_name;
		sstream.str(""); sstream.clear();
		MakeHTMLfromWU(conn, wu_part_name);
	}
	sat_file.close();
	
	std::cout << "*** done" << std::endl;
}

void MakeHTMLfromWU(MYSQL *conn, std::string wu_name_part)
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

	sstream << "<tr>" << std::endl << "<td> <b>" << mod_time_vec[0] << " UTC </b> </td>" << std::endl;
	sstream << "<td> <a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[0] <<
		"'>" << username_vec[0] << "</a>";
	if (teamname_vec[0] != "")
		sstream << " from " << teamname_vec[0];
	sstream << " /" << std::endl;
	sstream << "<a href = 'http://sat.isa.ru/pdsat/show_user.php?userid=" << userid_vec[1] <<
		"'>" << username_vec[1] << "</a>";
	if (teamname_vec[1] != "")
		sstream << " from " << teamname_vec[1];
	sstream << " </td>" << std::endl;
	sstream << "<td>Bivium</td>" << std::endl << "<td> </td>" << std::endl << "</tr>" << std::endl;

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