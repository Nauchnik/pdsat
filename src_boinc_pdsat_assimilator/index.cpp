#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "Mols.h"
//#include <omp.h>

using namespace std;

int main()
{
	ifstream ifile("sat_output");
	string str;
	getline(ifile,str);
	ifile.close();
	cout << str << endl;
	stringstream sstream;
	sstream << str;
	string sat_assign_str = "";
	while ((sat_assign_str.length() < 100) && (!sstream.eof()) )
		sstream >> sat_assign_str;

	int n = 10;
	int r = 2;
	MOLS mols(sat_assign_str, n, r, false);
	mols.print("squares");

	return 0;
}