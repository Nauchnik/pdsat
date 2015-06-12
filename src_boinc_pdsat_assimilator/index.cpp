#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "Mols.h"
#include <cmath>

using namespace std;

int main()
{
	ifstream ifile("sat_output");
	string str;
	int n = 10, r = 2;
	unsigned pair_count = 0;
	while (getline(ifile, str)) {
		std::cout << std::endl << "pair_count " << pair_count << std::endl;
		/*if ((pair_count == 1) || (pair_count == 2) || (pair_count == 3)) {
			pair_count++;
			continue;
		}*/
		cout << str << endl;
		stringstream sstream;
		sstream << str;
		string sat_assign_str = "";
		while ((sat_assign_str.length() < 100) && (!sstream.eof()))
			sstream >> sat_assign_str;

		if (sat_assign_str.size() == r*pow(n, 3)) {
			std::cout << "new pair" << std::endl;
			MOLS mols(sat_assign_str, 10, 2, false);
			//mols.printToCout();
		}
		pair_count++;
	}
	ifile.close();

	return 0;
}