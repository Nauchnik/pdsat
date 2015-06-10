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
	/*omp_set_num_threads(2);
	#pragma omp parallel
	for (long long i = 0; i < 1e50; i++) {
		int tnum = omp_get_thread_num();
		if (tnum != 0)
			system("pause");
		else
			exit(0);
	}*/

	bool isAdding = true;
	std::vector<double> cws_vi{5,6,3,1};
	for (unsigned cws_vi_index = 0; cws_vi_index < cws_vi.size() - 1; cws_vi_index++)
		if (cws_vi[cws_vi_index] < cws_vi[cws_vi_index + 1])
			isAdding = false;
	std::cout << isAdding;

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