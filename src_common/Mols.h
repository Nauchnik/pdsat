#ifndef MOLS_H
#define MOLS_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;

class LS{
private:

	int indexof(int *a, int c){
		int res = -1;
		for (int k = 0; k<order; k++){
			if (a[k] == c){ res = k; break; }
		}
		return res;
	}
	int* lselems;
	int order;
	int &item(int i, int j);
	int * row(int i);
	int * column(int j);
	int * maindiag();
	int * subdiag();
	bool checkonlyone(int * a);
	void print(int *a){
		cout << endl;
		for (int i = 0; i<order; i++){
			cout << a[i] << " ";
		}
		cout << endl;
	}
	bool checkrow(int i)
	{
		int * t = row(i);
		return checkonlyone(t);
	}
	bool checkcolumn(int j){
		int * t = column(j);
		return checkonlyone(t);
	}
	bool checkmaindiag(){
		int * t = maindiag();
		return checkonlyone(t);
	}
	bool checksubdiag(){
		int * t = subdiag();

		return checkonlyone(t);
	}

	void transprows(int i1, int i2);
	void transpcolumns(int j1, int j2);
public:

	LS(string s, int n);
	LS(int *ss, int n);
	int it(int i, int j);
	bool parse(int * ss);
	bool check(bool diag);
	void normalizefirstrow();
	void normalizefirstcolumn();
	void reorder();
	string tostring();
};


class MOLS {
private:
	int order;
	int count;
	bool checkort(int ia, int ib);
	int checkortinc(int ia, int ib);
	int checkortinc_markings(int ia, int ib);
	int checkortinc_markings_tofile(int ia, int ib, string filename);

public:
	vector<int>markings;
	int printincompleteassumption(const char * fn);
	void print_assumptions_for_incidency_matrix(vector<int> IM, const char *fn);
	vector<int> extractassumptions(int numofsquare, int firstrow, int lastrow, int firstcolumn, int lastcolumn);
	vector<LS> Squares;
	MOLS();
	~MOLS();
	vector<int> incompleteassumptionmarkings(vector<int> current_markings);
	vector<int> incompleteassumption(vector<int> current_markings);

	void print_incomplete_assumption_tofile(vector<int> current_markings, string filename);
	void exportLStofile(const char * fn);
	MOLS(string s, int n, int r, bool wspaces);
	MOLS(vector<int> ss, int n, int r);
	void import(const char *fn, int n, int r);	
	int ortogonalitycheck();
	int ortogonalitycheck_withmarkings();
	int ortogonalitycheck_withmarkings_tofile(string filename);
	int check(bool diag);
	void reorder();
	void print(const char * fn);
};

#endif