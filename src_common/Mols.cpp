#include "MOLS.h"

int ijktov(int i, int j, int k, int n){
	return i*n*n + j*n + k + 1;
}

int strtoi(string s){
	int x = atoi(s.c_str());
	return x;
}

string inttostr(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}
LS::LS(int * b, int n){
	lselems = b;
	order = n;
	cout << endl;
//	report << endl;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			cout << b[i*order + j] << " ";
			//report << b[i*order + j] << " ";
		}
		cout << endl;
	//	report << endl;
	}
}
LS::LS(string s, int n){
	order = n;
	int *a = new (int[n*n*n]);
	int *b = new (int[n*n]);
	string t;
	int i1 = 0;
	int i2 = 0;
	int k = 0;
	int g = 0;
	for (unsigned i = 0; i<s.length(); i++){
		if (s[i] == ' ') {
			i2 = i;
			t = s.substr(i1, i2 - i1);
			//cout <<t<<" ";
			a[k] = strtoi(t);
			if (a[k]>0) {/*cout<<k%n<< " ";*/b[g] = k%n; g++; }

			k++;
			//if ((k>0)&&((k%(n*n)==0))) cout<<endl;
			//if (k==n*n*n) cout <<endl<<endl;
			i1 = i2;
		}
	}
	lselems = b;
}
string LS::tostring(){
	string a;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			a = a + inttostr(item(i, j)) + ' ';
		}
		a = a + '\n';
	}
	return a;
}

int * LS::row(int i)
{
	int * a;
	a = new (int[order]);
	for (int k = 0; k<order; k++){
		a[k] = lselems[i*order + k];
	}
	return a;
}

int * LS::column(int j)
{
	int * a;
	a = new (int[order]);
	for (int k = 0; k<order; k++){
		a[k] = lselems[k*order + j];
	}
	return a;
}
int * LS::maindiag()
{
	int * a;
	a = new (int[order]);
	for (int k = 0; k<order; k++){
		a[k] = lselems[k*order + k];
	}
	return a;
}
int * LS::subdiag()
{
	int * a;
	a = new (int[order]);
	for (int k = 0; k<order; k++){
		a[k] = lselems[k*order + order - k - 1];
	}

	return a;
}

int &LS::item(int i, int j){
	return lselems[i*order + j];
}
int LS::it(int i, int j){
	int r = lselems[i*order + j];
	return r;
}

bool LS::checkonlyone(int * a){
	bool res = true;
	int *cc = new (int[order]);
	for (int k = 0; k<order; k++){
		cc[k] = 0;
	}
	for (int k = 0; k<order; k++){
		cc[a[k]]++;
		if (cc[a[k]]>1) { res = false; }
	}
	return res;
}

void LS::transprows(int i1, int i2){
	int *t = row(i1);
	for (int k = 0; k<order; k++){
		item(i1, k) = item(i2, k);
	}
	for (int k = 0; k<order; k++){
		item(i2, k) = t[k];
	}
}

void LS::transpcolumns(int j1, int j2){
	int *t = column(j1);
	for (int k = 0; k<order; k++){
		item(k, j1) = item(k, j2);
	}
	for (int k = 0; k<order; k++){
		item(k, j2) = t[k];
	}
}
bool LS::check(bool diag){

	bool a1 = true;
	for (int k = 0; k<order; k++){
		if (checkrow(k) == false){
			a1 = false;
			cout << "row " << k << " fail" << endl;
//			report << "row " << k << " fail" << endl;
		}
	}
	bool a2 = true;
	for (int k = 0; k<order; k++){
		if (checkcolumn(k) == false){
			a1 = false;
			cout << "column " << k << " fail" << endl;
//			report << "column " << k << " fail" << endl;
		}
	}
	bool a3 = true;
	if (diag == true){
		{
			if (checkmaindiag() == false){
				a3 = false;
				cout << "maindiag fail" << endl;
//				report << "maindiag fail" << endl;
			}
			if (checksubdiag() == false){
				a3 = false;
				cout << "subdiag fail" << endl;
//				report << "subdiag fail" << endl;
			}
		}
	}
	bool res = a1&&a2&&a3;
	return res;
}
void LS::normalizefirstrow(){
	int * r1 = row(1);
	for (int k = 0; k<order; k++){
		if (r1[k] != k){
			int tt1 = indexof(r1, k);
			transpcolumns(tt1, k);
			r1[tt1] = r1[k];
			r1[k] = k;
		}
	}
}
void LS::normalizefirstcolumn(){
	int * c1 = column(0);
	for (int k = 0; k<order; k++){
		if (c1[k] != k){
			int tt1 = indexof(c1, k);
			transprows(tt1, k);
			c1[tt1] = c1[k];
			c1[k] = k;
		}
	}
}
void LS::reorder(){
	int *r1 = row(0);
	int *f = new(int[order]);
	for (int k = 0; k<order; k++){
		f[r1[k]] = k;
	}
	print(r1);
	print(f);

	int *t = new (int[order*order]);
	for (int k = 0; k<order*order; k++){
		t[k] = f[lselems[k]];
		//	cout <<endl<<"("<<lselems[k]<<")->"<<f[lselems[k]]<<endl;
	}
	for (int k = 0; k<order*order; k++){
		lselems[k] = t[k];
	}
}

void MOLS::reorder(){
	for (int i = 0; i<count; i++){
		Squares[i].reorder();
	}
}
MOLS::MOLS()
{
	order = 0;
	count = 0;
}
MOLS::~MOLS(){
	//out.close();
}

vector<int> MOLS::extractassumptions(int numofsquare, int firstrow, int lastrow, int firstcolumn, int lastcolumn){
	vector<int> a;
	for (int i = firstrow; i<lastrow; i++){
		for (int j = firstcolumn; j<lastcolumn; j++){
			int t = Squares[numofsquare].it(i, j);
			int r = numofsquare *order*order*order + i*order*order + j*order + t + 1;
			a.push_back(r);
		}
	}
	return a;
}
void MOLS::import(const char*fn, int n, int r){
	ifstream myfile(fn);
	int *b;
	order = n;
	count = r;
	b = new (int[n*n*r]);
	int sc = 0;
	if (myfile.is_open())
	{
		while (!myfile.eof())//( myfile.good() )
		{
			string line = "";
			getline(myfile, line);
			for (unsigned i = 0; i<line.length(); i++){
				if ((line[i] >= '0') && (line[i]<'0' + 10)){
					b[sc] = line[i] - '0';
					sc++;
				}
			}
		}
		for (int k = 0; k<r; k++){
			int *z = new (int[n*n]);
			for (int t = 0; t<n*n; t++){
				z[t] = b[k*n*n + t];
			}
			LS h = LS(z, order);
			Squares.push_back(h);
		}
	}

}

struct index_int{
	int i;
	int j;
};

int MOLS::printincompleteassumption(const char * fn){

	LS a = Squares[1];
	LS b = Squares[2];
	int count = 0;
	cout << endl << "PAIRS" << endl;
	int * r = new (int[order*order]);
	vector<int> pairs;
	for (int k = 0; k<order*order; k++){ r[k] = 0; }
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			int h = a.it(i, j)*order + b.it(i, j);
			cout << h << " ";
			pairs.push_back(h);
			r[h]++;
			if (r[h] == 1){ count++; }
		}
		cout << endl;
	}
	vector<int> indexesweneed;
	cout << endl;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			int h = a.it(i, j)*order + b.it(i, j);
			if (r[h] != 1) {
				pairs[i*order + j] = 0;
				indexesweneed.push_back(0);
				cout << 0 << " ";
			}
			if (r[h] == 1) {
				index_int a;
				a.i = i;
				a.j = j;
				indexesweneed.push_back(1);
				cout << 1 << " ";
				pairs[i*order + j] = 0;
			}
		}
		cout << endl;
	}
	print_assumptions_for_incidency_matrix(indexesweneed, fn);
	return count;
}
vector<int> MOLS::incompleteassumption(vector<int> current_markings){
	// returns values of markings to be put to CNF
	vector<int> result;
	for (unsigned i = 0; i<current_markings.size(); i++){
		if (current_markings[i] == 1){
			int temp_var = 3 * order*order*order + i + 1;
			result.push_back(temp_var);
		}
	}	
	//using current_markings returns values of existing square
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			if ((current_markings[i*order + j] == 1) && (markings[i*order + j] == 1)){
				for (int u = 0; u<count; u++){
					int temp_var = u*order*order*order + i*order*order + j*order + Squares[u].it(i, j) + 1;
					result.push_back(temp_var);
				}
			}
		}
	}
	return result;
}
void MOLS::print_incomplete_assumption_tofile(vector<int> current_markings, string filename){
	ofstream out;
	out.open(filename.c_str());
	out << "Current Markings" << endl << endl;
	for (int i = 0; i < order; i++){
		for (int j = 0; j < order; j++){
			out << current_markings[i*order + j] << " ";
		}
		out << endl;
	}
	out << endl << endl << "Squares assumptions" << endl << endl;
	for (int k = 0; k < count; k++){
		for (int i = 0; i<order; i++){
			for (int j = 0; j<order; j++){
				if ((current_markings[i*order + j] == 1) && (markings[i*order + j] == 1)){
					out << Squares[k].it(i, j) << " ";
				}
				else { out << "x "; }

			}
			out << endl;
		}
	}
	out.close();
}

vector<int> MOLS::incompleteassumptionmarkings(vector<int> current_markings){
	// returns values of markings to be put to CNF
	cout << endl;
	vector<int> result;
	for (unsigned i = 0; i<current_markings.size(); i++){
		if (current_markings[i] == 1){
			int temp_var = 3 * order*order*order + i + 1;
			result.push_back(temp_var);
		}
	}
	return result;
}
void MOLS::print_assumptions_for_incidency_matrix(vector<int> IM, const char *fn){
	ofstream out(fn);
	//out.open(fn);
	cout << "Assumptions" << endl;
	//	out<<"Assumptions"<<endl;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			if (IM[i*order + j] == 1){
				for (int r = 0; r<3; r++){
					cout << r* order * order * order + ijktov(i, j, Squares[r].it(i, j), order) << " 0" << endl;
					out << r* order * order * order + ijktov(i, j, Squares[r].it(i, j), order) << " 0" << endl;
				}
			}
		}
	}
	out.close();
}
bool MOLS::checkort(int ia, int ib){
	bool f = true;
	LS a = Squares[ia];
	LS b = Squares[ib];
	int * r = new (int[order*order]);
	for (int k = 0; k<order*order; k++){ r[k] = 0; }
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			int h = a.it(i, j)*order + b.it(i, j);
			r[h]++;
			if (r[h]>1){
				f = false;
				cout << "orthogonality check failed" << endl;
//				report << "orthogonality check failed" << endl;
			}
		}
	}
	return f;
}
int MOLS::checkortinc(int ia, int ib){
	LS a = Squares[ia];
	LS b = Squares[ib];
	int count = 0;
	int * r = new (int[order*order]);
	for (int k = 0; k<order*order; k++){ r[k] = 0; }
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			int h = a.it(i, j)*order + b.it(i, j);
			r[h]++;
			if (r[h] == 1){ count++; }
		}
	}

	if (count<order*order){
		cout << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
//		report << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
		for (int i = 0; i<order; i++){
			for (int j = 0; j<order; j++){
				cout << r[i*order + j] << " ";
//				report << r[i*order + j] << " ";
			}
			cout << endl;
//			report << endl;
		}
	}
	return count;
}

int MOLS::ortogonalitycheck(){
	int res = 100;
	for (int i = 0; i<count; i++){
		for (int j = i + 1; j<count; j++){
			int r = checkortinc(i, j);
			if (r<order*order){
				cout << "Squares " << i + 1 << " and " << j + 1 << " are orthogonal in " << r << " entries" << endl;
//				report << "Squares " << i + 1 << " and " << j + 1 << " are orthogonal in " << r << " entries" << endl;

			}
			if (res>r){ res = r; }
		}
	}
	return res;
}

int MOLS::checkortinc_markings(int ia, int ib){
	LS a = Squares[ia];
	LS b = Squares[ib];
	int count = 0;
	int * r = new (int[order*order]);
	for (int k = 0; k<order*order; k++){ r[k] = 0; }
	int markings_numberofones = 0;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			if (markings[i*order + j] == 1){
				markings_numberofones++;
				int h = a.it(i, j)*order + b.it(i, j);
				r[h]++;
				if (r[h] == 1){ count++; }
			}
		}
	}

	//if (count<markings_numberofones){
	cout << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
//	report << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			cout << r[i*order + j] << " ";
//			report << r[i*order + j] << " ";
		}
		cout << endl;
//		report << endl;
	}
	//}
	return count;
}
int MOLS::checkortinc_markings_tofile(int ia, int ib, string filename){
	ofstream out;
	out.open(filename.c_str(), ios::app);
	LS a = Squares[ia];
	LS b = Squares[ib];
	int count = 0;
	int * r = new (int[order*order]);
	for (int k = 0; k<order*order; k++){ r[k] = 0; }
	int markings_numberofones = 0;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			if (markings[i*order + j] == 1){
				markings_numberofones++;
				int h = a.it(i, j)*order + b.it(i, j);
				r[h]++;
				if (r[h] == 1){ count++; }
			}
		}
	}

	//if (count<markings_numberofones){
	cout << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
//	report << endl << "Orthogonality Matrix for squares " << ia + 1 << " and " << ib + 1 << endl;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			out << r[i*order + j] << " ";
		}
		out << endl;
	}
	//}
	out.close();
	return count;
}
int MOLS::ortogonalitycheck_withmarkings(){
	int res = 100;
	for (int i = 0; i<count; i++){
		for (int j = i + 1; j<count; j++){
			int r = checkortinc_markings(i, j);
			if (r<order*order){
				cout << "Squares " << i + 1 << " and " << j + 1 << " are orthogonal in " << r << " entries (according to markings)" << endl;
//				report << "Squares " << i + 1 << " and " << j + 1 << " are orthogonal in " << r << " entries (according to markings)" << endl;
			}
			if (res>r){ res = r; }
		}
	}
	return res;

}


int MOLS::ortogonalitycheck_withmarkings_tofile(string filename){
	ofstream out;
	out.open(filename.c_str(), ios::app);
	for (unsigned i = 0; i<Squares.size(); i++){
		out << Squares[i].tostring() << endl;
	}
	out << endl << "Markings" << endl;

	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			out << markings[i*order + j] << " ";
		}
		out << endl;
	}
	out << endl;
	int res = order*order;
	for (int i = 0; i<count; i++){
		for (int j = i + 1; j<count; j++){
			int r = checkortinc_markings(i, j);
			if (r<order*order){
				out << "Squares " << i + 1 << " and " << j + 1 << " are orthogonal in " << r << " entries (according to markings)" << endl;
			}
			if (res>r){ res = r; }
		}
	}

	out.close();
	return res;

}

int MOLS::check(bool diag){
	bool f = true;
	int r = 0;
	for (int i = 0; i<count; i++){
		LS a = Squares[i];
		bool r = a.check(diag);
		if (r == false){
			f = false;
			cout << "square " << i << " is wrong" << endl;
//			report << "square " << i << " is wrong" << endl;
		}
	}

	if (f == true){
		r = ortogonalitycheck();
		return r;
	}
	else {
		return 0;
	}
}

MOLS::MOLS(string s, int n, int r, bool wspaces){
	int *a = new (int[n*n*n*n*n]);
	int *b = new (int[n*n*r]);
	string t;
	int i1 = 0;
	int i2 = 0;
	int k = 0;
	int g = 0;

	if (wspaces == true){
		for (unsigned i = 0; i<s.size(); i++){
			if (s[i] == ' ')
			{
				i2 = i;
				t = s.substr(i1, i2 - i1);
				//cout <<t<<" ";
				a[k] = strtoi(t);
				k++;
				i1 = i2;
			}
		}
	}
	else {
		for (unsigned i = 0; i<s.length(); i++){
			//a[i]=s[i]-'0';
			if (s[i] == '0'){ a[i] = 0; }
			else{ a[i] = 1; }
		}
	}
	for (int i = 0; i<n*n*n*r; i++){
		if (a[i]>0) {
			b[g] = i%n;
			g++;
		}
	}

	cout << endl;
	order = n;
	count = r;
	for (int k = 0; k<r; k++){
		int *z = new (int[n*n]);
		for (int t = 0; t<n*n; t++){
			z[t] = b[k*n*n + t];
		}
		LS h = LS(z, order);
		Squares.push_back(h);
	}
	cout << endl << "Markings:" << endl;
	for (int i = n*n*n*r; i<n*n*n*r + n*n; i++){
		cout << (a[i]>0) << " ";
		if (a[i]>0){
			markings.push_back(1);
		}
		else{
			markings.push_back(0);
		}
		if ((i + 1) % n == 0) cout << endl;
	}

}

MOLS::MOLS(vector<int> ss, int n, int r){
	int *b = new (int[n*n*r]);
	int g = 0;

	for (int i = 0; i<n*n*n*r; i++){
		if (ss[i]>0) {
			b[g] = i%n;
			g++;
		}
	}
	cout << endl;
	order = n;
	count = r;
	for (int k = 0; k<r; k++){
		int *z = new (int[n*n]);
		for (int t = 0; t<n*n; t++){
			z[t] = b[k*n*n + t];
		}
		LS h = LS(z, order);
		Squares.push_back(h);
	}
	cout << endl << "Markings:" << endl;
	for (int i = n*n*n*r; i<n*n*n*r + n*n; i++){
		cout << (ss[i]>0) << " ";
		if (ss[i]>0){
			markings.push_back(1);
		}
		else{
			markings.push_back(0);
		}
		if ((i + 1) % n == 0) cout << endl;
	}

}
void MOLS::print(const char * fn){
	ofstream myfile(fn);
	if (myfile.is_open()){
		for (int k = 0; k<count; k++){
			myfile << Squares[k].tostring() << '\n';
			myfile << "reordered version:\n";
			Squares[k].reorder();
			myfile << Squares[k].tostring() << '\n';
		}
	}
}
void MOLS::exportLStofile(const char * fn){
	ofstream myfile(fn);
	if (myfile.is_open()){
		for (int i = 0; i<count; i++){
			LS a = Squares[i];
			for (int l = 0; l<order; l++){
				for (int j = 0; j<order; j++){
					int tmp = a.it(l, j);
					myfile << order*order*order*i + order*order*l + order*j + tmp + 1 << " 0" << endl;
				}
			}

		}
	}
	myfile.close();
}


bool eq(LS a, LS b, int order){
	bool res = true;
	for (int i = 0; i<order; i++){
		for (int j = 0; j<order; j++){
			if (a.it(i, j) != b.it(i, j)){ res = false; }
		}
	}
	return res;
}
struct ijk{
	int sq;
	int i;
	int j;
	int k;
	int phase;

};

struct less_than_ijk
{
	inline bool operator() (const ijk& st1, const ijk& st2)
	{
		bool r;
		if (st1.sq<st2.sq){ r = true; }
		else{
			if (st1.sq>st2.sq) { r = false; }
			else {
				if (st1.i<st2.i){ r = true; }
				else if (st1.i>st2.i){ r = false; }
				else {
					if (st1.j<st2.j){ r = true; }
					else if (st1.j>st2.j){ r = false; }
					else{
						if (st1.k<st2.k){ r = true; }
						else if (st1.k >= st2.k){ r = false; }
					}
				}
			}
		}
		return r;
	}
};

