//#define len 128;
#include <vector>

using namespace std;

const int regAlen = 19;
const int regBlen = 22;
const int regClen = 23;
const int regDlen = 17;
const int sbit1 = 3;
const int sbit2 = 7;
const int sbit3 = 10;

bool regA[regAlen];
bool regB[regBlen];
bool regC[regClen];
bool regD[regDlen];

class a5_2 {
public:
	a5_2();
	bool shift_rslosA();
	bool shift_rslosB();
	bool shift_rslosC();
	bool shift_rslosD();
	bool majority( bool A, bool B, bool C );
	bool getNextBit();
	void setKey( const std::vector<bool> &key );
	void getState( vector<bool> &state );
	unsigned sum_reg_len;
};

a5_2 :: a5_2 ()
{
	sum_reg_len = regAlen + regBlen + regClen + regDlen;
}

bool a5_2 :: shift_rslosA(){
	bool x = regA[18];
	bool y = regA[18]^regA[17]^regA[16]^regA[13];
	for(int j = 18; j > 0; j=j-1){
		regA[j] = regA[j-1];
	}
	regA[0] = y;
	return x;
}

bool a5_2 :: shift_rslosB(){
	bool x = regB[21];
	bool y = regB[21]^regB[20];
	for(int j = 21; j > 0; j=j-1){
		regB[j] = regB[j-1];
	}
	regB[0] = y;
	return x;
}

bool a5_2 :: shift_rslosC(){
	bool x = regC[22];
	bool y = regC[22]^regC[21]^regC[20]^regC[7];
	for(int j = 22; j > 0; j=j-1){
		regC[j] = regC[j-1];
	}
	regC[0] = y;
	return x;
}

bool a5_2 :: shift_rslosD(){
	bool x = regD[16];
	bool y = regD[16]^regD[11];
	for(int j = 22; j > 0; j=j-1){
		regD[j] = regD[j-1];
	}
	regD[0] = y;
	return x;
}

bool a5_2 :: majority(bool A, bool B, bool C){
    return A&B|A&C|B&C;
}

bool a5_2 :: getNextBit() {
	bool maj_regD = majority(regD[sbit1], regD[sbit2], regD[sbit3]);
	if(!(maj_regD ^ regD[sbit3])) shift_rslosA();
	if(!(maj_regD ^ regD[sbit1])) shift_rslosB();
	if(!(maj_regD ^ regD[sbit2])) shift_rslosC();				
	
	// majority functions of registers
	bool maj_regA = majority(regA[12], regA[14], regA[15]);
	bool maj_regB = majority(regB[9], regB[13], regB[16]);
	bool maj_regC = majority(regC[13], regC[16], regC[18]);
	bool reg_xor = regA[18] ^ regB[21] ^ regC[22];
	
	return reg_xor ^ maj_regA ^ maj_regB ^ maj_regC;
}

void a5_2 :: setKey( const std::vector<bool> &key )
{
	if ( key.size() < regAlen + regBlen + regClen + regDlen ) {
		cerr << "key.size() < sum length of registers a5/2" << endl;
		exit(1);
	}
	unsigned i = 0;
	for(; i < regAlen && i < key.size(); i++)
		regA[i] = key[i];
	for(; i < regBlen && i < key.size(); i++)
		regB[i] = key[i];
	for(; i < regClen && i < key.size(); i++)
		regC[i] = key[i];
	for(; i < regDlen && i < key.size(); i++)
		regD[i] = key[i];
}

void a5_2 :: getState( vector<bool> &state )
{
	state.clear();
	for(unsigned i = 0; i < regAlen; i++)
		state.push_back( regA[i] );
	for(unsigned i = 0; i < regBlen; i++)
		state.push_back( regB[i] );
	for(unsigned i = 0; i < regClen; i++)
		state.push_back( regC[i] );
	for(unsigned i = 0; i < regDlen; i++)
		state.push_back( regD[i] );
}