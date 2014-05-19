define len 128;

bool regA[19];
bool regB[22];
bool regC[23];
bool regD[17];

bool result[len];

class a5_2{

bit shift_rslosA(){
	bit x = regA[18];
	bit y = regA[18]^regA[17]^regA[16]^regA[13];
	for(int j = 18; j > 0; j=j-1){
		regA[j] = regA[j-1];
	}
	regA[0] = y;
	return x;
}

bit shift_rslosB(){
	bit x = regB[21];
	bit y = regB[21]^regB[20];
	for(int j = 21; j > 0; j=j-1){
		regB[j] = regB[j-1];
	}
	regB[0] = y;
	return x;
}

bit shift_rslosC(){
	bit x = regC[22];
	bit y = regC[22]^regC[21]^regC[20]^regC[7];
	for(int j = 22; j > 0; j=j-1){
		regC[j] = regC[j-1];
	}
	regC[0] = y;
	return x;
}

bit shift_rslosD(){
	bit x = regD[16];
	bit y = regD[16]^regD[11];
	for(int j = 22; j > 0; j=j-1){
		regD[j] = regD[j-1];
	}
	regD[0] = y;
	return x;
}

bit majority(bit A, bit B, bit C){
    return A&B|A&C|B&C;
}

}

/*
void main(){

	// номера битов синхронизации
	int sbit1 = 3;
	int sbit2 = 7;
	int sbit3 = 10;

	for(int i = 0; i < len; i = i + 1)
	{
		bit maj_regD = majority(regD[sbit1], regD[sbit2], regD[sbit3]);
		if(!(maj_regD ^ regD[sbit3])) shift_rslosA();
		if(!(maj_regD ^ regD[sbit1])) shift_rslosB();
		if(!(maj_regD ^ regD[sbit2])) shift_rslosC();				

		// мажоритарные функции регистров
		__mem bit maj_regA = majority(regA[12], regA[14], regA[15]);
		__mem bit maj_regB = majority(regB[9], regB[13], regB[16]);
		__mem bit maj_regC = majority(regC[13], regC[16], regC[18]);
		__mem bit reg_xor = regA[18] ^ regB[21] ^ regC[22];
		
		result[i] = reg_xor ^ maj_regA ^ maj_regB ^ maj_regC;
	}
}
*/
