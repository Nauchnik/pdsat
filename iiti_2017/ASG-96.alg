define len 200;

__in bit regA[31];
__in bit regB[32];
__in bit regC[33];
__out bit result[len];

void shift_regA()
{
	bit y = regA[30]^regA[6];
	for (int j = 30; j > 0; j = j-1)
	{
		regA[j] = regA[j-1];
	}
	regA[0] = y;
}

void shift_regB()
{
	bit y = regB[31]^regB[6]^regB[4]^regB[2]^regB[1]^regB[0];
	for (int j = 31; j > 0; j = j-1)
	{
		regB[j] = regB[j-1];
	}
	regB[0] = y;
}

void shift_regC()
{
	bit y = regC[32]^regC[15]^regC[3]^regC[0];
	for (int j = 32; j > 0; j = j-1)
	{
		regC[j] = regC[j-1];
	}
	regC[0] = y;
}

void main()
{
	for (int i = 0; i < len; i = i+1)
	{
		shift_regA();

		if (regA[30]&1)
		{
			shift_regB();
		}

		if (!(regA[30]&1))
		{
			shift_regC();
		}

		result[i] = regB[31]^regC[32];
	}
}