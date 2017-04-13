define len 144;

__in bit regA[23];
__in bit regB[24];
__in bit regC[25];

__out bit result[len];

bit stream[len];

void shift_regA()
{
	bit temp = regA[22] ^ regA[21] ^ regA[19] ^ regA[17];
	
	for (int j = 22; j > 0; j = j - 1)
	{
		regA[j] = regA[j - 1];
	}
	
	regA[0] = temp;
}

void shift_regB()
{
	bit temp = regB[0] ^ regB[1] ^ regB[8] ^ regB[9] ^ regB[15] ^ (regB[7] & regB[18]);

	for (int j = 23; j > 0; j = j - 1)
	{
		regB[j] = regB[j - 1];
	}

	regB[0] = temp;
}

void shift_regC()
{
	bit temp = regC[0] ^ regC[4] ^ regC[16] ^ (regC[1] & regC[21] & regC[23]);

	for (int j = 24; j > 0; j = j - 1)
	{
		regC[j] = regC[j - 1];
	}

	regC[0] = temp;
}

void main()
{
	int i = 0;

	for (i = 0; i < len; i = i + 1)
	{
		shift_regA();

		if (regA[22] & 1)
		{
			shift_regB();
		}

		if (regA[22] ^ 1)
		{
			shift_regC();
		}

		stream[i] = regA[22]  ^ regB[23] ^ regC[24];
	}

	for (i = 0; i < len; i = i + 1)
	{
		result[i] = stream[i];
	}
}
