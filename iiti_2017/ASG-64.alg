define len 64;

__in bit regA[19];
__in bit regB[22];
__in bit regC[23];

__out bit result[len];



bit stream[len];



void shift_regA()
{
	bit temp = regA[18] ^ regA[17] ^ regA[16] ^ regA[13];
	
	for (int j = 18; j > 0; j = j - 1)
	{
		regA[j] = regA[j - 1];
	}
	
	regA[0] = temp;
}

void shift_regB()
{
	bit temp = regB[21] ^ regB[20];

	for (int j = 21; j > 0; j = j - 1)
	{
		regB[j] = regB[j - 1];
	}

	regB[0] = temp;
}

void shift_regC()
{
	bit temp = regC[22] ^ regC[21] ^ regC[20] ^ regC[7];

	for (int j = 22; j > 0; j = j - 1)
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

		if (regA[18] & 1)
		{
			shift_regB();
		}

		if (regA[18] ^ 1)
		{
			shift_regC();
		}

		stream[i] = regB[21] ^ regC[22];
	}

	for (i = 0; i < len; i = i + 1)
	{
		result[i] = stream[i];
	}
}