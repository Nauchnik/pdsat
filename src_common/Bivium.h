#ifndef BIVIUM_H
#define BIVIUM_H

#include "IStreamCipher.h"
//#include "stdafx.h"

/*
* ������ ���������� ������������� ������ ��� ��������� �������� ������ - ��� (����, �������� �����).
* ����������� ���� ������������ � ������ ��� ��������� ������.
*/ 
class Bivium: public IStreamCipher{
public:
	Bivium();
	void reset();
	int getRegSize();
	// ��������� ����� � ����������������� �������
	void setKey(const std::vector<bool>& key);
	void setIV(std::vector<bool>& iv);
	size_t getKeySize() const;
	size_t getIVSize() const;
	// ������������� ����������
	void init();
	// �������� ��������� ����������
	//bool isReady();
	// ��������� ���������� ���� ��������� ������
	bool getNextBit();
	// ��������� �������� ������� 
	bool getStreamBit(std::vector<bool>& vec, size_t size);
	// �������� ������� ��������� ���������
	bool getRegisterState(std::vector<bool>& vec);
	// for debug
	void printReg(std::ostream& out);
private:
	void shift_regs();
private:
	// ������ ������� ������� - 177 ������� �����
	static const size_t reg_size_ = 177;
	static const size_t lenA = 93;
	static const size_t lenB = 84;
	static const size_t key_size_ = 80;
	static const size_t iv_size_ = 80;

	bool m_reg[reg_size_];
	bool m_key[key_size_];
	bool m_iv[iv_size_];
};

#endif // BIVIUM_H