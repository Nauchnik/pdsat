//#pragma once 
#ifndef ISTREAM_CIPHER_H
#define ISTREAM_CIPHER_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>

class IStreamCipher 
{
public:
	virtual ~IStreamCipher(){}
	// ��������� 
	virtual void reset() = 0;
	virtual void setKey(const std::vector<bool>& key) = 0;
	virtual void setIV(std::vector<bool>& iv) = 0;
	// ������������� ����������
	virtual void init() = 0;
	// ��������� ���������� ���� ��������� ������
	virtual bool getNextBit() = 0;
	// ��������� �������� ������� 
	virtual bool getStreamBit(std::vector<bool>& vec, size_t size) = 0;
	// �������� ������� ��������� ���������
	virtual bool getRegisterState(std::vector<bool>& vec) = 0;
	// ������� ���������� ��������� � ���������� �����
	virtual void printReg(std::ostream& out) = 0;
	virtual size_t getKeySize() const = 0;
	virtual size_t getIVSize() const = 0;
};

#endif //ISTREAM_CIPHER_H