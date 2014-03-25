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
	// интерфейс 
	virtual void reset() = 0;
	virtual void setKey(const std::vector<bool>& key) = 0;
	virtual void setIV(std::vector<bool>& iv) = 0;
	// инициализация генератора
	virtual void init() = 0;
	// генерация следующего бита ключевого потока
	virtual bool getNextBit() = 0;
	// генерация битового вектора 
	virtual bool getStreamBit(std::vector<bool>& vec, size_t size) = 0;
	// получить текущее состояние регистров
	virtual bool getRegisterState(std::vector<bool>& vec) = 0;
	// вывести содержимое регистров в переданный поток
	virtual void printReg(std::ostream& out) = 0;
	virtual size_t getKeySize() const = 0;
	virtual size_t getIVSize() const = 0;
};

#endif //ISTREAM_CIPHER_H