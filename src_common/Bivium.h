#ifndef BIVIUM_H
#define BIVIUM_H

#include "IStreamCipher.h"
//#include "stdafx.h"

/*
* Данная реализация предназначена только для генерации тестовых данных - пар (ключ, ключевой поток).
* Построенные пары поставляются в шаблон для получения тестов.
*/ 
class Bivium: public IStreamCipher{
public:
	Bivium();
	void reset();
	int getRegSize();
	// установка ключа и инициализирующего вектора
	void setKey(const std::vector<bool>& key);
	void setIV(std::vector<bool>& iv);
	size_t getKeySize() const;
	size_t getIVSize() const;
	// инициализация генератора
	void init();
	// проверка состояния генератора
	//bool isReady();
	// генерация следующего бита ключевого потока
	bool getNextBit();
	// генерация битового вектора 
	bool getStreamBit(std::vector<bool>& vec, size_t size);
	// получить текущее состояние регистров
	bool getRegisterState(std::vector<bool>& vec);
	// for debug
	void printReg(std::ostream& out);
private:
	void shift_regs();
private:
	// размер рабочей области - 177 битовых ячеек
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