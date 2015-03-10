#include "Bivium.h"

Bivium::Bivium()
{}

void Bivium::reset()
{
	size_t i;
	for(i = 0; i < reg_size_; i++)
		m_reg[i] = false;
	for(i = 0; i < key_size_; i++)
		m_key[i] = false;
	for(i = 0; i < iv_size_; i++)
		m_iv[i] = false;
}

void Bivium::setKey(const std::vector<bool> &key)
{
	size_t i = 0;
	for(; i < key_size_ && i < key.size(); i++)
		m_key[i] = key[i];
	for(; i < key_size_; i++)
		m_key[i] = false;
}

void Bivium::setIV(std::vector<bool> &iv)
{
	size_t i = 0;
	for(; i < iv_size_ && i < iv.size(); i++)
		m_iv[i] = iv[i];
	for(; i < iv_size_; i++)
		m_iv[i] = false;
}

size_t Bivium::getKeySize() const
{
	return key_size_;
}

size_t Bivium::getIVSize() const
{
	return iv_size_;
}

void Bivium::init()
{
	// сброс старых значений, если они были
	for(size_t i = 0; i < reg_size_; i++)
		m_reg[i] = false;

	for(size_t i = 0; i < key_size_; i++)
		m_reg[i] = m_key[i];

	for(size_t i = 0; i < iv_size_; i++)
		m_reg[lenA + i] = m_iv[i];

	m_reg[reg_size_ - 3] = true;
	m_reg[reg_size_ - 2] = true;
	m_reg[reg_size_ - 1] = true;

	for(size_t i = 0; i < 4*reg_size_; i++)
	{
		shift_regs();
	//	std::cout << "Reg" << i << ": "; 
	//	printReg(std::cout);
	}
}

void Bivium::shift_regs()
{
	bool t1 = m_reg[65]  ^ m_reg[90]  & m_reg[91]  ^ m_reg[92]  ^ m_reg[170];
	bool t2 = m_reg[161] ^ m_reg[174] & m_reg[175] ^ m_reg[176] ^ m_reg[68];

	for(size_t i = lenA - 1; i > 0; i--)
		m_reg[i] = m_reg[i-1];
	m_reg[0] = t2;

	for(size_t i = lenB + lenA - 1; i > lenA; i--)
		m_reg[i] = m_reg[i-1];
	m_reg[lenA] = t1;
}

bool Bivium::getNextBit()
{
	bool t1 = m_reg[65]^m_reg[92];
	bool t2 = m_reg[161]^m_reg[176];
	shift_regs();
	return t1^t2;
}

bool Bivium::getStreamBit(std::vector<bool>& vec, size_t size)
{
	vec.clear();
	for(size_t i = 0; i < size; i++)
		vec.push_back(getNextBit());
	return true;
}

bool Bivium::getRegisterState(std::vector<bool>& vec)
{
	vec.clear();
	for(size_t i = 0; i < reg_size_; i++)
		vec.push_back(m_reg[i]);
	return true;
}

void Bivium::printReg(std::ostream& out)
{
	for(size_t i = 0; i < reg_size_; i++)
		out << (m_reg[i] ? "1," : "0,");
	out << std::endl;
}

int Bivium::getRegSize()
{
	return reg_size_;
}
