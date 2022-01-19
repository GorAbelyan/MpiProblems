#include "Complex.h"

Complex::Complex(double real, double imag)
{
	this->real = real;
	this->imag = imag;
}
Complex Complex::operator+(Complex const& obj)
{
	Complex res;
	res.real = real + obj.real;
	res.imag = imag + obj.imag;
	return res;
}


