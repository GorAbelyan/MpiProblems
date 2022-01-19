#pragma once
class Complex
{
public:
	Complex(double real = 0.0, double imag = 0.0);
	Complex operator + (Complex const& obj);
	double real, imag;
};

