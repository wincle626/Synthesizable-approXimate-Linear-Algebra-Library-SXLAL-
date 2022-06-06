/*
 * mathfunc.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_MATH_MATHFUNC_HPP_
#define HEADER_MATH_MATHFUNC_HPP_

#include "data.hpp"
#include <cmath>

template<class T>
void complex_add(Complex<T> a, Complex<T> b, Complex<T> &c){
	c.real = a.real + b.real;
	c.imag = a.imag + b.imag;
}

template<class T>
void complex_sub(Complex<T> a, Complex<T> b, Complex<T> &c){
	c.real = a.real - b.real;
	c.imag = a.imag - b.imag;
}

template<class T>
void complex_mul(Complex<T> a, Complex<T> b, Complex<T> &c){
	c.real = a.real * b.real - a.imag * b.imag;
	c.imag = a.real * b.imag + a.imag * b.real;
}

template<class T>
void complex_div(Complex<T> a, Complex<T> b, Complex<T> &c){
	double x = a.real*b.real + a.imag*b.imag;
	double y = a.imag*b.real - a.real*b.imag;
	double z = b.real*b.real + b.imag*b.imag;
	c.real = x/z;
	c.imag = y/z;
}

template<class T>
void complex_conj(Complex<T>a, Complex<T> &b){
	b.real = a.real;
	b.imag = a.imag * -1;
}

template<class T>
void complex_abs(Complex<T> a, T &b){
	b = std::sqrt(a.real * a.real + a.imag * a.imag);
}

template<class T>
void complex_abs2(Complex<T> a, T &b){
	b = a.real * a.real + a.imag * a.imag;
}

template<class T>
void exp2complex(T phi, Complex<T> &a){
	a.real = std::cos(phi);
	a.imag = std::sin(phi);
}


#endif /* HEADER_MATH_MATHFUNC_HPP_ */
