/*
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_DATA_HPP_
#define SRC_DATA_HPP_

// Data structure
#define ROW 4 // matrix row number
#define COL 4 // matrix column number
#define DIAG (COL<ROW ? COL : ROW) // diagonal matrix size
#define DIAG_VALUE 50 // diagonal matrix value scale
#define DIAG_RATIO 0.5 // diagonal matrix sparse ratio
#define SPARSE_RATIO 0.5 // vector sparse ratio

#define SIZE (COL>ROW ? COL : ROW) // square matrix size

#define ROW1 230 // left multiply matrix row number
#define COL1 220 // left multiply matrix column number
#define ROW2 COL1 // right multiply matrix row number
#define COL2 240 // right multiply matrix column number

// Data type
#define FLOAT_SIZE RAND_MAX // floating point fraction scale
#define INTEGER_SCALE 10 // floating point integer scale

const double PI = 3.141592653589793238460; // PI constant

// Complex data structure
template<typename T>
struct Complex{
	T real;
	T imag;
};
template<typename T>
using cmplx = Complex<T>;

// Complex data class
template<class T>
class CMPLX{
public:
	T real;
	T imag;
	CMPLX(T r = 0) {real = r; imag = 0;}
	CMPLX(T r = 0, T i =0)  {real = r;   imag = i;}

	CMPLX<T> operator = (CMPLX<T> const &obj) {
		CMPLX<T> res;
        res.real = obj.real;
        res.imag = obj.imag;
        return res;
   }

	CMPLX<T> operator + (CMPLX<T> const &obj) {
		CMPLX<T> res;
        res.real = real + obj.real;
        res.imag = imag + obj.imag;
        return res;
   }

	CMPLX<T> operator - (CMPLX<T> const &obj) {
		CMPLX<T> res;
        res.real = real - obj.real;
        res.imag = imag - obj.imag;
        return res;
   }

	CMPLX<T> operator * (CMPLX<T> const &obj) {
		CMPLX<T> res;
        res.real = real * obj.real - imag * obj.imag;
        res.imag = real * obj.imag + imag * obj.real;
        return res;
   }

	CMPLX<T> operator / (CMPLX<T> const &obj) {
		CMPLX<T> res;
		T x = real*obj.real + imag*obj.imag;
		T y = imag*obj.real - real*obj.imag;
		T z = obj.real*obj.real + obj.imag*obj.imag;
		res.real = x/z;
		res.imag = y/z;
        return res;
   }

};

#endif /* SRC_DATA_HPP_ */
