/*
 * fft.hpp
 *
 *  Created on: 24 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_FFT_HPP_
#define SRC_FFT_HPP_

#include "common.hpp"
#include "data.hpp"
#include "fpt_algebra.hpp"
#include <fftw3.h>

/*void complex_add(Complex<double> a, Complex<double> b, Complex<double> &c);
void complex_sub(Complex<double> a, Complex<double> b, Complex<double> &c);
void complex_mul(Complex<double> a, Complex<double> b, Complex<double> &c);
void complex_div(Complex<double> a, Complex<double> b, Complex<double> &c);
void complex_abs(Complex<double> a, double &b);
void complex_abs2(Complex<double> a, double &b);
void exp2complex(double phi, Complex<double> &a);
unsigned int bitreverse(unsigned int num);
void swapvecelement(Complex<double> Vec[DIAG], int i, int j);
void bitreverse_reorder_1d(Complex<double> Vec[DIAG]);
void fftw3_rfft_1d(double IN[DIAG], double OUT[DIAG]);
void fftw3_cfft_1d(Complex<double> IN[DIAG], Complex<double> OUT[DIAG]);
void fftw3_cifft_1d(Complex<double> IN[DIAG], Complex<double> OUT[DIAG]);
void fftw3_cfft_2d(Complex<double> IN[DIAG][DIAG], Complex<double> OUT[DIAG][DIAG]);
void fftw3_cifft_2d(Complex<double> IN[DIAG][DIAG], Complex<double> OUT[DIAG][DIAG]);
std::vector<std::complex<double>> fft_recursive_cpp(std::vector<std::complex<double>> &x);
void fft_radix2_iterative_cpp_1d(std::vector<std::complex<double>> x,
					   std::vector<std::complex<double>> &y,
					   int log2n);
void ifft_radix2_iterative_cpp_1d(std::vector<std::complex<double>> x,
						std::vector<std::complex<double>> &y,
						int log2num);
std::vector<std::complex<double>> fft_radix2_recursive_cpp_1d(std::vector<std::complex<double>> &x);
std::vector<std::complex<double>> ifft_radix2_recursive_cpp_1d(std::vector<std::complex<double>> &x);
void fft_radix2_iterative_c_1d(Complex<double> x[DIAG],
					 Complex<double> y[DIAG],
					 int log2n);
void ifft_radix2_iterative_c_1d(Complex<double> x[DIAG],
					  Complex<double> y[DIAG],
					  int log2n);
void fft_bluestein_iterative_c_1d(Complex<double> x[DIAG],
						Complex<double> y[DIAG]);
void ifft_bluestein_iterative_c_1d(Complex<double> x[DIAG],
						Complex<double> y[DIAG]);
void fft_radix2_iterative_c_2d(Complex<double> x[DIAG][DIAG],
				 	 Complex<double> y[DIAG][DIAG],
					 int log2n);
void ifft_radix2_iterative_c_2d(Complex<double> x[DIAG][DIAG],
					  Complex<double> y[DIAG][DIAG],
					  int log2n);
*/

class FFT_C{
public:

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

	template<class T>
	T bitreverse(T num, unsigned int bitsnum)
	{
	    unsigned int  NO_OF_BITS = bitsnum;
	    T reverse_num = 0;
	    unsigned int i;
	    for (i = 0; i < NO_OF_BITS; i++)
	    {
	        if((num & (1 << i)))
	           reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
	   }
	    return reverse_num;
	}

	// Iterative FFT function to compute the DFT
	// of given coefficient vector
	template<class T, int LEN>
	void fft_radix2_iterative_c_1d(Complex<T> x[LEN],
					 	 	Complex<T> y[LEN],
							int log2n){
	    int n = LEN;

	    // bit reversal of the given array
	    for (int i = 0; i < n; ++i) {
	        int rev = bitreverse<unsigned int>(i, log2n);
	        y[i].real = x[rev].real;
	        y[i].imag = x[rev].imag;
	    }

	    // j is iota
	    cmplx<T> tmp;
	    for (int s = 1; s <= log2n; ++s) {
	        int m = 1 << s; // 2 power s
	        int m2 = m >> 1; // m2 = m/2 -1
	        cmplx<T> w;
	        w.real = 1; w.imag = 0;

	        // principle root of nth complex
	        // root of unity.
	        cmplx<T> wm;
	        wm.real = (T) std::cos(PI / m2);
	        wm.imag = (T) std::sin(PI / m2);
	        for (int j = 0; j < m2; ++j) {
	            for (int k = j; k < n; k += m) {

	                // t = twiddle factor
	            	cmplx<T> t;
	            	complex_mul<T>(w, y[k+m2], t);
	            	cmplx<T> u;
	            	u.real = y[k].real; u.imag = y[k].imag;

	                // similar calculating y[k]
	            	complex_add<T>(u, t, tmp);
	            	y[k].real = tmp.real;
	            	y[k].imag = tmp.imag;

	                // similar calculating y[k+n/2]
	            	complex_sub<T>(u, t, tmp);
	            	y[k + m2].real = tmp.real;
	            	y[k + m2].imag = tmp.imag;
	            }
	            complex_mul<T>(w, wm, tmp);
	            w.real = tmp.real; w.imag = tmp.imag;
	        }
	    }
	}

	template<class T, int LEN>
	void ifft_radix2_iterative_c_1d(Complex<T> x[LEN],
						  Complex<T> y[LEN],
						  int log2n)
	{
	    // conjugate the complex numbers
		for(int i=0;i<LEN;i++)
			complex_conj<T>(x[i],x[i]);

	    // forward fft
	    fft_radix2_iterative_c_1d<T, LEN>( x, y, log2n );

	    // conjugate the complex numbers again
		for(int i=0;i<LEN;i++)
			complex_conj<T>(y[i],y[i]);

	    // scale the numbers
		for(int i=0;i<LEN;i++){
			y[i].real /= LEN; y[i].imag /= LEN;
		}
	}

	// Iterative FFT function to compute the DFT
	// of given coefficient vector
	template<class T, int M, int N>
	void fft_radix2_iterative_c_2d(Complex<double> x[M][N],
					 	 Complex<double> y[M][N],
						 int log2n){
		// clone input
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				y[i][j].real = x[i][j].real;
				y[i][j].imag = x[i][j].imag;
			}
		}

		Complex<T> vec1[N];
		Complex<T> vec_fft1[N];
		Complex<T> vec2[M];
		Complex<T> vec_fft2[M];
		// do the row based fft
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				vec1[j] = y[i][j];
			}
			fft_radix2_iterative_c_1d<T, N>(vec1,vec_fft1,log2n);
			for(int j=0;j<N;j++){
				y[i][j] = vec_fft1[j];
			}
		}
		// do the column based fft
		for(int i=0;i<N;i++){
			for(int j=0;j<M;j++){
				vec2[j] = y[j][i];
			}
			fft_radix2_iterative_c_1d<T, M>(vec2,vec_fft2,log2n);
			for(int j=0;j<M;j++){
				y[j][i] = vec_fft2[j];
			}
		}
	}

	template<class T, int M, int N>
	void ifft_radix2_iterative_c_2d(Complex<double> x[M][N],
						  Complex<double> y[M][N],
						  int log2n){
	    // conjugate the complex numbers
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				complex_conj<T>(x[i][j],x[i][j]);

	    // forward fft
	    fft_radix2_iterative_c_2d<T, M, N>( x, y, log2n );

	    // conjugate the complex numbers again
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				complex_conj<T>(y[i][j],y[i][j]);

	    // scale the numbers
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				y[i][j].real /= (T)(M*N); y[i][j].imag /= (T)(M*N);
			}
		}

	}

	template<class T, int LEN, int LEN2P, int LOG2N>
	void fft_bluestein_iterative_c_1d(Complex<T> x[LEN],
									  Complex<T> y[LEN]){

		T cos_table[LEN];
		T sin_table[LEN];

		// find the closest power of 2 value
		/*size_t m = 1;
		while (m / 2 <= LEN) {
			m *= 2;
		}
		int log2n = (int)std::log2(m);*/

		Complex<T> a[LEN2P], b[LEN2P];
		Complex<T> ao[LEN2P], aoo[LEN2P], bo[LEN2P];
		for(size_t i=0;i<LEN2P;i++){
			a[i].real = 0;
			a[i].imag = 0;
			b[i].real = 0;
			b[i].imag = 0;
		}

		// Trignometric tables
		for (size_t i = 0; i < LEN; i++) {
			unsigned long long temp = (unsigned long long)i * i;
			temp %= (unsigned long long)LEN * 2;
			float angle =  M_PI * temp / LEN;
			// Less accurate version if long long is unavailable: double angle = M_PI * i * i / n;
			cos_table[i] = (T)std::cos(angle);
			sin_table[i] = (T)std::sin(angle);
		}

		// Temporary vectors and preprocessing
		for (size_t i = 0; i < LEN; i++) {
			a[i].real =  x[i].real * cos_table[i] + x[i].imag * sin_table[i];
			a[i].imag = -x[i].real * sin_table[i] + x[i].imag * cos_table[i];
		}
		b[0].real = cos_table[0];
		b[0].imag = sin_table[0];
		for (size_t i = 1; i < LEN; i++) {
			b[i].real = cos_table[i];
			b[i].imag = sin_table[i];
			b[LEN2P - i].real = cos_table[i];
			b[LEN2P - i].imag = sin_table[i];
		}

		fft_radix2_iterative_c_1d<T, LEN2P>(a, ao, LOG2N);
		fft_radix2_iterative_c_1d<T, LEN2P>(b, bo, LOG2N);

		for (size_t i = 0; i < LEN2P; i++) {
			T temp = ao[i].real * bo[i].real - ao[i].imag * bo[i].imag;
			ao[i].imag = ao[i].imag * bo[i].real + ao[i].real * bo[i].imag;
			ao[i].real = temp;
		}
		ifft_radix2_iterative_c_1d<T, LEN2P>(ao, aoo, LOG2N);

		// Postprocessing
		for (size_t i = 0; i < LEN; i++) {
			y[i].real =  aoo[i].real * cos_table[i] + aoo[i].imag * sin_table[i];
			y[i].imag = -aoo[i].real * sin_table[i] + aoo[i].imag * cos_table[i];
		}

	}

	template<class T, int LEN, int LEN2P, int LOG2N>
	void ifft_bluestein_iterative_c_1d(Complex<T> x[LEN],
							Complex<T> y[LEN]){

	    // conjugate the complex numbers
		for(int i=0;i<LEN;i++)
			complex_conj(x[i],x[i]);

	    // forward fft
		fft_bluestein_iterative_c_1d<T, LEN, LEN2P, LOG2N>( x, y );

	    // conjugate the complex numbers again
		for(int i=0;i<LEN;i++)
			complex_conj(y[i],y[i]);

	    // scale the numbers
		for(int i=0;i<LEN;i++){
			y[i].real /= (T)LEN; y[i].imag /= (T)LEN;
		}
	}

	// Iterative FFT function to compute the DFT
	// of given coefficient vector
	template<class T, int M, int N, int M2P, int N2P, int L2NR, int L2NC>
	void fft_bluestein_iterative_c_2d(Complex<T> x[M][N],
					 	 Complex<T> y[M][N]){
		// clone input
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				y[i][j].real = x[i][j].real;
				y[i][j].imag = x[i][j].imag;
			}
		}

		Complex<T> vec1[N];
		Complex<T> vec_fft1[N];
		Complex<T> vec2[M];
		Complex<T> vec_fft2[M];
		// do the row based fft
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				vec1[j] = y[i][j];
			}
			fft_bluestein_iterative_c_1d<T, N, N2P, L2NC>(vec1,vec_fft1);
			for(int j=0;j<N;j++){
				y[i][j] = vec_fft1[j];
			}
		}
		// do the column based fft
		for(int i=0;i<N;i++){
			for(int j=0;j<M;j++){
				vec2[j] = y[j][i];
			}
			fft_bluestein_iterative_c_1d<T, M, M2P, L2NR>(vec2,vec_fft2);
			for(int j=0;j<M;j++){
				y[j][i] = vec_fft2[j];
			}
		}
	}

	template<class T, int M, int N, int M2P, int N2P, int L2NR, int L2NC>
	void ifft_bluestein_iterative_c_2d(Complex<T> x[M][N],
						  Complex<T> y[M][N]){
	    // conjugate the complex numbers
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				complex_conj<T>(x[i][j],x[i][j]);

	    // forward fft
		fft_bluestein_iterative_c_2d<T, M, N, M2P, N2P, L2NR, L2NC>( x, y );

	    // conjugate the complex numbers again
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				complex_conj<T>(y[i][j],y[i][j]);

	    // scale the numbers
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				y[i][j].real /= (M*N); y[i][j].imag /= (M*N);
			}
		}

	}

};

class FFT_CPP{
public:

	template<class T>
	T bitreverse(T num, size_t bitsnum)
	{
		size_t  NO_OF_BITS = bitsnum;
	    T reverse_num = 0;
	    size_t i;
	    for (i = 0; i < NO_OF_BITS; i++)
	    {
	        if((num & (1 << i)))
	           reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
	   }
	    return reverse_num;
	}

	// Cooley–Tukey FFT (in-place, divide-and-conquer)
	// Higher memory requirements and redundancy although more intuitive
	// Recursive function of FFT
	template<class T>
	void fft_radix2_recursive_cpp_1d(
			std::vector<std::complex<T>> &x,
			std::vector<std::complex<T>> &y)
	{
	    int n = x.size();

	    // if input contains just one element
	    if (n == 1)
	        y = std::vector<std::complex<T>>(1, x[0]);
	    else{
			// For storing n complex nth roots of unity
			std::vector<std::complex<T>> w(n);
			for (int i = 0; i < n/2; i++) {
				T alpha = 2 * M_PI * i / n;
				w[i] = std::complex<T>(cos(alpha), sin(alpha));
			}

			std::vector<std::complex<T>> A0(n / 2), A1(n / 2);
			for (int i = 0; i < n / 2; i++) {

				// even indexed coefficients
				A0[i] = x[i * 2];

				// odd indexed coefficients
				A1[i] = x[i * 2 + 1];
			}

			// Recursive call for even indexed coefficients
			std::vector<std::complex<T>> y0;
			fft_radix2_recursive_cpp_1d<T>(A0, y0);

			// Recursive call for odd indexed coefficients
			std::vector<std::complex<T>> y1;
			fft_radix2_recursive_cpp_1d<T>(A1, y1);

			// for storing values of y0, y1, y2, ..., yn-1.
			std::vector<std::complex<T>> y(n);

			for (int k = 0; k < n / 2; k++) {
				y[k] = y0[k] + w[k] * y1[k];
				y[k + n / 2] = y0[k] - w[k] * y1[k];
			}
	    }

	}

	// Iterative FFT function to compute the DFT
	// of given coefficient vector
	template<class T>
	void fft_radix2_iterative_cpp_1d(std::vector<std::complex<T>> x,
						   std::vector<std::complex<T>> &y,
						   size_t log2n)
	{
		size_t n = x.size();

	    // bit reversal of the given array
	    for (size_t i = 0; i < n; ++i) {
	    	size_t rev = bitreverse<size_t>(i, log2n);
	        y[i] = x[rev];
	    }

	    // j is iota
	    const std::complex<T> J(0,1);
	    for (size_t s = 1; s <= log2n; ++s) {
	    	size_t m = 1 << s; // 2 power s
	    	size_t m2 = m >> 1; // m2 = m/2 -1
	        std::complex<T> w(1,0);

	        // principle root of nth complex
	        // root of unity.
	        std::complex<T> wm = std::exp(J * (PI / m2));
	        for (size_t j = 0; j < m2; ++j) {
	            for (size_t k = j; k < n; k += m) {

	                // t = twiddle factor
	            	std::complex<T> t = w * y[k + m2];
	            	std::complex<T> u = y[k];

	                // similar calculating y[k]
	                y[k] = u + t;

	                // similar calculating y[k+n/2]
	                y[k + m2] = u - t;
	            }
	            w *= wm;
	        }
	    }
	}

	template<class T>
	void ifft_radix2_recursive_cpp_1d(
			std::vector<std::complex<T>> &x,
			std::vector<std::complex<T>> &y){

		y.resize(x.size());

	    // conjugate the complex numbers
		for(int i=0;i<(int)x.size();i++)
			x[i] = std::conj(x[i]);

	    // forward fft
		fft_radix2_recursive_cpp_1d<T>( x, y );

	    // conjugate the complex numbers again
		for(int i=0;i<(int)x.size();i++)
			y[i] = std::conj(y[i]);

	    // scale the numbers
		for(int i=0;i<(int)x.size();i++)
			y[i] /= y.size();

	}

	template<class T>
	void ifft_radix2_iterative_cpp_1d(std::vector<std::complex<T>> x,
							std::vector<std::complex<T>> &y,
							int log2num)
	{
	    // conjugate the complex numbers
		for(int i=0;i<(int)x.size();i++)
			x[i] = std::conj(x[i]);

	    // forward fft
	    fft_radix2_iterative_cpp_1d<T>( x, y, log2num );

	    // conjugate the complex numbers again
		for(int i=0;i<(int)x.size();i++)
			y[i] = std::conj(y[i]);

	    // scale the numbers
		for(int i=0;i<(int)x.size();i++)
			y[i] /= y.size();
	}
};

class FFT_FFTW3{
public:

	/*
	 * Real-to-Real Transform Kinds:
	FFTW_R2HC computes a real-input DFT with output in “halfcomplex” format, i.e. real and imaginary parts
	for a transform of size n stored as: r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1 (Logical N=n, inverse is FFTW_HC2R.)

	FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)

	FFTW_DHT computes a discrete Hartley transform. (Logical N=n, inverse is FFTW_DHT.)

	FFTW_REDFT00 computes an REDFT00 transform, i.e. a DCT-I. (Logical N=2*(n-1), inverse is FFTW_REDFT00.)

	FFTW_REDFT10 computes an REDFT10 transform, i.e. a DCT-II (sometimes called “the” DCT).
	(Logical N=2*n, inverse is FFTW_REDFT01.)

	FFTW_REDFT01 computes an REDFT01 transform, i.e. a DCT-III (sometimes called “the” IDCT, being the inverse of DCT-II).
	 (Logical N=2*n, inverse is FFTW_REDFT=10.)

	FFTW_REDFT11 computes an REDFT11 transform, i.e. a DCT-IV. (Logical N=2*n, inverse is FFTW_REDFT11.)

	FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)

	FFTW_RODFT10 computes an RODFT10 transform, i.e. a DST-II. (Logical N=2*n, inverse is FFTW_RODFT01.)

	FFTW_RODFT01 computes an RODFT01 transform, i.e. a DST-III. (Logical N=2*n, inverse is FFTW_RODFT=10.)

	FFTW_RODFT11 computes an RODFT11 transform, i.e. a DST-IV. (Logical N=2*n, inverse is FFTW_RODFT11.)
	*/
	/*
	 * Planner Flags
	FFTW_ESTIMATE specifies that, instead of actual measurements of different algorithms, a simple heuristic is used
	to pick a (probably sub-optimal) plan quickly. With this flag, the input/output arrays are not overwritten during planning.

	FFTW_MEASURE tells FFTW to find an optimized plan by actually computing several FFTs and measuring their execution time.
	Depending on your machine, this can take some time (often a few seconds). FFTW_MEASURE is the default planning option.

	FFTW_PATIENT is like FFTW_MEASURE, but considers a wider range of algorithms and often produces a “more optimal” plan
	(especially for large transforms), but at the expense of several times longer planning time (especially for large transforms).

	FFTW_EXHAUSTIVE is like FFTW_PATIENT, but considers an even wider range of algorithms, including many that we think are unlikely
	to be fast, to produce the most optimal plan but with a substantially increased planning time.

	FFTW_WISDOM_ONLY is a special planning mode in which the plan is only created if wisdom is available for the given problem, and
	otherwise a NULL plan is returned. This can be combined with other flags, e.g. ‘FFTW_WISDOM_ONLY | FFTW_PATIENT’ creates a plan
	only if wisdom is available that was created in FFTW_PATIENT or FFTW_EXHAUSTIVE mode. The FFTW_WISDOM_ONLY flag is intended for
	users who need to detect whether wisdom is available; for example, if wisdom is not available one may wish to allocate new arrays
	for planning so that user data is not overwritten.
	 */
	template<class T, int LEN>
	void fftw3_rfft_1d(T IN[LEN], T OUT[LEN]){

		fftw_plan plan = fftw_plan_r2r_1d(LEN, IN, OUT, FFTW_RODFT10, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		fftw_cleanup();

	}
	/*
	 * The aliases FFTW_FORWARD and FFTW_BACKWARD are provided, where FFTW_FORWARD stands for -1
	 * Note also that we use the standard “in-order” output ordering—the k-th output corresponds
	 * to the frequency k/n (or k/T, where T is your total sampling period). For those who like
	 * to think in terms of positive and negative frequencies, this means that the positive
	 * frequencies are stored in the first half of the output and the negative frequencies are
	 * stored in backwards order in the second half of the output. (The frequency -k/n is the same
	 *  as the frequency (n-k)/n.)
	*/
	template<class T, int LEN>
	void fftw3_cfft_1d(Complex<T> IN[LEN], Complex<T> OUT[LEN]){
		fftw_complex in[LEN];
		fftw_complex out[LEN];
	//	in = (fftw_complex *) malloc(sizeof(fftw_complex)*LEN);
	//	out = (fftw_complex *) malloc(sizeof(fftw_complex)*LEN);
		for(int i=0;i<LEN;i++){
			in[i][0] = IN[i].real;
			in[i][1] = IN[i].imag;
			//std::cout << in[i][0] << "," << in[i][1] << "\n";
		}
		fftw_plan plan = fftw_plan_dft_1d(LEN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<LEN;i++){
			OUT[i].real = out[i][0];
			OUT[i].imag = out[i][1];
			//std::cout << out[i][0] << "," << out[i][1] << "\n";
		}
		fftw_destroy_plan(plan);
		fftw_cleanup();

	}

	template<class T, int LEN>
	void fftw3_cifft_1d(Complex<T> IN[LEN], Complex<T> OUT[LEN]){
		fftw_complex in[LEN];
		fftw_complex out[LEN];
	//	in = (fftw_complex *) malloc(sizeof(fftw_complex)*LEN);
	//	out = (fftw_complex *) malloc(sizeof(fftw_complex)*LEN);
		for(int i=0;i<LEN;i++){
			in[i][0] = IN[i].real;
			in[i][1] = IN[i].imag;
			//std::cout << in[i][0] << "," << in[i][1] << "\n";
		}
		fftw_plan plan = fftw_plan_dft_1d(LEN, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<LEN;i++){
			OUT[i].real = out[i][0];
			OUT[i].imag = out[i][1];
			//std::cout << out[i][0] << "," << out[i][1] << "\n";
		}
		fftw_destroy_plan(plan);
		fftw_cleanup();

	}

	template<class T, int M, int N>
	void fftw3_cfft_2d(Complex<T> IN[M][N], Complex<T> OUT[M][N]){
		fftw_complex *in;
		fftw_complex *out;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				in[i*M+j][0] = IN[i][j].real;
				in[i*M+j][1] = IN[i][j].imag;
			}
		}
		fftw_plan plan = fftw_plan_dft_2d(M, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				OUT[i][j].real = out[i*M+j][0];
				OUT[i][j].imag = out[i*M+j][1];
			}
		}
		fftw_destroy_plan(plan);
		fftw_cleanup();
	}
	template<class T, int M, int N>
	void fftw3_cfft_2d(T IN, T OUT){
		fftw_complex *in;
		fftw_complex *out;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				in[i*M+j][0] = IN(i,j).real();
				in[i*M+j][1] = IN(i,j).imag();
			}
		}
		fftw_plan_with_nthreads(omp_get_max_threads());
		fftw_plan plan = fftw_plan_dft_2d(M, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				OUT(i,j).real(out[i*M+j][0]);
				OUT(i,j).imag(out[i*M+j][1]);
			}
		}
		fftw_destroy_plan(plan);
		fftw_cleanup_threads();
		fftw_cleanup();
	}

	template<class T, int M, int N>
	void fftw3_cifft_2d(Complex<T> IN[M][N], Complex<T> OUT[M][N]){
		fftw_complex *in;
		fftw_complex *out;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				in[i*M+j][0] = IN[i][j].real;
				in[i*M+j][1] = IN[i][j].imag;
			}
		}
		fftw_plan plan = fftw_plan_dft_2d(M, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				OUT[i][j].real = out[i*M+j][0]/(M*N);
				OUT[i][j].imag = out[i*M+j][1]/(M*N);
			}
		}
		fftw_destroy_plan(plan);
		fftw_cleanup();
	}
	template<class T, int M, int N>
	void fftw3_cifft_2d(T IN, T OUT){
		fftw_complex *in;
		fftw_complex *out;
		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M * N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				in[i*M+j][0] = IN(i,j).real();
				in[i*M+j][1] = IN(i,j).imag();
			}
		}
		fftw_plan_with_nthreads(omp_get_max_threads());
		fftw_plan plan = fftw_plan_dft_2d(M, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				OUT(i,j).real(out[i*M+j][0]/(M*N));
				OUT(i,j).imag(out[i*M+j][1]/(M*N));
			}
		}
		fftw_destroy_plan(plan);
		fftw_cleanup_threads();
		fftw_cleanup();
	}


};

#endif /* SRC_FFT_HPP_ */
