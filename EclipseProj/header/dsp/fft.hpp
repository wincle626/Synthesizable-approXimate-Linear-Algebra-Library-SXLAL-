/*
 * fft.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_DSP_FFT_HPP_
#define HEADER_DSP_FFT_HPP_

#include "mathfunc.hpp"

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
void FFT_radix2_iterative_c_1d(Complex<T> x[LEN],
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
void iFFT_radix2_iterative_c_1d(Complex<T> x[LEN],
					  Complex<T> y[LEN],
					  int log2n)
{
	// conjugate the complex numbers
	for(int i=0;i<LEN;i++)
		complex_conj<T>(x[i],x[i]);

	// forward fft
	FFT_radix2_iterative_c_1d<T, LEN>( x, y, log2n );

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
void FFT_radix2_iterative_c_2d(Complex<double> x[M][N],
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
		FFT_radix2_iterative_c_1d<T, N>(vec1,vec_fft1,log2n);
		for(int j=0;j<N;j++){
			y[i][j] = vec_fft1[j];
		}
	}
	// do the column based fft
	for(int i=0;i<N;i++){
		for(int j=0;j<M;j++){
			vec2[j] = y[j][i];
		}
		FFT_radix2_iterative_c_1d<T, M>(vec2,vec_fft2,log2n);
		for(int j=0;j<M;j++){
			y[j][i] = vec_fft2[j];
		}
	}
}

template<class T, int M, int N>
void iFFT_radix2_iterative_c_2d(Complex<double> x[M][N],
					  Complex<double> y[M][N],
					  int log2n){
	// conjugate the complex numbers
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			complex_conj<T>(x[i][j],x[i][j]);

	// forward fft
	FFT_radix2_iterative_c_2d<T, M, N>( x, y, log2n );

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
void FFT_bluestein_iterative_c_1d(Complex<T> x[LEN],
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
	for(int i=0;i<LEN2P;i++){
		a[i].real = 0;
		a[i].imag = 0;
		b[i].real = 0;
		b[i].imag = 0;
	}

	// Trignometric tables
	for (int i = 0; i < LEN; i++) {
		unsigned long long temp = (unsigned long long)i * i;
		temp %= (unsigned long long)LEN * 2;
		float angle =  M_PI * temp / LEN;
		// Less accurate version if long long is unavailable: double angle = M_PI * i * i / n;
		cos_table[i] = (T)std::cos(angle);
		sin_table[i] = (T)std::sin(angle);
	}

	// Temporary vectors and preprocessing
	for (int i = 0; i < LEN; i++) {
		a[i].real =  x[i].real * cos_table[i] + x[i].imag * sin_table[i];
		a[i].imag = -x[i].real * sin_table[i] + x[i].imag * cos_table[i];
	}
	b[0].real = cos_table[0];
	b[0].imag = sin_table[0];
	for (int i = 1; i < LEN; i++) {
		b[i].real = cos_table[i];
		b[i].imag = sin_table[i];
		b[LEN2P - i].real = cos_table[i];
		b[LEN2P - i].imag = sin_table[i];
	}

	FFT_radix2_iterative_c_1d<T, LEN2P>(a, ao, LOG2N);
	FFT_radix2_iterative_c_1d<T, LEN2P>(b, bo, LOG2N);

	for (int i = 0; i < LEN2P; i++) {
		T temp = ao[i].real * bo[i].real - ao[i].imag * bo[i].imag;
		ao[i].imag = ao[i].imag * bo[i].real + ao[i].real * bo[i].imag;
		ao[i].real = temp;
	}
	iFFT_radix2_iterative_c_1d<T, LEN2P>(ao, aoo, LOG2N);

	// Postprocessing
	for (int i = 0; i < LEN; i++) {
		y[i].real =  aoo[i].real * cos_table[i] + aoo[i].imag * sin_table[i];
		y[i].imag = -aoo[i].real * sin_table[i] + aoo[i].imag * cos_table[i];
	}

}

template<class T, int LEN, int LEN2P, int LOG2N>
void iFFT_bluestein_iterative_c_1d(Complex<T> x[LEN],
						Complex<T> y[LEN]){

	// conjugate the complex numbers
	for(int i=0;i<LEN;i++)
		complex_conj(x[i],x[i]);

	// forward fft
	FFT_bluestein_iterative_c_1d<T, LEN, LEN2P, LOG2N>( x, y );

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
void FFT_bluestein_iterative_c_2d(Complex<T> x[M][N],
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
		FFT_bluestein_iterative_c_1d<T, N, N2P, L2NC>(vec1,vec_fft1);
		for(int j=0;j<N;j++){
			y[i][j] = vec_fft1[j];
		}
	}
	// do the column based fft
	for(int i=0;i<N;i++){
		for(int j=0;j<M;j++){
			vec2[j] = y[j][i];
		}
		FFT_bluestein_iterative_c_1d<T, M, M2P, L2NR>(vec2,vec_fft2);
		for(int j=0;j<M;j++){
			y[j][i] = vec_fft2[j];
		}
	}
}

template<class T, int M, int N, int M2P, int N2P, int L2NR, int L2NC>
void iFFT_bluestein_iterative_c_2d(Complex<T> x[M][N],
					  Complex<T> y[M][N]){
	// conjugate the complex numbers
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			complex_conj<T>(x[i][j],x[i][j]);

	// forward fft
	FFT_bluestein_iterative_c_2d<T, M, N, M2P, N2P, L2NR, L2NC>( x, y );

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



#endif /* HEADER_DSP_FFT_HPP_ */
