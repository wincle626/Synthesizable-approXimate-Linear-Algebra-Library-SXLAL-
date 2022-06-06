/*
 * fft.cpp
 *
 *  Created on: 24 Jul 2019
 *      Author: yunwu
 */
/*
#include "fft.hpp"
void complex_add(Complex<double> a, Complex<double> b, Complex<double> &c){
	c.real = a.real + b.real;
	c.imag = a.imag + b.imag;
}

void complex_sub(Complex<double> a, Complex<double> b, Complex<double> &c){
	c.real = a.real - b.real;
	c.imag = a.imag - b.imag;
}

void complex_mul(Complex<double> a, Complex<double> b, Complex<double> &c){
	c.real = a.real * b.real - a.imag * b.imag;
	c.imag = a.real * b.imag + a.imag * b.real;
}

void complex_div(Complex<double> a, Complex<double> b, Complex<double> &c){
	double x = a.real*b.real + a.imag*b.imag;
	double y = a.imag*b.real - a.real*b.imag;
	double z = b.real*b.real + b.imag*b.imag;
	c.real = x/z;
	c.imag = y/z;
}

void complex_conj(Complex<double>a, Complex<double> &b){
	b.real = a.real;
	b.imag = a.imag * -1;
}

void complex_abs(Complex<double> a, double &b){
	b = std::sqrt(a.real * a.real + a.imag * a.imag);
}

void complex_abs2(Complex<double> a, double &b){
	b = a.real * a.real + a.imag * a.imag;
}

void exp2complex(double phi, Complex<double> &a){
	a.real = std::cos(phi);
	a.imag = std::sin(phi);
}

template<class T>
void subvec(T *vec, T *subvec,
			int size, int subsize,
			int start, int step){
	for(int i=start,j=0;i<size&&j<subsize;i+=step,j++)
		subvec[j] = vec[i];
}


// * Real-to-Real Transform Kinds:
//FFTW_R2HC computes a real-input DFT with output in “halfcomplex” format, i.e. real and imaginary parts
//for a transform of size n stored as: r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1 (Logical N=n, inverse is FFTW_HC2R.)
//
//FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)
//
//FFTW_DHT computes a discrete Hartley transform. (Logical N=n, inverse is FFTW_DHT.)
//
//FFTW_REDFT00 computes an REDFT00 transform, i.e. a DCT-I. (Logical N=2*(n-1), inverse is FFTW_REDFT00.)
//
//FFTW_REDFT10 computes an REDFT10 transform, i.e. a DCT-II (sometimes called “the” DCT).
//(Logical N=2*n, inverse is FFTW_REDFT01.)
//
//FFTW_REDFT01 computes an REDFT01 transform, i.e. a DCT-III (sometimes called “the” IDCT, being the inverse of DCT-II).
// (Logical N=2*n, inverse is FFTW_REDFT=10.)
//
//FFTW_REDFT11 computes an REDFT11 transform, i.e. a DCT-IV. (Logical N=2*n, inverse is FFTW_REDFT11.)
//
//FFTW_RODFT00 computes an RODFT00 transform, i.e. a DST-I. (Logical N=2*(n+1), inverse is FFTW_RODFT00.)
//
//FFTW_RODFT10 computes an RODFT10 transform, i.e. a DST-II. (Logical N=2*n, inverse is FFTW_RODFT01.)
//
//FFTW_RODFT01 computes an RODFT01 transform, i.e. a DST-III. (Logical N=2*n, inverse is FFTW_RODFT=10.)
//
//FFTW_RODFT11 computes an RODFT11 transform, i.e. a DST-IV. (Logical N=2*n, inverse is FFTW_RODFT11.)
//
// * Planner Flags
//FFTW_ESTIMATE specifies that, instead of actual measurements of different algorithms, a simple heuristic is used
//to pick a (probably sub-optimal) plan quickly. With this flag, the input/output arrays are not overwritten during planning.
//
//FFTW_MEASURE tells FFTW to find an optimized plan by actually computing several FFTs and measuring their execution time.
//Depending on your machine, this can take some time (often a few seconds). FFTW_MEASURE is the default planning option.
//
//FFTW_PATIENT is like FFTW_MEASURE, but considers a wider range of algorithms and often produces a “more optimal” plan
//(especially for large transforms), but at the expense of several times longer planning time (especially for large transforms).
//
//FFTW_EXHAUSTIVE is like FFTW_PATIENT, but considers an even wider range of algorithms, including many that we think are unlikely
//to be fast, to produce the most optimal plan but with a substantially increased planning time.
//
//FFTW_WISDOM_ONLY is a special planning mode in which the plan is only created if wisdom is available for the given problem, and
//otherwise a NULL plan is returned. This can be combined with other flags, e.g. ‘FFTW_WISDOM_ONLY | FFTW_PATIENT’ creates a plan
//only if wisdom is available that was created in FFTW_PATIENT or FFTW_EXHAUSTIVE mode. The FFTW_WISDOM_ONLY flag is intended for
//users who need to detect whether wisdom is available; for example, if wisdom is not available one may wish to allocate new arrays
//for planning so that user data is not overwritten.

void fftw3_rfft_1d(double IN[DIAG], double OUT[DIAG]){

	fftw_plan plan = fftw_plan_r2r_1d(DIAG, IN, OUT, FFTW_RODFT10, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	fftw_cleanup();

}

// bitsnum is normally the log2 of the size
unsigned int bitreverse1(unsigned int num, unsigned int bitsnum)
{
    unsigned int  NO_OF_BITS = bitsnum;
    unsigned int reverse_num = 0;
    unsigned int i;
    for (i = 0; i < NO_OF_BITS; i++)
    {
        if((num & (1 << i)))
           reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
   }
    return reverse_num;
}
unsigned int bitreverse2(unsigned int num, unsigned int bitsnum){
	unsigned int count = bitsnum;
	unsigned int reverse_num = num;
	num >>= 1;
	while(num)
	{
		reverse_num <<= 1;
		reverse_num |= num & 1;
		num >>= 1;
		count--;
	}
	reverse_num <<= count;
	return reverse_num;
}
void swapvecelement(Complex<double> Vec[DIAG], int i, int j){
	Complex<double> tmp;
	tmp.real = Vec[i].real;
	tmp.imag = Vec[i].imag;
	Vec[i].real = Vec[j].real;
	Vec[i].imag = Vec[j].imag;
	Vec[j].real = tmp.real;
	Vec[j].imag = tmp.imag;
}
void bitreverse_reorder_1d(Complex<double> Vec[DIAG]){
	for(int i=0;i<DIAG/2;i++){
		unsigned int bitsnum = (unsigned int) std::log2(DIAG);
		unsigned int bitrevind = bitreverse1(i, bitsnum);
		//std::cout << "swap " << i << " and " << bitrevind << std::endl;
		swapvecelement(Vec, i, bitrevind);
	}
}

// * The aliases FFTW_FORWARD and FFTW_BACKWARD are provided, where FFTW_FORWARD stands for -1
// * Note also that we use the standard “in-order” output ordering—the k-th output corresponds
// * to the frequency k/n (or k/T, where T is your total sampling period). For those who like
// * to think in terms of positive and negative frequencies, this means that the positive
// * frequencies are stored in the first half of the output and the negative frequencies are
// * stored in backwards order in the second half of the output. (The frequency -k/n is the same
// *  as the frequency (n-k)/n.)
void fftw3_cfft_1d(Complex<double> IN[DIAG], Complex<double> OUT[DIAG]){
	fftw_complex in[DIAG];
	fftw_complex out[DIAG];
//	in = (fftw_complex *) malloc(sizeof(fftw_complex)*DIAG);
//	out = (fftw_complex *) malloc(sizeof(fftw_complex)*DIAG);
	for(int i=0;i<DIAG;i++){
		in[i][0] = IN[i].real;
		in[i][1] = IN[i].imag;
		//std::cout << in[i][0] << "," << in[i][1] << "\n";
	}
	fftw_plan plan = fftw_plan_dft_1d(DIAG, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	for(int i=0;i<DIAG;i++){
		OUT[i].real = out[i][0];
		OUT[i].imag = out[i][1];
		//std::cout << out[i][0] << "," << out[i][1] << "\n";
	}
	fftw_destroy_plan(plan);
	fftw_cleanup();

}
void fftw3_cifft_1d(Complex<double> IN[DIAG], Complex<double> OUT[DIAG]){
	fftw_complex in[DIAG];
	fftw_complex out[DIAG];
//	in = (fftw_complex *) malloc(sizeof(fftw_complex)*DIAG);
//	out = (fftw_complex *) malloc(sizeof(fftw_complex)*DIAG);
	for(int i=0;i<DIAG;i++){
		in[i][0] = IN[i].real;
		in[i][1] = IN[i].imag;
		//std::cout << in[i][0] << "," << in[i][1] << "\n";
	}
	fftw_plan plan = fftw_plan_dft_1d(DIAG, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	for(int i=0;i<DIAG;i++){
		OUT[i].real = out[i][0];
		OUT[i].imag = out[i][1];
		//std::cout << out[i][0] << "," << out[i][1] << "\n";
	}
	fftw_destroy_plan(plan);
	fftw_cleanup();

}
void fftw3_cfft_2d(Complex<double> IN[DIAG][DIAG], Complex<double> OUT[DIAG][DIAG]){
	fftw_complex *in;
	fftw_complex *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DIAG * DIAG);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DIAG * DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			in[i*DIAG+j][0] = IN[i][j].real;
			in[i*DIAG+j][1] = IN[i][j].imag;
		}
	}
	fftw_plan plan = fftw_plan_dft_2d(DIAG, DIAG, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			OUT[i][j].real = out[i*DIAG+j][0];
			OUT[i][j].imag = out[i*DIAG+j][1];
		}
	}
	fftw_destroy_plan(plan);
	fftw_cleanup();
}
void fftw3_cifft_2d(Complex<double> IN[DIAG][DIAG], Complex<double> OUT[DIAG][DIAG]){
	fftw_complex *in;
	fftw_complex *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DIAG * DIAG);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DIAG * DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			in[i*DIAG+j][0] = IN[i][j].real;
			in[i*DIAG+j][1] = IN[i][j].imag;
		}
	}
	fftw_plan plan = fftw_plan_dft_2d(DIAG, DIAG, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			OUT[i][j].real = out[i*DIAG+j][0]/(DIAG*DIAG);
			OUT[i][j].imag = out[i*DIAG+j][1]/(DIAG*DIAG);
		}
	}
	fftw_destroy_plan(plan);
	fftw_cleanup();
}



// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
// Recursive function of FFT
std::vector<std::complex<double>> fft_radix2_recursive_cpp_1d(std::vector<std::complex<double>> &x)
{
    int n = x.size();

    // if input contains just one element
    if (n == 1)
        return std::vector<std::complex<double>>(1, x[0]);

    // For storing n complex nth roots of unity
    std::vector<std::complex<double>> w(n);
    for (int i = 0; i < n/2; i++) {
        double alpha = 2 * M_PI * i / n;
        w[i] = std::complex<double>(cos(alpha), sin(alpha));
    }

    std::vector<std::complex<double>> A0(n / 2), A1(n / 2);
    for (int i = 0; i < n / 2; i++) {

        // even indexed coefficients
        A0[i] = x[i * 2];

        // odd indexed coefficients
        A1[i] = x[i * 2 + 1];
    }

    // Recursive call for even indexed coefficients
    std::vector<std::complex<double>> y0 = fft_radix2_recursive_cpp_1d(A0);

    // Recursive call for odd indexed coefficients
    std::vector<std::complex<double>> y1 = fft_radix2_recursive_cpp_1d(A1);

    // for storing values of y0, y1, y2, ..., yn-1.
    std::vector<std::complex<double>> y(n);

    for (int k = 0; k < n / 2; k++) {
        y[k] = y0[k] + w[k] * y1[k];
        y[k + n / 2] = y0[k] - w[k] * y1[k];
    }
    return y;
}

// Iterative FFT function to compute the DFT
// of given coefficient vector
void fft_radix2_iterative_cpp_1d(std::vector<std::complex<double>> x,
					   std::vector<std::complex<double>> &y,
					   int log2n)
{
    int n = x.size();

    // bit reversal of the given array
    for (int i = 0; i < n; ++i) {
        int rev = bitreverse1(i, log2n);
        y[i] = x[rev];
    }

    // j is iota
    const std::complex<double> J(0,1);
    for (int s = 1; s <= log2n; ++s) {
        int m = 1 << s; // 2 power s
        int m2 = m >> 1; // m2 = m/2 -1
        std::complex<double> w(1,0);

        // principle root of nth complex
        // root of unity.
        std::complex<double> wm = std::exp(J * (PI / m2));
        for (int j = 0; j < m2; ++j) {
            for (int k = j; k < n; k += m) {

                // t = twiddle factor
            	std::complex<double> t = w * y[k + m2];
            	std::complex<double> u = y[k];

                // similar calculating y[k]
                y[k] = u + t;

                // similar calculating y[k+n/2]
                y[k + m2] = u - t;
            }
            w *= wm;
        }
    }
}

std::vector<std::complex<double>> ifft_radix2_recursive_cpp_1d(std::vector<std::complex<double>> &x){

	std::vector<std::complex<double>> y;
	y.resize(x.size());

    // conjugate the complex numbers
	for(int i=0;i<(int)x.size();i++)
		x[i] = std::conj(x[i]);

    // forward fft
	y = fft_radix2_recursive_cpp_1d( x );

    // conjugate the complex numbers again
	for(int i=0;i<(int)x.size();i++)
		y[i] = std::conj(y[i]);

    // scale the numbers
	for(int i=0;i<(int)x.size();i++)
		y[i] /= y.size();

	return y;
}

void ifft_radix2_iterative_cpp_1d(std::vector<std::complex<double>> x,
						std::vector<std::complex<double>> &y,
						int log2num)
{
    // conjugate the complex numbers
	for(int i=0;i<(int)x.size();i++)
		x[i] = std::conj(x[i]);

    // forward fft
    fft_radix2_iterative_cpp_1d( x, y, log2num );

    // conjugate the complex numbers again
	for(int i=0;i<(int)x.size();i++)
		y[i] = std::conj(y[i]);

    // scale the numbers
	for(int i=0;i<(int)x.size();i++)
		y[i] /= y.size();
}

// Iterative FFT function to compute the DFT
// of given coefficient vector
void fft_radix2_iterative_c_1d(Complex<double> x[DIAG],
				 	 Complex<double> y[DIAG],
					 int log2n){

    int n = (int) std::pow(2,log2n);

    // bit reversal of the given array
    for (int i = 0; i < n; ++i) {
        int rev = bitreverse1(i, log2n);
        y[i].real = x[rev].real;
        y[i].imag = x[rev].imag;
    }

    // j is iota
    const std::complex<double> J(0,1);
    cmplx<double> tmp;
    for (int s = 1; s <= log2n; ++s) {
        int m = 1 << s; // 2 power s
        int m2 = m >> 1; // m2 = m/2 -1
        std::complex<double> w1(1,0);
        cmplx<double> w;
        w.real = 1; w.imag = 0;

        // principle root of nth complex
        // root of unity.
        cmplx<double> wm;
        wm.real = std::cos(PI / m2);
        wm.imag = std::sin(PI / m2);
        for (int j = 0; j < m2; ++j) {
            for (int k = j; k < n; k += m) {

                // t = twiddle factor
            	cmplx<double> t;
            	complex_mul(w, y[k+m2], t);
            	cmplx<double> u;
            	u.real = y[k].real; u.imag = y[k].imag;

                // similar calculating y[k]
            	complex_add(u, t, tmp);
            	y[k].real = tmp.real;
            	y[k].imag = tmp.imag;

                // similar calculating y[k+n/2]
            	complex_sub(u, t, tmp);
            	y[k + m2].real = tmp.real;
            	y[k + m2].imag = tmp.imag;
            }
            complex_mul(w, wm, tmp);
            w.real = tmp.real; w.imag = tmp.imag;
        }
    }
}

void ifft_radix2_iterative_c_1d(Complex<double> x[DIAG],
					  Complex<double> y[DIAG],
					  int log2n)
{
    int n = (int) std::pow(2,log2n);

    // conjugate the complex numbers
	for(int i=0;i<n;i++)
		complex_conj(x[i],x[i]);

    // forward fft
    fft_radix2_iterative_c_1d( x, y, log2n );

    // conjugate the complex numbers again
	for(int i=0;i<n;i++)
		complex_conj(y[i],y[i]);

    // scale the numbers
	for(int i=0;i<n;i++){
		y[i].real /= n; y[i].imag /= n;
	}
}

void fft_bluestein_iterative_c_1d(Complex<double> x[DIAG],
						Complex<double> y[DIAG]){

	double cos_table[DIAG];
	double sin_table[DIAG];

	// find the closest power of 2 value
	size_t m = 1;
	while (m / 2 <= DIAG) {
		m *= 2;
	}
	int log2n = (int)std::log2(m);

	Complex<double> a[m], b[m];
	Complex<double> ao[m], aoo[m], bo[m];
	for(size_t i=0;i<m;i++){
		a[i].real = 0;
		a[i].imag = 0;
		b[i].real = 0;
		b[i].imag = 0;
	}

	// Trignometric tables
	for (size_t i = 0; i < DIAG; i++) {
		unsigned long long temp = (unsigned long long)i * i;
		temp %= (unsigned long long)DIAG * 2;
		double angle = M_PI * temp / DIAG;
		// Less accurate version if long long is unavailable: double angle = M_PI * i * i / n;
		cos_table[i] = cos(angle);
		sin_table[i] = sin(angle);
	}

	// Temporary vectors and preprocessing
	for (size_t i = 0; i < DIAG; i++) {
		a[i].real =  x[i].real * cos_table[i] + x[i].imag * sin_table[i];
		a[i].imag = -x[i].real * sin_table[i] + x[i].imag * cos_table[i];
	}
	b[0].real = cos_table[0];
	b[0].imag = sin_table[0];
	for (size_t i = 1; i < DIAG; i++) {
		b[i].real = cos_table[i];
		b[i].imag = sin_table[i];
		b[m - i].real = cos_table[i];
		b[m - i].imag = sin_table[i];
	}

	fft_radix2_iterative_c_1d(a, ao, log2n);
	fft_radix2_iterative_c_1d(b, bo, log2n);

	for (size_t i = 0; i < m; i++) {
		double temp = ao[i].real * bo[i].real - ao[i].imag * bo[i].imag;
		ao[i].imag = ao[i].imag * bo[i].real + ao[i].real * bo[i].imag;
		ao[i].real = temp;
	}
	ifft_radix2_iterative_c_1d(ao, aoo, log2n);

	// Postprocessing
	for (size_t i = 0; i < DIAG; i++) {
		y[i].real =  aoo[i].real * cos_table[i] + aoo[i].imag * sin_table[i];
		y[i].imag = -aoo[i].real * sin_table[i] + aoo[i].imag * cos_table[i];
	}

}

void ifft_bluestein_iterative_c_1d(Complex<double> x[DIAG],
						Complex<double> y[DIAG]){

    int n = DIAG;

    // conjugate the complex numbers
	for(int i=0;i<n;i++)
		complex_conj(x[i],x[i]);

    // forward fft
	fft_bluestein_iterative_c_1d( x, y );

    // conjugate the complex numbers again
	for(int i=0;i<n;i++)
		complex_conj(y[i],y[i]);

    // scale the numbers
	for(int i=0;i<n;i++){
		y[i].real /= n; y[i].imag /= n;
	}
}

// Iterative FFT function to compute the DFT
// of given coefficient vector
void fft_radix2_iterative_c_2d(Complex<double> x[DIAG][DIAG],
				 	 Complex<double> y[DIAG][DIAG],
					 int log2n){
	// clone input
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			y[i][j].real = x[i][j].real;
			y[i][j].imag = x[i][j].imag;
		}
	}

	Complex<double> vec[DIAG];
	Complex<double> vec_fft[DIAG];
	// do the row based fft
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			vec[j] = y[i][j];
		}
		fft_radix2_iterative_c_1d(vec,vec_fft,log2n);
		for(int j=0;j<DIAG;j++){
			y[i][j] = vec_fft[j];
		}
	}
	// do the column based fft
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			vec[j] = y[j][i];
		}
		fft_radix2_iterative_c_1d(vec,vec_fft,log2n);
		for(int j=0;j<DIAG;j++){
			y[j][i] = vec_fft[j];
		}
	}
}

void ifft_radix2_iterative_c_2d(Complex<double> x[DIAG][DIAG],
					  Complex<double> y[DIAG][DIAG],
					  int log2n){
    // conjugate the complex numbers
	for(int i=0;i<DIAG;i++)
		for(int j=0;j<DIAG;j++)
			complex_conj(x[i][j],x[i][j]);

    // forward fft
    fft_radix2_iterative_c_2d( x, y, log2n );

    // conjugate the complex numbers again
	for(int i=0;i<DIAG;i++)
		for(int j=0;j<DIAG;j++)
			complex_conj(y[i][j],y[i][j]);

    // scale the numbers
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			y[i][j].real /= (DIAG*DIAG); y[i][j].imag /= (DIAG*DIAG);
		}
	}

}

// Iterative FFT function to compute the DFT
// of given coefficient vector
void fft_bluestein_iterative_c_2d(Complex<double> x[DIAG][DIAG],
				 	 Complex<double> y[DIAG][DIAG]){
	// clone input
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			y[i][j].real = x[i][j].real;
			y[i][j].imag = x[i][j].imag;
		}
	}

	Complex<double> vec[DIAG];
	Complex<double> vec_fft[DIAG];
	// do the row based fft
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			vec[j] = y[i][j];
		}
		fft_bluestein_iterative_c_1d(vec,vec_fft);
		for(int j=0;j<DIAG;j++){
			y[i][j] = vec_fft[j];
		}
	}
	// do the column based fft
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			vec[j] = y[j][i];
		}
		fft_bluestein_iterative_c_1d(vec,vec_fft);
		for(int j=0;j<DIAG;j++){
			y[j][i] = vec_fft[j];
		}
	}
}

void ifft_bluestein_iterative_c_2d(Complex<double> x[DIAG][DIAG],
					  Complex<double> y[DIAG][DIAG]){
    // conjugate the complex numbers
	for(int i=0;i<DIAG;i++)
		for(int j=0;j<DIAG;j++)
			complex_conj(x[i][j],x[i][j]);

    // forward fft
    fft_bluestein_iterative_c_2d( x, y );

    // conjugate the complex numbers again
	for(int i=0;i<DIAG;i++)
		for(int j=0;j<DIAG;j++)
			complex_conj(y[i][j],y[i][j]);

    // scale the numbers
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			y[i][j].real /= (DIAG*DIAG); y[i][j].imag /= (DIAG*DIAG);
		}
	}

}
*/
