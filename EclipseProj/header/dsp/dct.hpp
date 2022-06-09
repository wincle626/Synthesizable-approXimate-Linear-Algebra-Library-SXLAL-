/*
 * dct.hpp discrete cosine transform (DCT)
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_DSP_DCT_HPP_
#define HEADER_DSP_DCT_HPP_

#include <cmath>

#define pi 3.14159265359

// DCT-1D Type-I
template<class T, int N>
void DCT_1D_TYPE1(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/(N-1));
	T b = (T) pi/(N-1);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		T kd1 = (i==1||i==N)?1:0;
		T c = (T) 1/std::sqrt(1+kd1);
		for(j=1;j<=N;j++){
			T kd2 = (j==1||j==N)?1:0;
			T d = (T) 1/std::sqrt(1+kd2);
			T e = (T) std::cos(b*(j-1)*(i-1));
			sum += x[j-1]*c*d*e;
		}
		dctx[i-1] = a*sum;
	}

}

// DCT-1D Type-II
template<class T, int N>
void DCT_1D_TYPE2(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*2);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		T kd1 = (i==1)?1:0;
		T c = (T) 1/std::sqrt(1+kd1);
		for(j=1;j<=N;j++){
			T d = (T) std::cos(b*(2*j-1)*(i-1));
			sum += x[j-1]*c*d;
		}
		dctx[i-1] = a*sum;
	}

}

// DCT-1D Type-III
template<class T, int N>
void DCT_1D_TYPE3(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*2);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		for(j=1;j<=N;j++){
			T kd2 = (j==1)?1:0;
			T c = (T) 1/std::sqrt(1+kd2);
			T d = (T) std::cos(b*(j-1)*(2*i-1));
			sum += x[j-1]*c*d;
		}
		dctx[i-1] = a*sum;
	}

}

// DCT-1D Type-IV
template<class T, int N>
void DCT_1D_TYPE4(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*4);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		for(j=1;j<=N;j++){
			T c = (T) std::cos(b*(2*j-1)*(2*i-1));
			sum += x[j-1]*c;
		}
		dctx[i-1] = a*sum;
	}

}

// DCT-2D
template<class T, int M, int N>
void DCT_2D(T x[M][N], T dctx[M][N]){

	int p, q, m, n;
	T a = (T) std::sqrt(1/M);
	T b = (T) std::sqrt(2/M);
	T c = (T) std::sqrt(1/N);
	T d = (T) std::sqrt(2/N);
	for(p=0;p<M;p++){
		T ap = (p==0)?a:b;
		for(q=0;q<N;q++){
			T sum = (T) 0;
			T aq = (q==0)?c:d;
			for(m=0;m<M;m++){
				T e = std::cos((2*m+1)*pi*p/(2*M));
				for(n=0;n<N;n++){
					T f = std::cos((2*n+1)*pi*q/(2*N));
					sum += x[m][n]*e*f;
				}
			}
			dctx[p][q] = ap*aq*sum;
		}
	}
}

// Inverse DCT-1D Type-I
template<class T, int N>
void iDCT_1D_TYPE1(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/(N-1));
	T b = (T) pi/(N-1);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		T kd1 = (i==1||i==N)?1:0;
		T c = (T) 1/std::sqrt(1+kd1);
		for(j=1;j<=N;j++){
			T kd2 = (j==1||j==N)?1:0;
			T d = (T) 1/std::sqrt(1+kd2);
			T e = (T) std::cos(b*(j-1)*(i-1));
			sum += x[j-1]*c*d*e;
		}
		dctx[i-1] = a*sum;
	}

}

// Inverse DCT-1D Type-II
template<class T, int N>
void iDCT_1D_TYPE2(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*2);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		for(j=1;j<=N;j++){
			T kd2 = (j==1)?1:0;
			T c = (T) 1/std::sqrt(1+kd2);
			T d = (T) std::cos(b*(j-1)*(2*i-1));
			sum += x[j-1]*c*d;
		}
		dctx[i-1] = a*sum;
	}

}

// Inverse DCT-1D Type-III
template<class T, int N>
void iDCT_1D_TYPE3(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*2);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		T kd1 = (i==1)?1:0;
		T c = (T) 1/std::sqrt(1+kd1);
		for(j=1;j<=N;j++){
			T d = (T) std::cos(b*(2*j-1)*(i-1));
			sum += x[j-1]*c*d;
		}
		dctx[i-1] = a*sum;
	}

}

// Inverse DCT-1D Type-IV
template<class T, int N>
void iDCT_1D_TYPE4(T x[N], T dctx[N]){

	int i, j;
	T a = (T) std::sqrt(2/N);
	T b = (T) pi/(N*4);
	for(i=1;i<=N;i++){
		T sum = (T) 0;
		for(j=1;j<=N;j++){
			T c = (T) std::cos(b*(2*j-1)*(2*i-1));
			sum += x[j-1]*c;
		}
		dctx[i-1] = a*sum;
	}

}


// Inverse DCT-2D
template<class T, int M, int N>
void iDCT_2D(T x[M][N], T dctx[M][N]){

	int p, q, m, n;
	T a = (T) std::sqrt(1/M);
	T b = (T) std::sqrt(2/M);
	T c = (T) std::sqrt(1/N);
	T d = (T) std::sqrt(2/N);
	for(m=0;m<M;m++){
		for(n=0;n<N;n++){
			T sum = (T) 0;
			for(p=0;p<M;p++){
				T ap = (p==0)?a:b;
				T e = std::cos((2*m+1)*pi*p/(2*M));
				for(q=0;q<N;q++){
					T aq = (q==0)?c:d;
					T f = std::cos((2*n+1)*pi*q/(2*N));
					sum += ap*aq*x[p][q]*e*f;
				}
			}
			dctx[m][n] = sum;
		}
	}
}


#endif /* HEADER_DSP_DCT_HPP_ */
