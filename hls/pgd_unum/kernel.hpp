
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "ap_int.h"
#include "ap_fixed.h"
#include "softposit.hpp"

#define ROW 1024 // matrix row number
#define COL 1024 // matrix column number
#define DIAG (COL<ROW ? COL : ROW) // diagonal matrix size

#define SOFT_POSIT_P1RECISION 3// SoftPosit precision 0:posit32; 1:posit16; 2:posit8; 3:positx
#define TOTALBITS 12

#define BOX_CONST 10 // Box constraint on the variable


template<class T, int M, T spfuncsub(T, T)>
void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncsub(V1[i], V2[i]);
	}
}
template<class T, int M, T spfuncsub(T, T, int), int BITS>
void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncsub(V1[i], V2[i], BITS);
	}
}

template<class T, int M, T spfuncmul(T, T)>
void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncmul(V1[i], S);
	}
}
template<class T, int M, T spfuncmul(T, T, int), int BITS>
void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncmul(V1[i], S, BITS);
	}
}

template<class T, int M, bool spfunclt(T, T)>
void VEC_SCALAR_MAX( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = !spfunclt(V1[i], S) ? V1[i] : S;
	}
}

template<class T, int M, bool spfunclt(T, T)>
void VEC_SCALAR_MIN( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfunclt(V1[i], S) ? V1[i] : S;
	}
}

template<class T, int M, T spfuncadd(T, T)>
void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncadd(V1[i], V2[i]);
	}
}
template<class T, int M, T spfuncadd(T, T, int), int BITS>
void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncadd(V1[i], V2[i], BITS);
	}
}

template<class T, int M, int N, T spfunc(double),
		T spfuncmul(T, T), T spfuncadd(T, T)>
void MAT_VEC_MUL(T A[M][N],
									  T B[N],
									  T C[M]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
		T tmp = spfunc(0);
		C[i] = spfunc(0);
		for ( int j=0; j<N; j++ ){
			tmp = spfuncmul(A[i][j], B[j]);
			C[i] = spfuncadd(C[i], tmp);
		}
	}
}
template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
		T spfuncmul(T, T, int), T spfuncadd(T, T, int),
		int BITS>
void MAT_VEC_MUL(T A[M][N],
									  T B[N],
									  T C[M]){
#pragma HLS INLINE off
	T zero = spfunc1(spfunc(0, BITS), BITS);
	for ( int i=0; i<M; i++ ){
		T tmp = zero;
		C[i] = zero;
		for ( int j=0; j<N; j++ ){
			tmp = spfuncmul(A[i][j], B[j], BITS);
			C[i] = spfuncadd(C[i], tmp, BITS);
		}
	}
}

template<class T, int M>
void VEC_EQ( T V1[M], T V2[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V2[i] = V1[i];
	}
}

template<class T, int M, T spfuncsub(T, T)>
void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncsub(V1[i], S);
	}
}
template<class T, int M, T spfuncsub(T, T, int), int BITS>
void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncsub(V1[i], S, BITS);
	}
}

template<class T, int M, T spfuncadd(T, T)>
void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncadd(V1[i], S);
	}
}
template<class T, int M, T spfuncadd(T, T, int), int BITS>
void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncadd(V1[i], S, BITS);
	}
}

template<class T, int M, T spfunc(double), T spfuncmul(T, T)>
void VEC_MINUS( T V1[M], T V2[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V2[i] = spfuncmul(spfunc(-1),V1[i]);
	}
}
template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), T spfuncmul(T, T, int), int BITS>
void VEC_MINUS( T V1[M], T V2[M] ){
#pragma HLS INLINE off
	T minusone = spfunc1(spfunc(-1, BITS), BITS);
	for( int i=0; i<M; i++ ){
		V2[i] = spfuncmul(minusone,V1[i], BITS);
	}
}

template<class T, unsigned int M, unsigned int N>
void MAT_TRANS(T Mat[M][N], T MatT[N][M]){
#pragma HLS INLINE off
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<N;j++){
			MatT[i][j] = Mat[j][i];
		}
	}
}

template<class T, int M, int N, int P, T spfunc(double),
		T spfuncmul(T, T), T spfuncadd(T, T)>
void MAT_MUL(T A[M][N],
			T B[N][P],
			T C[M][P]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
		for ( int j=0; j<P; j++ ){
			T tmp = spfunc(0);
			C[i][j] = spfunc(0);
			for ( int k=0; k<N; k++ ){
				tmp = spfuncmul(A[i][k], B[k][j]);
				C[i][j] = spfuncadd(C[i][j], tmp);
			}
		}
	}
}
template<class T, class T1, int M, int N, int P, T1 spfunc(double, int), T spfunc1(T1, int),
		T spfuncmul(T, T, int), T spfuncadd(T, T, int), int BITS>
void MAT_MUL(T A[M][N],
											    T B[N][P],
												T C[M][P]){
#pragma HLS INLINE off
	T zero = spfunc1(spfunc(0, BITS), BITS);
	for ( int i=0; i<M; i++ ){
		for ( int j=0; j<P; j++ ){
			T tmp = zero;
			C[i][j] = zero;
			for ( int k=0; k<N; k++ ){
				tmp = spfuncmul(A[i][k], B[k][j], BITS);
				C[i][j] = spfuncadd(C[i][j], tmp, BITS);
			}
		}
	}
}

template<class T, int M, int N, T spfunc(double)>
void IDENDTITY_MAT( T A[M][N] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		for( int j=0; j<N; j++ ){
			if(i==j)
				A[i][j] = spfunc(1);
			else
				A[i][j] = spfunc(0);
		}
	}
}
template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
void IDENDTITY_MAT( T A[M][N] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		for( int j=0; j<N; j++ ){
			if(i==j)
				A[i][j] = spfunc1(spfunc(1, BITS), BITS);
			else
				A[i][j] = spfunc1(spfunc(0, BITS), BITS);
		}
	}
}

template<class T, int M, int N,
	 	 	 T spfuncmul(T, T)>
void MAT_SCALAR_DOTMUL(T MatA[M][N],
					   T B,
					   T MatC[M][N]){
#pragma HLS INLINE off
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			MatC[i][j] = spfuncmul(MatA[i][j], B);
}
template<class T, int M, int N,
	 	 	 T spfuncmul(T, T, int), int BITS>
void MAT_SCALAR_DOTMUL(T MatA[M][N],
					   T B,
					   T MatC[M][N]){
#pragma HLS INLINE off
	for(int i=0;i<M;i++)
		for(int j=0;j<N;j++)
			MatC[i][j] = spfuncmul(MatA[i][j], B, BITS);
}

template<class T, int M, int N, T spfuncadd(T, T)>
void MAT_ADD(T A[M][N],
								      T B[M][N],
								      T C[M][N]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
		for ( int j=0; j<N; j++ ){
			C[i][j] = spfuncadd(A[i][j], B[i][j]);
		}
	}
}
template<class T, int M, int N, T spfuncadd(T, T, int), int BITS>
void MAT_ADD(T A[M][N],
								      T B[M][N],
								      T C[M][N]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
		for ( int j=0; j<N; j++ ){
			C[i][j] = spfuncadd(A[i][j], B[i][j], BITS);
		}
	}
}

template<class T, int M, T spfunc(double)>
void ZEROS_VEC( T V[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V[i] = spfunc(0);
	}
}
template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
void ZEROS_VEC( T V[M] ){
#pragma HLS INLINE off
	T zero = spfunc1(spfunc(0, BITS), BITS);
	for( int i=0; i<M; i++ ){
		V[i] = zero;
	}
}

template<class T, int M, T spfunc(double),
		T spfuncmul(T, T), T spfuncadd(T, T),
		T spfuncsqrt(T)>
void VEC_NORM( T V1[M], T &S ){
#pragma HLS INLINE off
	T tmp = spfunc(0);
	S = spfunc(0);
	for( int i=0; i<M; i++ ){
		tmp = spfuncmul(V1[i], V1[i]);
		S = spfuncadd(S, tmp);
	}
	S = spfuncsqrt(S);
}
template<class T, class T1, int M, T1 spfunc(double, int),
		T spfunc1(T1, int),T spfuncmul(T, T, int), T spfuncadd(T, T, int),
		T spfuncsqrt(T, int), int BITS>
void VEC_NORM( T V1[M], T &S ){
#pragma HLS INLINE off
	T tmp = spfunc1(spfunc(0, BITS), BITS);
	S = spfunc1(spfunc(0, BITS), BITS);
	for( int i=0; i<M; i++ ){
		tmp = spfuncmul(V1[i], V1[i], BITS);
		S = spfuncadd(S, tmp, BITS);
	}
	S = spfuncsqrt(S, BITS);
}

template<class T, int M, T spfuncdiv(T,T)>
void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncdiv(V1[i], V2[i]);
	}
}
template<class T, int M, T spfuncdiv(T,T,int), int BITS>
void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = spfuncdiv(V1[i], V2[i], BITS);
	}
}

template<class T, unsigned int M, T spfunc(double),
		 T spfuncadd(T, T), T spfuncsub(T, T),
		 T spfuncmul(T, T), T spfuncdiv(T, T),
		 T spfuncsqrt(T), bool spfunclt(T, T)>
void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){
#pragma HLS INLINE off

	// copy the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = spfunc(0);
		}
	}

	// decomposition in-place
	for(unsigned int j=0;j<M;j++){
		// compute the diagonal element
		T LL = Mat[j][j];
		for(unsigned int k=0;k<j;k++){
			T tmp = spfuncmul(MatL[j][k],MatL[j][k]);
			LL = spfuncsub(LL, tmp);
//			if(spfunclt(LL,spfunc(0))){
//				exit(-1);
//			}
		}
		MatL[j][j] = spfuncsqrt(LL);

		// compute the non-diagonal element
		T inv = spfuncdiv(spfunc(1), MatL[j][j]);
		for(unsigned int i=j+1;i<M;i++){
			LL = Mat[i][j];
//				std::cout << LL;
			for(unsigned int k=0;k<j;k++){
				T tmp = spfuncmul(MatL[i][k], MatL[j][k]);
				LL = spfuncsub(LL, tmp);
			}
			MatL[i][j] = spfuncmul(LL, inv);
		}
	}

	// transpose L to get U
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatU[i][j] = MatL[j][i];
		}
	}
}
template<class T, class T1, unsigned int M, T1 spfunc(double, int), T spfunc1(T1, int),
		 T spfuncadd(T, T, int), T spfuncsub(T, T, int),
		 T spfuncmul(T, T, int), T spfuncdiv(T, T, int),
		 T spfuncsqrt(T, int), bool spfunclt(T, T), int BITS>
void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){
#pragma HLS INLINE off

	T zero = spfunc1(spfunc(0, BITS), BITS);
	T one = spfunc1(spfunc(1, BITS), BITS);

	// copy the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = zero;
		}
	}

	// decomposition in-place
	for(unsigned int j=0;j<M;j++){
		// compute the diagonal element
		T LL = Mat[j][j];
		for(unsigned int k=0;k<j;k++){
			T tmp = spfuncmul(MatL[j][k],MatL[j][k], BITS);
			LL = spfuncsub(LL, tmp, BITS);
//			if(spfunclt(LL,zero)){
//				exit(-1);
//			}
		}
		MatL[j][j] = spfuncsqrt(LL, BITS);

		// compute the non-diagonal element
		T inv = spfuncdiv(one, MatL[j][j], BITS);
		for(unsigned int i=j+1;i<M;i++){
			LL = Mat[i][j];
//				std::cout << LL;
			for(unsigned int k=0;k<j;k++){
				T tmp = spfuncmul(MatL[i][k], MatL[j][k], BITS);
				LL = spfuncsub(LL, tmp, BITS);
			}
			MatL[i][j] = spfuncmul(LL, inv, BITS);
		}
	}

	// transpose L to get U
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatU[i][j] = MatL[j][i];
		}
	}
}

template<class T, int M, int N, T spfunc(double),
		T spfuncmul(T,T), T spfuncadd(T,T),
		T spfuncsub(T,T), T spfuncsqrt(T),
		T spfuncdiv(T,T), bool spfunceq(T, T)>
void QRD_HH(T Mat[M][N],
		    T MatQ[M][N],
		    T MatR[N][N]){
#pragma HLS INLINE off
	int i,j,k;
	//R=A;
	for(j=0;j<M;j++){
		for(i=0;i<N;i++){
			MatR[j][i] = Mat[j][i];
			//Q=eye(m);
			if(i==j)
				MatQ[j][i] = spfunc(1);
			else
				MatQ[j][i] = spfunc(0);
		}
	}

	T g, s;
	T x[M], v[M], w[M], u[N];
	T tmp1[M][N], tmp2[M][N];
	for(k=0;k<M-1;k++){
		// x=zeros(m,1);
		for(j=0;j<M;j++){
			x[j] = spfunc(0);
		}
		ZEROS_VEC<T, M, spfunc>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_NORM<T,M,spfunc,spfuncmul,spfuncadd,spfuncsqrt>(x, g);
		// v=x;
		VEC_EQ<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = spfuncadd(x[k], g);
		//s=norm(v);
		VEC_NORM<T, M,spfunc,spfuncmul,spfuncadd,spfuncsqrt>(v, s);
		if(!spfunceq(s, spfunc(0))){
			// w=v/s;
			VEC_DIV<T, M, spfuncdiv>(v, s, w);
			// u=2*R'*w;
			for(i=0;i<N;i++){
				u[i] = spfunc(0);
				for(j=0;j<M;j++){
					u[i] = spfuncadd(u[i], spfuncmul(spfuncmul(spfunc(2), MatR[j][i]), w[j]));
				}
			}
			// R=R-w*u'; %Product HR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatR[j][i] = spfuncsub(MatR[j][i], spfuncmul(w[j], u[i]));
				}
			}
			// Q=Q-2*Q*w*w'; %Product QR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					tmp1[j][i] = spfuncmul(w[j], w[i]);
				}
			}
			MAT_MUL<T, M, N, N, spfunc, spfuncmul, spfuncadd>(MatQ, tmp1, tmp2);
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatQ[j][i] = spfuncsub(MatQ[j][i], spfuncmul(spfunc(2), tmp2[j][i]));
				}
			}
		}
	}

}
template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
		T spfuncmul(T,T, int), T spfuncadd(T,T, int),
		T spfuncsub(T,T, int), T spfuncsqrt(T, int),
		T spfuncdiv(T,T, int), bool spfunceq(T, T), int BITS>
void QRD_HH(T Mat[M][N],
		    T MatQ[M][N],
		    T MatR[N][N]){
#pragma HLS INLINE off
	int i,j,k;
	//R=A;
	T zero = spfunc1(spfunc(0, BITS), BITS);
	T one = spfunc1(spfunc(1, BITS), BITS);
	T two = spfunc1(spfunc(2, BITS), BITS);
	for(j=0;j<M;j++){
		for(i=0;i<N;i++){
			MatR[j][i] = Mat[j][i];
			//Q=eye(m);
			if(i==j)
				MatQ[j][i] = one;
			else
				MatQ[j][i] = zero;
		}
	}

	T g, s;
	T x[M], v[M], w[M], u[N];
	T tmp1[M][N], tmp2[M][N];
	for(k=0;k<M-1;k++){
		// x=zeros(m,1);
		for(j=0;j<M;j++){
			x[j] = zero;
		}
		ZEROS_VEC<T, T1, M, spfunc, spfunc1, BITS>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_NORM<T,T1,M,spfunc,spfunc1,spfuncmul,spfuncadd,spfuncsqrt,BITS>(x, g);
		// v=x;
		VEC_EQ<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = spfuncadd(x[k], g, BITS);
		//s=norm(v);
		VEC_NORM<T,T1,M,spfunc,spfunc1,spfuncmul,spfuncadd,spfuncsqrt,BITS>(v, s);
		if(!spfunceq(s, zero)){
			// w=v/s;
			VEC_DIV<T, M, spfuncdiv, BITS>(v, s, w);
			// u=2*R'*w;
			for(i=0;i<N;i++){
				u[i] = zero;
				for(j=0;j<M;j++){
					u[i] = spfuncadd(u[i], spfuncmul(spfuncmul(two, MatR[j][i], BITS), w[j], BITS), BITS);
				}
			}
			// R=R-w*u'; %Product HR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatR[j][i] = spfuncsub(MatR[j][i], spfuncmul(w[j], u[i], BITS), BITS);
				}
			}
			// Q=Q-2*Q*w*w'; %Product QR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					tmp1[j][i] = spfuncmul(w[j], w[i], BITS);
				}
			}
			MAT_MUL<T, T1, M, N, N, spfunc, spfunc1, spfuncmul, spfuncadd, BITS>(MatQ, tmp1, tmp2);
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatQ[j][i] = spfuncsub(MatQ[j][i], spfuncmul(two, tmp2[j][i], BITS), BITS);
				}
			}
		}
	}

}

template<class T, int M, T spfunc(double),
		T spfuncmul(T, T), T spfuncadd(T, T),
		T spfuncsub(T, T), T spfuncdiv(T, T)>
void UPTRIANGULARMATINV(T R[M][M],T Ri[M][M]){
#pragma HLS INLINE off
	int i=0,j=0,k=0;
	// R inversion
	for(i=0;i<M;i++){
		for(j=0;j<M;j++){
			Ri[i][j]=spfunc(0);
		}
	}
	for(i=0;i<M;i++){
		Ri[i][i]=spfuncdiv(spfunc(1),R[i][i]);
		for(j=i+1;j<M;j++){
			for(k=0;k<=j-1;k++){
				Ri[i][j] = spfuncsub(Ri[i][j], spfuncdiv(spfuncmul(Ri[i][k], R[k][j]), R[j][j]));
			}
		}
	}
}
template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
		T spfuncmul(T, T, int), T spfuncadd(T, T, int),
		T spfuncsub(T, T, int), T spfuncdiv(T, T, int), int BITS>
void UPTRIANGULARMATINV(T R[M][M],T Ri[M][M]){
#pragma HLS INLINE off

	T zero = spfunc1(spfunc(0, BITS), BITS);
	T one = spfunc1(spfunc(1, BITS), BITS);

	int i=0,j=0,k=0;
	// R inversion
	for(i=0;i<M;i++){
		for(j=0;j<M;j++){
			Ri[i][j]=zero;
		}
	}
	for(i=0;i<M;i++){
		Ri[i][i]=spfuncdiv(one,R[i][i], BITS);
		for(j=i+1;j<M;j++){
			for(k=0;k<=j-1;k++){
				Ri[i][j] = spfuncsub(Ri[i][j], spfuncdiv(spfuncmul(Ri[i][k], R[k][j], BITS), R[j][j], BITS), BITS);
			}
		}
	}
}
template<class T, int M>
void ORTHOGONALMATINV(T Q[M][M],T Qi[M][M]){
#pragma HLS INLINE off
	int i=0,j=0;
	// Q inversion
	for(i=0;i<M;i++){
		for(j=0;j<M;j++){
			Qi[i][j] = Q[j][i];
		}
	}
}
template<class T, int M, T spfunc(double),
		T spfuncmul(T,T), T spfuncadd(T,T),
		T spfuncsub(T,T), T spfuncsqrt(T),
		T spfuncdiv(T,T), bool spfunceq(T,T)>
void MAT_QRINV(T A[M][M], T B[M][M]){
#pragma HLS INLINE off
	T Q[M][M], R[M][M];
	T Qi[M][M], Ri[M][M];
	QRD_HH<T, M, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq>(A, Q, R);
	UPTRIANGULARMATINV<T, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncdiv>(R, Ri);
	ORTHOGONALMATINV<T, M>(Q, Qi);
}
template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
		T spfuncmul(T,T, int), T spfuncadd(T,T, int),
		T spfuncsub(T,T, int), T spfuncsqrt(T, int),
		T spfuncdiv(T,T, int), bool spfunceq(T,T), int BITS>
void MAT_QRINV(T A[M][M], T B[M][M]){
#pragma HLS INLINE off
	T Q[M][M], R[M][M];
	T Qi[M][M], Ri[M][M];
	QRD_HH<T, T1, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq, BITS>(A, Q, R);
	UPTRIANGULARMATINV<T, T1, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncdiv, BITS>(R, Ri);
	ORTHOGONALMATINV<T, M>(Q, Qi);
	MAT_MUL<T, T1, M, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, BITS>(Ri, Qi, B);
}
