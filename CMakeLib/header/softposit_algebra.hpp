/*
 * softposit_algebra.hpp
 *
 *  Created on: Aug 22, 2019
 *      Author: yunwu
 */

#ifndef SRC_SOFTPOSIT_ALGEBRA_HPP_
#define SRC_SOFTPOSIT_ALGEBRA_HPP_

#include "data.hpp"
#include "common.hpp"
#include "softposit.hpp"

class SoftPosit_Algebra{
public:


	// Generate all zero matrix
	template<class T, int M, int N, T spfunc(double)>
	void ZEROS_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = spfunc(0);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double, int), int BITS>
	void ZEROS_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = spfunc(0, BITS);
			}
		}
	}

	// Generate all one matrix
	template<class T, int M, int N, T spfunc(double)>
	void ONES_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = spfunc(1);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double, int), int BITS>
	void ONES_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = spfunc(1, BITS);
			}
		}
	}

	// Generate identity matrix
	template<class T, int M, int N, T spfunc(double)>
	void IDENDTITY_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				if(i==j)
					A[i][j] = spfunc(1);
				else
					A[i][j] = spfunc(0);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double)>
	void IDENDTITY_MAT( T **A ){
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
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				if(i==j)
					A[i][j] = spfunc1(spfunc(1, BITS), BITS);
				else
					A[i][j] = spfunc1(spfunc(0, BITS), BITS);
			}
		}
	}
	template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
	void IDENDTITY_MAT( T **A ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				if(i==j)
					A[i][j] = spfunc1(spfunc(1, BITS), BITS);
				else
					A[i][j] = spfunc1(spfunc(0, BITS), BITS);
			}
		}
	}

	// Generate random matrix
	template<class T, int M, int N, T spfunc(double)>
	void RND_MAT( T A[M][N] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
				double rndnum = (double)( rnd ) / FLOAT_SIZE;
				A[i][j] = spfunc(rndnum);
			}
		}
	}

	// Generate random matrix
	template<class T, int M, int N, T spfunc(double)>
	void RND_DIAGMAT( T A[M][N] ){
		srand (time(NULL));
		int sparse_num = floor( DIAG_RATIO * M );
		for( int i=0; i<M; i++){
			for( int j=0; j<N; j++ ){
				if( i<sparse_num && i==j){
					int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
					double rndnum = (double)( rnd ) / FLOAT_SIZE;
					A[i][j] = spfunc(rndnum);
				}else
					A[i][j] = spfunc(0);
			}
		}
	}

	// Generate all zero vector
	template<class T, int M, T spfunc(double)>
	void ZEROS_VEC( T V[M] ){
		for( int i=0; i<M; i++ ){
			V[i] = spfunc(0);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
	void ZEROS_VEC( T V[M] ){
		T zero = spfunc1(spfunc(0, BITS), BITS);
		for( int i=0; i<M; i++ ){
			V[i] = zero;
		}
	}

	// Generate all one vector
	template<class T, int M, T spfunc(double)>
	void ONES_VEC( T V[M] ){
		for( int i=0; i<M; i++ ){
			V[i] = spfunc(1);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
	void ONES_VEC( T V[M] ){
		T one = spfunc1(spfunc(1, BITS), BITS);
		for( int i=0; i<M; i++ ){
			V[i] = one;
		}
	}

	// Generate random vector
	template<class T, int M, T spfunc(double)>
	void RND_VEC( T V[M] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
			double rndnum = (double)( rnd ) / FLOAT_SIZE;
			V[i] = spfunc(rndnum);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), int BITS>
	void RND_VEC( T V[M] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
			double rndnum = (double)( rnd ) / FLOAT_SIZE;
			V[i] = spfunc1(spfunc(rndnum, BITS), BITS);
		}
	}

	// Transfer vector values
	template<class T, int M>
	void VEC_EQ( T V1[M], T V2[M] ){
		for( int i=0; i<M; i++ ){
			V2[i] = V1[i];
		}
	}
	template<class T1, class T2, int M>
	void VEC_EQ( T1 V1[M], T2 V2[M] ){
		for( int i=0; i<M; i++ ){
			V2[i] = (T2) V1[i];
		}
	}

	// reverse the sign of vector elements
	template<class T, int M, T spfunc(double), T spfuncmul(T, T)>
	void VEC_MINUS( T V1[M], T V2[M] ){
		for( int i=0; i<M; i++ ){
			V2[i] = spfuncmul(spfunc(-1),V1[i]);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int), T spfuncmul(T, T, int), int BITS>
	void VEC_MINUS( T V1[M], T V2[M] ){
		T minusone = spfunc1(spfunc(-1, BITS), BITS);
		for( int i=0; i<M; i++ ){
			V2[i] = spfuncmul(minusone,V1[i], BITS);
		}
	}

	// Vector addition
	template<class T, int M, T spfuncadd(T, T)>
	void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncadd(V1[i], V2[i]);
		}
	}
	template<class T, int M, T spfuncadd(T, T, int), int BITS>
	void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncadd(V1[i], V2[i], BITS);
		}
	}
	template<class T1, class T2, class T3, int M, T3 spfuncadd(T3, T3)>
	void VEC_ADD( T1 V1[M], T2 V2[M], T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncadd(V1[i], V2[i]);
		}
	}

	// Vector subtraction
	template<class T, int M, T spfuncsub(T, T)>
	void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncsub(V1[i], V2[i]);
		}
	}
	template<class T, int M, T spfuncsub(T, T, int), int BITS>
	void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncsub(V1[i], V2[i], BITS);
		}
	}

	template<class T, unsigned int M, unsigned int N>
	void MAT_TRANS(T Mat[M][N], T MatT[N][M]){
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<N;j++){
				MatT[i][j] = Mat[j][i];
			}
		}
	}
	template<class T, unsigned int M, unsigned int N>
	void MAT_TRANS(T **Mat, T **MatT){
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<N;j++){
				MatT[i][j] = Mat[j][i];
			}
		}
	}

	// Vector multiplication
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void VEC_MUL_2SCALAR( T V1[M], T V2[M], T S ){
		T tmp = spfunc(0);
		S = spfunc(0);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V2[i]);
			S = spfuncadd(S, tmp);
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_MUL_2SCALAR( T1 V1[M], T2 V2[M], T3 S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V2[i];
		}
	}
	template<class T, int M, int N, T spfuncmul(T, T)>
	void VEC_MUL_2MATRIX( T V1[M], T V2[N], T M3[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				M3[i][j] = spfuncmul(V1[i], V2[j]);
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void VEC_MUL_2MATRIX( T1 V1[M], T2 V2[N], T3 M3[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				M3[i][j] = V1[i] * V2[j];
			}
		}
	}

	// Vector division
	template<class T, int M, T spfuncdiv(T,T)>
	void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncdiv(V1[i], V2[i]);
		}
	}
	template<class T, int M, T spfuncdiv(T,T,int), int BITS>
	void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncdiv(V1[i], V2[i], BITS);
		}
	}
	template<class T, int M, T spfuncdiv(T,T)>
	void VEC_DIV( T V1[M], T v, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncdiv(V1[i], v);
		}
	}
	template<class T, int M, T spfuncdiv(T,T,int), int BITS>
	void VEC_DIV( T V1[M], T v, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncdiv(V1[i], v, BITS);
		}
	}

	// Vector norm
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T),
			T spfuncsqrt(T)>
	void VEC_NORM( T V1[M], T &S ){
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
		T tmp = spfunc1(spfunc(0, BITS), BITS);
		S = spfunc1(spfunc(0, BITS), BITS);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V1[i], BITS);
			S = spfuncadd(S, tmp, BITS);
		}
		S = spfuncsqrt(S, BITS);
	}
	template<class T1, class T2, int M>
	void VEC_NORM( T1 V1[M], T2 &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
		double tmp = (double) S;
		S = (T2) std::sqrt(S);
	}

	// Vector square norm
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void VEC_NORM2( T V1[M], T &S ){
		T tmp = spfunc(0);
		S = spfunc(0);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V1[i]);
			S = spfuncadd(S, tmp);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int), int BITS>
	void VEC_NORM2( T V1[M], T &S ){
		T tmp = spfunc1(spfunc(0, BITS), BITS);
		S = spfunc1(spfunc(0, BITS), BITS);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V1[i], BITS);
			S = spfuncadd(S, tmp, BITS);
		}
	}
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void VEC_SQUARE_NORM( T V1[M], T &S ){
		T tmp = spfunc(0);
		S = spfunc(0);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V1[i]);
			S = spfuncadd(S, tmp);
		}
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int), int BITS>
	void VEC_SQUARE_NORM( T V1[M], T &S ){
		T tmp = spfunc1(spfunc(0, BITS), BITS);
		S = spfunc1(spfunc(0, BITS), BITS);
		for( int i=0; i<M; i++ ){
			tmp = spfuncmul(V1[i], V1[i], BITS);
			S = spfuncadd(S, tmp, BITS);
		}
	}
	template<class T1, class T2, int M>
	void VEC_SQUARE_NORM( T1 V1[M], T2 &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
	}

	// Vector scalar addition
	template<class T, int M, T spfuncadd(T, T)>
	void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
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
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_ADD( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] + S;
		}
	}

	// Vector scalar subtraction
	template<class T, int M, T spfuncsub(T, T)>
	void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncsub(V1[i], S);
		}
	}
	template<class T, int M, T spfuncsub(T, T, int), int BITS>
	void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncsub(V1[i], S, BITS);
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_SUB( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] - S;
		}
	}

	// Vector scalar multiplication
	template<class T, int M, T spfuncmul(T, T)>
	void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncmul(V1[i], S);
		}
	}
	template<class T, int M, T spfuncmul(T, T, int), int BITS>
	void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncmul(V1[i], S, BITS);
		}
	}
	template<class T, int M, T spfuncmul(T, T)>
	void VEC_SCALAR_MUL( T S, T V1[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncmul(V1[i], S);
		}
	}
	template<class T, int M, T spfuncmul(T, T, int), int BITS>
	void VEC_SCALAR_MUL( T S, T V1[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfuncmul(V1[i], S, BITS);
		}
	}

	// Vector scalar compare
	template<class T, int M, bool spfunclt(T, T)>
	void VEC_SCALAR_MIN( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfunclt(V1[i], S) ? V1[i] : S;
		}
	}
	template<class T, int M, bool spfunclt(T, T)>
	void VEC_SCALAR_MIN( T S, T V1[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = spfunclt(V1[i], S) ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MIN( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] < S ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MIN( T2 S, T1 V1[M], T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] < S ? V1[i] : S;
		}
	}
	template<class T, int M, bool spfunclt(T, T)>
	void VEC_SCALAR_MAX( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = !spfunclt(V1[i], S) ? V1[i] : S;
		}
	}
	template<class T, int M, bool spfunclt(T, T)>
	void VEC_SCALAR_MAX( T S, T V1[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = !spfunclt(V1[i], S) ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MAX( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] > S ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MAX( T2 S, T1 V1[M], T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] > S ? V1[i] : S;
		}
	}

	// The basic matrix addition with complexity of O(NN)
	template<class T, int M, int N, T spfuncadd(T, T)>
	void MAT_ADD(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = spfuncadd(A[i][j], B[i][j]);
			}
		}
	}
	template<class T, int M, int N, T spfuncadd(T, T)>
	void MAT_ADD(T **A,
									      T **B,
									      T **C){
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
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = spfuncadd(A[i][j], B[i][j], BITS);
			}
		}
	}
	template<class T, int M, int N, T spfuncadd(T, T, int), int BITS>
	void MAT_ADD(T **A,
									      T **B,
									      T **C){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = spfuncadd(A[i][j], B[i][j], BITS);
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void MAT_ADD(T1 A[M][N],
									      T2 B[M][N],
									      T3 C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] + B[i][j];
			}
		}
	}

	// The basic matrix subtraction with complexity of O(NN)
	template<class T, int M, int N, T spfuncsub(T, T)>
	void MAT_SUB(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = spfuncsub(A[i][j], B[i][j]);
			}
		}
	}
	template<class T, int M, int N, T spfuncsub(T, T, int), int BITS>
	void MAT_SUB(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = spfuncsub(A[i][j], B[i][j], BITS);
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void MAT_SUB(T1 A[M][N],
									      T2 B[M][N],
									      T3 C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] - B[i][j];
			}
		}
	}

	// Matrix Vector multiplication
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void MAT_VEC_MUL(T A[M][N],
										  T B[N],
										  T C[M]){
		for ( int i=0; i<M; i++ ){
			T tmp = spfunc(0);
			C[i] = spfunc(0);
			for ( int j=0; j<N; j++ ){
				tmp = spfuncmul(A[i][j], B[j]);
				C[i] = spfuncadd(C[i], tmp);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void MAT_VEC_MUL(T **A,
										  T *B,
										  T *C){
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
	template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int),
			int BITS>
	void MAT_VEC_MUL(T **A,
										  T *B,
										  T *C){
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
	template<class T1, class T2, class T3, int M, int N>
	void MAT_VEC_MUL(T1 A[M][N],
										  T2 B[N],
										  T3 C[M]){
		for ( int i=0; i<M; i++ ){
			C[i] = 0;
			for ( int j=0; j<N; j++ ){
				C[i] += A[i][j] * B[j];
			}
		}
	}
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void VEC_MAT_MUL(T A[M],
										  T B[M][N],
										  T C[N]){
		for ( int i=0; i<N; i++ ){
			T tmp = spfunc(0);
			C[i] = spfunc(0);
			for ( int j=0; j<M; j++ ){
				tmp = spfuncmul(A[j], B[j][i]);
				C[i] = spfuncadd(C[i], tmp);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int),
			int BITS>
	void VEC_MAT_MUL(T A[M],
										  T B[M][N],
										  T C[N]){
		T zero = spfunc1(spfunc(0, BITS), BITS);
		for ( int i=0; i<N; i++ ){
			T tmp = zero;
			C[i] = zero;
			for ( int j=0; j<M; j++ ){
				tmp = spfuncmul(A[j], B[j][i], BITS);
				C[i] = spfuncadd(C[i], tmp, BITS);
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void VEC_MAT_MUL(T1 A[M],
										  T2 B[M][N],
										  T3 C[N]){
		for ( int i=0; i<N; i++ ){
			C[i] = 0;
			for ( int j=0; j<M; j++ ){
				C[i] += A[j] * B[j][i];
			}
		}
	}

	// The basic matrix value transfer
	template<class T, int M, int N>
	void MAT_EQ(T A[M][N],
									      T B[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				B[i][j] = A[i][j];
			}
		}
	}
	template<class T1, class T2, int M, int N>
	void MAT_EQ(T1 A[M][N],
									      T2 B[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				B[i][j] = (T2) A[i][j];
			}
		}
	}

	// The basic arbitrary matrix multiplication with complexity of O(MNP)
	template<class T, int M, int N, int P, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void MAT_MUL(T A[M][N],
												    T B[N][P],
													T C[M][P]){
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
	template<class T, int M, int N, int P, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void MAT_MUL(T **A,
												    T **B,
													T **C){
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
	template<class T, class T1, int M, int N, int P, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int), int BITS>
	void MAT_MUL(T **A,
												    T **B,
													T **C){
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
	template<class T1, class T2, class T3, int M, int N, int P>
	void MAT_MUL(T1 A[M][N],
												    T2 B[N][P],
													T3 C[M][P]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<P; j++ ){
				C[i][j] = 0;
				for ( int k=0; k<N; k++ ){
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
	}

	// The basic square matrix multiplication with complexity of O(M^3)
	// (must be 2^X square matrix)
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T)>
	void SQUARE_MAT_MUL(T A[M][M],
												    T B[M][M],
													T C[M][M]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<M; j++ ){
				T tmp = spfunc(0);
				C[i][j] = spfunc(0);
				for ( int k=0; k<M; k++ ){
					tmp = spfuncmul(A[i][k], B[k][j]);
					C[i][j] = spfuncadd(C[i][j], tmp);
				}
			}
		}
	}
	template<class T, int M, T spfunc(double, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int),
			int BITS>
	void SQUARE_MAT_MUL(T A[M][M],
												    T B[M][M],
													T C[M][M]){
		T zero = spfunc1(spfunc(0, BITS), BITS);
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<M; j++ ){
				T tmp = zero;
				C[i][j] = zero;
				for ( int k=0; k<M; k++ ){
					tmp = spfuncmul(A[i][k], B[k][j], BITS);
					C[i][j] = spfuncadd(C[i][j], tmp, BITS);
				}
			}
		}
	}
	template<class T1, class T2, class T3, int M>
	void SQUARE_MAT_MUL(T1 A[M][M],
												    T2 B[M][M],
													T3 C[M][M]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<M; j++ ){
				C[i][j] = 0;
				for ( int k=0; k<M; k++ ){
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
	}

	// Find the nearest 2^X size
	int FIND_NEAREST_2SQUARE(int m, int n, int p){
		int max1 = m > n ? m : n;
		int max2 = max1 > p ? max1 : p;
		int result = 2;
		while ( result<max2 ){
			result = result * 2;
		}
		return result;
	}
	int FIND_NEAREST_2SQUARE(int m, int n){
		int max1 = m > n ? m : n;
		int result = 2;
		while ( result<max1 ){
			result = result * 2;
		}
		return result;
	}
	int FIND_NEAREST_2SQUARE(int m){
		int result = 2;
		while ( result<m ){
			result = result * 2;
		}
		return result;
	}

	// Extend general matrix to square, the square size must be largest
	template<class T, int M, int N, int S>
	void EXTEND_TO_SQUARE(T A[M][N],
													   T B[S][S]){
		for ( int i=0; i< S; i++){
			for ( int j=0; j<S; j++){
				if( i<M && j<N )
					B[i][j] = A[i][j];
				else
					B[i][j] = 0;
			}
		}
	}
	template<class T1, class T2, int M, int N, int S>
	void EXTEND_TO_SQUARE(T1 A[M][N],
													   T2 B[S][S]){
		for ( int i=0; i< S; i++){
			for ( int j=0; j<S; j++){
				if( i<M && j<N )
					B[i][j] = A[i][j];
				else
					B[i][j] = 0;
			}
		}
	}

	// The Strassen matrix multiplication with complexity of O(n^2.8074)
	template<class T, int M>
	void SQUARE_MAT_MUL_STRASSEN(int m,
											   	      T A[M][M],
													  T B[M][M],
													  T C[M][M]){

	    T A11[M][M], A12[M][M], A21[M][M], A22[M][M];
	    T B11[M][M], B12[M][M], B21[M][M], B22[M][M];
	    T C11[M][M], C12[M][M], C21[M][M], C22[M][M];
	    T M1[M][M], M2[M][M], M3[M][M], M4[M][M], M5[M][M], M6[M][M], M7[M][M];
	    T AA[M][M], BB[M][M];

	    if(m == 2) {  //2-order
	    	SQUARE_MAT_MUL<T, 2>( A, B, C );
	    } else {
			// Divide and conquer
			for(int i=0; i<m/2; i++) {
				for(int j=0; j<m/2; j++) {
					A11[i][j] = A[i][j];
					A12[i][j] = A[i][j+m/2];
					A21[i][j] = A[i+m/2][j];
					A22[i][j] = A[i+m/2][j+m/2];

					B11[i][j] = B[i][j];
					B12[i][j] = B[i][j+m/2];
					B21[i][j] = B[i+m/2][j];
					B22[i][j] = B[i+m/2][j+m/2];
				}
			}

			//Calculate M1 = (A0 + A3) × (B0 + B3)
			MAT_ADD<T, m/2>( A11, A22, AA );
			MAT_ADD<T, m/2>( B11, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M1 );

			//Calculate M2 = (A2 + A3) × B0
			MAT_ADD<T, m/2>( A21, A22, AA );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, B11, M2 );

			//Calculate M3 = A0 × (B1 - B3)
			MAT_SUB<T, m/2>( B12, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, A11, BB, M3 );

			//Calculate M4 = A3 × (B2 - B0)
			MAT_SUB<T, m/2>( B21, B11, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, A22, BB, M4 );

			//Calculate M5 = (A0 + A1) × B3
			MAT_ADD<T, m/2>( A11, A12, AA );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, B22, M5);

			//Calculate M6 = (A2 - A0) × (B0 + B1)
			MAT_SUB<T, m/2>( A21, A11, AA );
			MAT_ADD<T, m/2>( B11, B12, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M6 );

			//Calculate M7 = (A1 - A3) × (B2 + B3)
			MAT_SUB<T, m/2>( A12, A22, AA );
			MAT_ADD<T, m/2>( B21, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M7 );

			//Calculate C0 = M1 + M4 - M5 + M7
			MAT_ADD<T, m/2>( M1, M4, AA );
			MAT_ADD<T, m/2>( M5, M7, BB );
			MAT_SUB<T, m/2>( AA, BB, C11 );

			//Calculate C1 = M3 + M5
			MAT_ADD<T, m/2>( M3, M5, C12 );

			//Calculate C2 = M2 + M4
			MAT_ADD<T, m/2>( M2, M4, C21 );

			//Calculate C3 = M1 - M2 + M3 + M6
			MAT_SUB<T, m/2>( M1, M2, AA);
			MAT_ADD<T, m/2>( M3, M6, BB);
			MAT_ADD<T, m/2>( AA, BB, C22);

			//Set the result to C[M][M]
			for(int i=0; i<m/2; i++) {
				for(int j=0; j<m/2; j++) {
					C[i][j] = C11[i][j];
					C[i][j+m/2] = C12[i][j];
					C[i+m/2][j] = C21[i][j];
					C[i+m/2][j+m/2] = C22[i][j];
				}
			}
		}
	}

	// QR Decomposition
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T,T), T spfuncadd(T,T),
			T spfuncsub(T,T), T spfuncsqrt(T),
			T spfuncdiv(T,T)>
	void QRD_GS(T Mat[M][N],
			T MatQ[M][N],
			T MatR[N][N]){
		// Initialisation
		int i = 0,j = 0,k = 0;
		for(i=0;i<N;i++){
			for(j=0;j<M;j++){
				MatQ[j][i] = Mat[j][i];
			}
			for(j=0;j<N;j++){
				MatR[i][j] = spfunc(0);
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				T mul = spfuncmul(Mat[j][i], Mat[j][i]);
				norm[i] = spfuncadd(norm[i], mul);
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = spfuncsqrt(norm[i]);
			for(k=0;k<M;k++){
				MatQ[k][i]=spfuncdiv(MatQ[k][i], MatR[i][i]);
			}
			for(j=i+1;j<M;j++){
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = spfuncadd(MatR[i][j], spfuncmul(MatQ[k][i], MatQ[k][j]));
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = spfuncsub(MatQ[k][j], spfuncmul(MatQ[k][i], MatR[i][j]));
				}
				// update norm
				norm[j] = spfuncsub(norm[j], spfuncmul(MatR[i][j], MatR[i][j]));
			}
		}
	}
	template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T,T, int), T spfuncadd(T,T, int),
			T spfuncsub(T,T, int), T spfuncsqrt(T, int),
			T spfuncdiv(T,T, int), int BITS>
	void QRD_GS(T Mat[M][N],
			T MatQ[M][N],
			T MatR[N][N]){
		// Initialisation
		int i = 0,j = 0,k = 0;
		T zero = spfunc1(spfunc(0, BITS), BITS);
		for(i=0;i<N;i++){
			for(j=0;j<M;j++){
				MatQ[j][i] = Mat[j][i];
			}
			for(j=0;j<N;j++){
				MatR[i][j] = zero;
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				T mul = spfuncmul(Mat[j][i], Mat[j][i], BITS);
				norm[i] = spfuncadd(norm[i], mul, BITS);
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = spfuncsqrt(norm[i], BITS);
			for(k=0;k<M;k++){
				MatQ[k][i]=spfuncdiv(MatQ[k][i], MatR[i][i], BITS);
			}
			for(j=i+1;j<M;j++){
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = spfuncadd(MatR[i][j], spfuncmul(MatQ[k][i], MatQ[k][j], BITS), BITS);
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = spfuncsub(MatQ[k][j], spfuncmul(MatQ[k][i], MatR[i][j], BITS), BITS);
				}
				// update norm
				norm[j] = spfuncsub(norm[j], spfuncmul(MatR[i][j], MatR[i][j], BITS), BITS);
			}
		}
	}
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T,T), T spfuncadd(T,T),
			T spfuncsub(T,T), T spfuncdiv(T,T),
			T spfuncsqrt(T)>
	void QRD_MGS(T Mat[M][N],
				 T MatQ[M][N],
				 T MatR[N][N]){
		// Initialisation
		int i = 0,j = 0,k = 0;
		for(i=0;i<N;i++){
			for(j=0;j<M;j++){
				MatQ[j][i] = Mat[j][i];
			}
			for(j=0;j<N;j++){
				MatR[i][j] = spfunc(0);
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				T mul = spfuncmul(Mat[j][i], Mat[j][i]);
				norm[i] = spfuncadd(norm[i], mul);
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = spfuncsqrt(norm[i]);
			for(k=0;k<M;k++){
				MatQ[k][i] = spfuncdiv(MatQ[k][i],MatR[i][i]);
			}
			for(j=i+1;j<M;j++){
				float tmp;
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = spfuncadd(MatR[i][j], spfuncmul(MatQ[k][i], MatQ[k][j]));
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = spfuncsub(MatQ[k][j], spfuncmul(MatQ[k][i], MatR[i][j]));
				}
			// update norm: no update for QR_MGS
				// norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
			}
		}

	}
	template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T,T, int), T spfuncadd(T,T, int),
			T spfuncsub(T,T, int), T spfuncdiv(T,T, int),
			T spfuncsqrt(T, int), int BITS>
	void QRD_MGS(T Mat[M][N],
				 T MatQ[M][N],
				 T MatR[N][N]){
		// Initialisation
		int i = 0,j = 0,k = 0;
		T zero = spfunc1(spfunc(0, BITS), BITS);
		for(i=0;i<N;i++){
			for(j=0;j<M;j++){
				MatQ[j][i] = Mat[j][i];
			}
			for(j=0;j<N;j++){
				MatR[i][j] = zero;
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				T mul = spfuncmul(Mat[j][i], Mat[j][i], BITS);
				norm[i] = spfuncadd(norm[i], mul, BITS);
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = spfuncsqrt(norm[i], BITS);
			for(k=0;k<M;k++){
				MatQ[k][i] = spfuncdiv(MatQ[k][i],MatR[i][i], BITS);
			}
			for(j=i+1;j<M;j++){
				float tmp;
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = spfuncadd(MatR[i][j], spfuncmul(MatQ[k][i], MatQ[k][j], BITS), BITS);
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = spfuncsub(MatQ[k][j], spfuncmul(MatQ[k][i], MatR[i][j], BITS), BITS);
				}
			// update norm: no update for QR_MGS
				// norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
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
	template<class T, int M, int N, T spfunc(double),
			T spfuncmul(T,T), T spfuncadd(T,T),
			T spfuncsub(T,T), T spfuncsqrt(T),
			T spfuncdiv(T,T), bool spfunceq(T, T)>
	void QRD_HH(T **Mat,
			    T **MatQ,
			    T **MatR){
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
		T *x, *v, *w, *u;
		x = (T*) malloc(sizeof(T)*M);
		v = (T*) malloc(sizeof(T)*M);
		w = (T*) malloc(sizeof(T)*M);
		u = (T*) malloc(sizeof(T)*N);
		T **tmp1, **tmp2;
		tmp1 = (T**) malloc(sizeof(T*)*M);
		tmp2 = (T**) malloc(sizeof(T*)*M);
		for(int i=0;i<N;i++){
			tmp1[i] = (T*) malloc(sizeof(T)*N);
			tmp2[i] = (T*) malloc(sizeof(T)*N);
		}
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
	template<class T, class T1, int M, int N, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T,T, int), T spfuncadd(T,T, int),
			T spfuncsub(T,T, int), T spfuncsqrt(T, int),
			T spfuncdiv(T,T, int), bool spfunceq(T, T), int BITS>
	void QRD_HH(T **Mat,
			    T **MatQ,
			    T **MatR){
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
		T *x, *v, *w, *u;
		x = (T*) malloc(sizeof(T)*M);
		v = (T*) malloc(sizeof(T)*M);
		w = (T*) malloc(sizeof(T)*M);
		u = (T*) malloc(sizeof(T)*N);
		T **tmp1, **tmp2;
		tmp1 = (T**) malloc(sizeof(T*)*M);
		tmp2 = (T**) malloc(sizeof(T*)*M);
		for(int i=0;i<N;i++){
			tmp1[i] = (T*) malloc(sizeof(T)*N);
			tmp2[i] = (T*) malloc(sizeof(T)*N);
		}
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
	template<class T, int M, T spfunc(double),
			T spfuncmul(T, T), T spfuncadd(T, T),
			T spfuncsub(T, T), T spfuncdiv(T, T)>
	void UPTRIANGULARMATINV(T **R,T **Ri){
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
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T, T, int), T spfuncadd(T, T, int),
			T spfuncsub(T, T, int), T spfuncdiv(T, T, int), int BITS>
	void UPTRIANGULARMATINV(T **R,T **Ri){

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
		int i=0,j=0;
		// Q inversion
		for(i=0;i<M;i++){
			for(j=0;j<M;j++){
				Qi[i][j] = Q[j][i];
			}
		}
	}
	template<class T, int M>
	void ORTHOGONALMATINV(T **Q,T **Qi){
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
		T Q[M][M], R[M][M];
		T Qi[M][M], Ri[M][M];
		QRD_HH<T, M, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq>(A, Q, R);
		UPTRIANGULARMATINV<T, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncdiv>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, M, M, M, spfunc, spfuncmul, spfuncadd>(Ri, Qi, B);
	}
	template<class T, int M, T spfunc(double),
			T spfuncmul(T,T), T spfuncadd(T,T),
			T spfuncsub(T,T), T spfuncsqrt(T),
			T spfuncdiv(T,T), bool spfunceq(T,T)>
	void MAT_QRINV(T **A, T **B){
		T **Q, **R;
		T **Qi, **Ri;
		Q = (T**) malloc(sizeof(T*)*M);
		R = (T**) malloc(sizeof(T*)*M);
		Qi = (T**) malloc(sizeof(T*)*M);
		Ri = (T**) malloc(sizeof(T*)*M);
		for(int i=0;i<M;i++){
			Q[i] = (T*) malloc(sizeof(T)*M);
			R[i] = (T*) malloc(sizeof(T)*M);
			Qi[i] = (T*) malloc(sizeof(T)*M);
			Ri[i] = (T*) malloc(sizeof(T)*M);
		}
		QRD_HH<T, M, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq>(A, Q, R);
		UPTRIANGULARMATINV<T, M, spfunc, spfuncmul, spfuncadd, spfuncsub, spfuncdiv>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, M, M, M, spfunc, spfuncmul, spfuncadd>(Ri, Qi, B);
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T,T, int), T spfuncadd(T,T, int),
			T spfuncsub(T,T, int), T spfuncsqrt(T, int),
			T spfuncdiv(T,T, int), bool spfunceq(T,T), int BITS>
	void MAT_QRINV(T A[M][M], T B[M][M]){
		T Q[M][M], R[M][M];
		T Qi[M][M], Ri[M][M];
		QRD_HH<T, T1, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq, BITS>(A, Q, R);
		UPTRIANGULARMATINV<T, T1, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncdiv, BITS>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, T1, M, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, BITS>(Ri, Qi, B);
	}
	template<class T, class T1, int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncmul(T,T, int), T spfuncadd(T,T, int),
			T spfuncsub(T,T, int), T spfuncsqrt(T, int),
			T spfuncdiv(T,T, int), bool spfunceq(T,T), int BITS>
	void MAT_QRINV(T **A, T **B){
		T **Q, **R;
		T **Qi, **Ri;
		Q = (T**) malloc(sizeof(T*)*M);
		R = (T**) malloc(sizeof(T*)*M);
		Qi = (T**) malloc(sizeof(T*)*M);
		Ri = (T**) malloc(sizeof(T*)*M);
		for(int i=0;i<M;i++){
			Q[i] = (T*) malloc(sizeof(T)*M);
			R[i] = (T*) malloc(sizeof(T)*M);
			Qi[i] = (T*) malloc(sizeof(T)*M);
			Ri[i] = (T*) malloc(sizeof(T)*M);
		}
		QRD_HH<T, T1, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncsqrt, spfuncdiv, spfunceq, BITS>(A, Q, R);
		UPTRIANGULARMATINV<T, T1, M, spfunc, spfunc1, spfuncmul, spfuncadd, spfuncsub, spfuncdiv, BITS>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, T1, M, M, M, spfunc, spfunc1, spfuncmul, spfuncadd, BITS>(Ri, Qi, B);
	}


	// Cholesky–Banachiewicz (row based)
	// and Cholesky–Crout (column based)
	// LU decomposition
	template<class T, unsigned int M, T spfunc(double),
			 T spfuncadd(T, T), T spfuncsub(T, T),
			 T spfuncmul(T, T), T spfuncdiv(T, T),
			 T spfuncsqrt(T), bool spfunclt(T, T)>
	void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

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
				if(spfunclt(LL,spfunc(0))){
					exit(-1);
				}
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
	template<class T, unsigned int M, T spfunc(double),
			 T spfuncadd(T, T), T spfuncsub(T, T),
			 T spfuncmul(T, T), T spfuncdiv(T, T),
			 T spfuncsqrt(T), bool spfunclt(T, T)>
	void LU_CHOLBANACHROUT(T **Mat, T **MatL, T **MatU){

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
				if(spfunclt(LL,spfunc(0))){
					exit(-1);
				}
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
				if(spfunclt(LL,zero)){
					exit(-1);
				}
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
	template<class T, class T1, unsigned int M, T1 spfunc(double, int), T spfunc1(T1, int),
			 T spfuncadd(T, T, int), T spfuncsub(T, T, int),
			 T spfuncmul(T, T, int), T spfuncdiv(T, T, int),
			 T spfuncsqrt(T, int), bool spfunclt(T, T), int BITS>
	void LU_CHOLBANACHROUT(T **Mat, T **MatL, T **MatU){

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
				if(spfunclt(LL,zero)){
					exit(-1);
				}
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


	// Doolittle algorithm LU decomposition
	template<class T, unsigned int M, T spfunc(double), T spfuncsub(T, T),
			 T spfuncmul(T, T), T spfuncdiv(T, T)>
	void LU_DOOLITTLE(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

		// clean the matrix
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatL[i][j] = spfunc(0);
				MatU[i][j] = spfunc(0);
			}
		}

		// decomposition
		for(unsigned int i=0;i<M;i++){
			// upper triangular
			for(unsigned int k=0;k<M;k++){
				T tmp = Mat[i][k];
				for(unsigned int j=0;j<i;j++){
					T tmp1 = spfuncmul(MatL[i][j],MatU[j][k]);
					tmp = spfuncsub(tmp, tmp1);
				}
				MatU[i][k] = tmp;
			}
			// lower triangular
			MatL[i][i] = spfunc(1);
			for(unsigned int k=i+1;k<M;k++){
				T tmp = Mat[k][i];
				for(unsigned int j=0;j<i;j++){
					T tmp1 = spfuncmul(MatL[k][j],MatU[j][i]);
					tmp = spfuncsub(tmp, tmp1);
				}
				MatL[k][i] = spfuncdiv(tmp, MatU[i][i]);
			}
		}
	}
	template<class T, class T1, unsigned int M, T1 spfunc(double, int), T spfunc1(T1, int),
			T spfuncsub(T, T, int), T spfuncmul(T, T, int), T spfuncdiv(T, T, int), int BITS>
	void LU_DOOLITTLE(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

		T zero = spfunc1(spfunc(0, BITS), BITS);
		T one = spfunc1(spfunc(1, BITS), BITS);

		// clean the matrix
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatL[i][j] = zero;
				MatU[i][j] = zero;
			}
		}

		// decomposition
		for(unsigned int i=0;i<M;i++){
			// upper triangular
			for(unsigned int k=0;k<M;k++){
				T tmp = Mat[i][k];
				for(unsigned int j=0;j<i;j++){
					T tmp1 = spfuncmul(MatL[i][j],MatU[j][k], BITS);
					tmp = spfuncsub(tmp, tmp1, BITS);
				}
				MatU[i][k] = tmp;
			}
			// lower triangular
			MatL[i][i] = one;
			for(unsigned int k=i+1;k<M;k++){
				T tmp = Mat[k][i];
				for(unsigned int j=0;j<i;j++){
					T tmp1 = spfuncmul(MatL[k][j],MatU[j][i], BITS);
					tmp = spfuncsub(tmp, tmp1, BITS);
				}
				MatL[k][i] = spfuncdiv(tmp, MatU[i][i], BITS);
			}
		}
	}


	// Extract submatrix
	// input matrix must larger than output matrix
	// output matrix (N-M)*(L-K) begins at row M, col K
	template<class T, int R, int C,
			 int M, int N>
	void MAT_SUBMAT(T MatA[R][C], T MatB[M][N],
					int row, int col){
		for(int i=row, k=0;i<M+row;i++,k++){
			for(int j=col, l=0;j<N+col;j++,l++){
				MatB[k][l] = MatA[i][j];
			}
		}
	}
	template<class T, int R, int C,
			 int M, int N>
	void MAT_SUBMAT(T MatA[R][C], T MatB[M][N],
					int row, int col,
					int step_r, int step_c){
		for(int i=row, k=0;i<R&&k<M;i+=step_r,k++){
			for(int j=col, l=0;j<C&&l<N;j+=step_c,l++){
				MatB[k][l] = MatA[i][j];
			}
		}
	}

	// Matrix merge
	// MatA: (M,N), MatB (M,K), MatC: (M, N+K)
	template<class T, int M, int N, int K>
	void MAT_MERGE_H(T MatA[M][N],
					 T MatB[M][K],
					 T MatC[M][N+K]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j];
		for(int i=0;i<M;i++)
			for(int j=N,k=0;j<N+K&&k<K;j++,k++)
				MatC[i][j] = MatB[i][k];
	}
	// MatA: (M,K), MatB: (N,K), MatC: (M+N,K)
	template<class T, int M, int N, int K>
	void MAT_MERGE_V(T MatA[M][K],
					 T MatB[N][K],
					 T MatC[M+N][K]){
		for(int i=0;i<M;i++)
			for(int j=0;j<K;j++)
				MatC[i][j] = MatA[i][j];
		for(int i=M,k=0;i<M+K&&k<N;i++,k++)
			for(int j=0;j<K;j++)
				MatC[i][j] = MatB[k][j];
	}

	// Matrix Differential
	template<class T,
			 int M, int N,
			 int dim, int order,
			 T spfuncmul(T, T), T spfuncsub(T, T)>
	void MAT_DIFF(T MatA[M][N],
				  T MatB[M][N],
				  int dir){
		T MatTmp1[M][N];
		T MatTmp2[M][N];
		MAT_EQ<T, M, N>(MatA, MatTmp1);
		if(dim==1){
			for(int i=0;i<order;i++){
				for(int j=0;j<N;j++){
					for(int k=0;k<M-i-1;k++){
						MatTmp2[k][j] = spfuncsub(MatTmp1[k+1][j], MatTmp1[k][j]);
						MatTmp2[k][j] = spfuncmul(MatTmp2[k][j], dir);
					}
				}
				MAT_EQ<T, M, N>(MatTmp2, MatTmp1);
			}
			MAT_EQ<T, M, N>(MatTmp1, MatB);
		}else if(dim==2){
			for(int i=0;i<order;i++){
				for(int j=0;j<M;j++){
					for(int k=0;k<N-i-1;k++){
						MatTmp2[j][k] = spfuncsub(MatTmp1[j][k+1], MatTmp1[j][k]);
						MatTmp2[j][k] = spfuncmul(MatTmp2[j][k], dir);
					}
				}
				MAT_EQ<T, M, N>(MatTmp2, MatTmp1);
			}
			MAT_EQ<T, M, N>(MatTmp1, MatB);
		}else{
			std::cout << "only support 2 dimension !!" << std::endl;
			std::exit(0);
		}
	}

	// From real matrix to complex matrix
	template<class T, int M, int N>
	void MAT_REAL2COMPLEX(T MatA[M][N],
				 	 	  Complex<T> MatB[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatB[i][j].real = MatA[i][j];
				MatB[i][j].imag = 0;
			}
		}
	}

	// Extract real matrix from complex matrix
	template<class T, int M, int N>
	void MAT_COMPLEX_GETREAL(Complex<T> MatA[M][N],
							 T MatB[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatB[i][j] = MatA[i][j].real;
			}
		}
	}

	// Extract real submatrix from complex matrix
	template<class T, int M, int N, int P, int Q>
	void MAT_COMPLEX_GETREALSUBMAT(Complex<T> MatA[M][N],
							 T MatB[P][Q], int row, int col){
		for(int i=row;i<P+row;i++){
			for(int j=col;j<Q+col;j++){
				MatB[i][j] = MatA[i][j].real;
			}
		}
	}

	// Extract imag matrix from complex matrix
	template<class T, int M, int N>
	void MAT_COMPLEX_GETIMAG(Complex<T> MatA[M][N],
							 T MatB[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatB[i][j] = MatA[i][j].imag;
			}
		}
	}

	// Extract imag submatrix from complex matrix
	template<class T, int M, int N, int P, int Q>
	void MAT_COMPLEX_GETIMAGSUBMAT(Complex<T> MatA[M][N],
							 T MatB[P][Q], int row, int col){
		for(int i=row;i<P+row;i++){
			for(int j=col;j<Q+col;j++){
				MatB[i][j] = MatA[i][j].imag;
			}
		}
	}

	// Matrix compare scalar
	template<class T, int M, int N,
			 bool spfunclt(T, T)>
	void MAT_MAXCMP(T MatA[M][N],
					T B,
				    T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = !spfunclt(MatA[i][j], B) ? MatA[i][j] : B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_MAXCMP(T1 MatA[M][N],
					T2 B,
				    T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j]>=B ? MatA[i][j] : B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_MAXCMP(T2 B,
					T1 MatA[M][N],
				    T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j]>=B ? MatA[i][j] : B;
	}
	template<class T, int M, int N,
			 bool spfunclt(T, T)>
	void MAT_MINCMP(T MatA[M][N],
				    T B,
				    T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = !spfunclt(MatA[i][j], B) ? B : MatA[i][j];
	}
	template<class T1, class T2, int M, int N>
	void MAT_MINCMP(T1 MatA[M][N],
				    T2 B,
				    T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j]>=B ? B : MatA[i][j];
	}
	template<class T1, class T2, int M, int N>
	void MAT_MINCMP(T2 B,
					T1 MatA[M][N],
				    T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j]>=B ? B : MatA[i][j];
	}

	// Matrix sign
	template<class T, int M, int N, T spfunc(double),
			 bool spfunclt(T, T)>
	void MAT_SIGN(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfunclt(MatA[i][j], spfunc(0)) ?
						     1 : -1;
	}

	// Matrix absolute value
	template<class T, int M, int N, T spfunc(double),
		     bool spfunclt(T, T)>
	void MAT_ABS(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfunclt(MatA[i][j], spfunc(0)) ?
						     MatA[i][j] : -MatA[i][j];
	}

	// Matrix dot multiplication
	template<class T, int M, int N,
			 T spfuncmul(T, T)>
	void MAT_DOTSQUARE(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfuncmul(MatA[i][j], MatA[i][j]);
	}
	template<class T, int M, int N,
			 T spfuncmul(T, T, int), int BITS>
	void MAT_DOTSQUARE(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfuncmul(MatA[i][j], MatA[i][j], BITS);
	}

	// Matrix dot multiplication
	template<class T, int M, int N,
			 T spfuncmul(T, T)>
	void MAT_DOTMUL(T MatA[M][N],
				    T MatB[M][N],
					T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], MatB[i][j]);
	}
	template<class T, int M, int N,
			 T spfuncmul(T, T, int), int BITS>
	void MAT_DOTMUL(T MatA[M][N],
				    T MatB[M][N],
					T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], MatB[i][j], BITS);
	}

	// Matrix dot division
	template<class T, int M, int N>
	void MAT_DOTDIV(T MatA[M][N],
				    T MatB[M][N],
					T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] / MatB[i][j];
	}

	// Matrix dot division
	template<class T, int M, int N,
	 	 	 T spfuncdiv(T, T)>
	void MAT_COMPLEX_DOTDIV_REAL(Complex<T> MatA[M][N],
							 	 T MatB[M][N],
								 Complex<T> MatC[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatC[i][j].real = spfuncdiv(MatA[i][j].real, MatB[i][j]);
				MatC[i][j].imag = spfuncdiv(MatA[i][j].imag, MatB[i][j]);
			}
		}
	}
	template<class T, int M, int N,
	 	 	 T spfuncdiv(T, T, int), int BITS>
	void MAT_COMPLEX_DOTDIV_REAL(Complex<T> MatA[M][N],
							 	 T MatB[M][N],
								 Complex<T> MatC[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatC[i][j].real = spfuncdiv(MatA[i][j].real, MatB[i][j], BITS);
				MatC[i][j].imag = spfuncdiv(MatA[i][j].imag, MatB[i][j], BITS);
			}
		}
	}
	template<class T, int M, int N, int P, int Q>
	void MAT_COMPLEX_DOTDIV_REAL(Complex<T> MatA[M][N],
							 	 T MatB[P][Q],
								 Complex<T> MatC[M][N]){
		for(int i=0;i<P;i++){
			for(int j=0;j<Q;j++){
				MatC[i][j].real = MatA[i][j].real / MatB[i][j];
				MatC[i][j].imag = MatA[i][j].imag / MatB[i][j];
			}
		}
	}

	// Matrix dot inverse
	template<class T, int M, int N,
 	 	 	 T spfuncdiv(T, T)>
	void MAT_DOTINV(T MatA[M][N],
				    T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfuncdiv(1, MatA[i][j]);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncdiv(T, T, int), int BITS>
	void MAT_DOTINV(T MatA[M][N],
				    T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = spfuncdiv(1, MatA[i][j], BITS);
	}

	// Matrix dot multiply scalar
	template<class T, int M, int N,
	 	 	 T spfuncadd(T, T)>
	void MAT_SCALAR_DOTADD(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncadd(MatA[i][j], B);
	}
	template<class T, int M, int N,
	 	 	 T spfuncadd(T, T, int), int BITS>
	void MAT_SCALAR_DOTADD(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncadd(MatA[i][j], B, BITS);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncadd(T, T)>
	void MAT_SCALAR_DOTADD(T B,
						   T MatA[M][N],
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncadd(MatA[i][j], B);
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTADD(T1 MatA[M][N],
						   T2 B,
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] + B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTADD(T2 B,
						   T1 MatA[M][N],
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] + B;
	}

	// Matrix dot multiply scalar
	template<class T, int M, int N,
	 	 	 T spfuncsub(T, T)>
	void MAT_SCALAR_DOTSUB(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncsub(MatA[i][j], B);
	}
	template<class T, int M, int N,
	 	 	 T spfuncsub(T, T, int), int BITS>
	void MAT_SCALAR_DOTSUB(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncsub(MatA[i][j], B, BITS);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncsub(T, T)>
	void MAT_SCALAR_DOTSUB(T B,
						   T MatA[M][N],
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncsub(B, MatA[i][j]);
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTSUB(T1 MatA[M][N],
						   T2 B,
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] - B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTSUB(T2 B,
						   T1 MatA[M][N],
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = B - MatA[i][j];
	}

	// Matrix dot multiply scalar
	template<class T, int M, int N,
 	 	 	 T spfuncmul(T, T)>
	void MAT_SCALAR_DOTMUL(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], B);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncmul(T, T)>
	void MAT_SCALAR_DOTMUL(T **MatA,
						   T B,
						   T **MatC){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], B);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncmul(T, T, int), int BITS>
	void MAT_SCALAR_DOTMUL(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], B, BITS);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncmul(T, T, int), int BITS>
	void MAT_SCALAR_DOTMUL(T **MatA,
						   T B,
						   T **MatC){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], B, BITS);
	}
	template<class T, int M, int N,
	 	 	 T spfuncmul(T, T)>
	void MAT_SCALAR_DOTMUL(T B,
						   T MatA[M][N],
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncmul(MatA[i][j], B);
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTMUL(T1 MatA[M][N],
						   T2 B,
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] * B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTMUL(T2 B,
						   T1 MatA[M][N],
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] * B;
	}

	// Matrix dot division scalar
	template<class T, int M, int N,
 	 	 	 T spfuncdiv(T, T)>
	void MAT_SCALAR_DOTDIV(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncdiv(MatA[i][j], B);
	}
	template<class T, int M, int N,
 	 	 	 T spfuncdiv(T, T, int), int BITS>
	void MAT_SCALAR_DOTDIV(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncdiv(MatA[i][j], B, BITS);
	}
	template<class T, int M, int N,
	 	 	 T spfuncdiv(T, T)>
	void MAT_SCALAR_DOTDIV(T B,
						   T MatA[M][N],
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = spfuncdiv(B, MatA[i][j]);
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTDIV(T1 MatA[M][N],
						   T2 B,
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] / B;
	}
	template<class T1, class T2, int M, int N>
	void MAT_SCALAR_DOTDIV(T2 B,
						   T1 MatA[M][N],
						   T1 MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = B / MatA[i][j];
	}

	// Extract subvector
	template<class T, int M, int N>
	void VEC_SUBVEC(T VecA[M], T VecB[N], int start){
		for(int i=start,k=0;i<M&&k<N;i++,k++)
			VecB[k] = VecA[i];
	}
	template<class T, int M, int N>
	void VEC_SUBVEC(T VecA[M], T VecB[N],
					int start, int step){
		for(int i=start,k=0;i<M&&k<N;i+=step,k++)
			VecB[k] = VecA[i];
	}

	// Vector merge
	template<class T, int M, int N>
	void VEC_MERGE2VEC(T VecA[M],
				   T VecB[N],
				   T VecC[M+N]){
		for(int i=0;i<M;i++)
			VecC[i] = VecA[i];
		for(int i=M,k=0;i<N+M&&k<N;i++,k++)
			VecC[i] = VecB[i];

	}
	template<class T, int M>
	void VEC_MERGE2MAT(T VecA[M],
			   	   	   T VecB[M],
					   T MatC[2][M]){
		for(int i=0;i<M;i++){
			MatC[0][i] = VecA[i];
			MatC[1][i] = VecB[i];
		}
	}

	// Vector differential
	template<class T, int M,
			T spfuncsub(T, T)>
	void VEC_DIFF(T VecA[M],
				  T VecB[M],
				  int order){
		T VecTmp[M];
		VEC_EQ<T>(VecA, VecTmp);
		for(int i=0;i<order;i++){
			for(int j=0;j<M-i-1;j++)
				VecB[i] = spfuncsub(VecTmp[i+1], VecTmp[i]);
			VEC_EQ<T>(VecB, VecTmp);
		}
	}
	template<class T, int M,
			T spfuncsub(T, T, int), int BITS>
	void VEC_DIFF(T VecA[M],
				  T VecB[M],
				  int order){
		T VecTmp[M];
		VEC_EQ<T>(VecA, VecTmp);
		for(int i=0;i<order;i++){
			for(int j=0;j<M-i-1;j++)
				VecB[i] = spfuncsub(VecTmp[i+1], VecTmp[i], BITS);
			VEC_EQ<T>(VecB, VecTmp);
		}
	}

	// Vector absolute value
	template<class T, int M, T spfunc(double),
			 bool spfunclt(T, T)>
	void VEC_ABS(T VecA[M],
				 T VecB[M]){
		for(int i=0;i<M;i++)
			VecB[i] = !spfunclt(VecA[i], spfunc(0)) ? VecA[i] : -VecA[i];
	}

	// Vector dot division
	template<class T, int M,
			 T spfuncdiv(T, T)>
	void VEC_DOTDIV(T MatA[M],
				    T MatB[M],
					T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = spfuncdiv(MatA[i], MatB[i]);
	}
	template<class T, int M,
			 T spfuncdiv(T, T, int), int BITS>
	void VEC_DOTDIV(T MatA[M],
				    T MatB[M],
					T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = spfuncdiv(MatA[i], MatB[i], BITS);
	}

	// Vector dot division scalar
	template<class EigenT, class T, int M,
	 	 	 T spfuncdiv(T, T)>
	void VEC_SCALAR_DOTDIV(T MatA[M],
						   T B,
						   T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = spfuncdiv(MatA[i], B);
	}
	template<class EigenT, class T, int M,
	 	 	 T spfuncdiv(T, T, int), int BITS>
	void VEC_SCALAR_DOTDIV(T MatA[M],
						   T B,
						   T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = spfuncdiv(MatA[i], B, BITS);
	}
	template<class EigenT, class T, int M,
	 	 	 T spfuncdiv(T, T)>
	void VEC_SCALAR_DOTDIV(T A,
						   T MatB[M],
						   T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = spfuncdiv(A, MatB[i]);
	}


private:

protected:
};


#endif /* SRC_SOFTPOSIT_ALGEBRA_HPP_ */
