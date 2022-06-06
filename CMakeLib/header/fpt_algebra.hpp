/*
 *	This is light floating point algebra library
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_FPT_ALGEBRA_HPP_
#define SRC_FPT_ALGEBRA_HPP_

#include "data.hpp"
#include "common.hpp"

class Float_Point_Algebra{

public:

	// Generate all zero matrix
	template<class T, int M, int N>
	void ZEROS_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = 0;
			}
		}
	}

	// Generate all one matrix
	template<class T, int M, int N>
	void ONES_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				A[i][j] = 1;
			}
		}
	}

	// Generate identity matrix
	template<class T, int M, int N>
	void IDENDTITY_MAT( T A[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				if(i==j)
					A[i][j] = 1;
				else
					A[i][j] = 0;
			}
		}
	}
	template<class T, int M, int N>
	void IDENDTITY_MAT( T **A ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				if(i==j)
					A[i][j] = 1;
				else
					A[i][j] = 0;
			}
		}
	}

	// Generate random matrix
	template<class T, int M, int N>
	void RND_MAT( T A[M][N] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
				double rndnum = (double) INTEGER_SCALE * rnd ;
				A[i][j] = (T) rndnum;
			}
		}
	}

	// Generate random matrix
	template<class T, int M, int N>
	void RND_DIAGMAT( T A[M][N] ){
		srand (time(NULL));
		int sparse_num = floor( DIAG_RATIO * M );
		for( int i=0; i<M; i++){
			for( int j=0; j<N; j++ ){
				if( i<sparse_num && i==j){
					double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
					double rndnum = (double) INTEGER_SCALE * rnd ;
					A[i][j] = (T) rndnum;
				}else
					A[i][j] = 0;
			}
		}
	}

	// Generate all zero vector
	template<class T, int M>
	void ZEROS_VEC( T V[M] ){
		for( int i=0; i<M; i++ ){
			V[i] = 0;
		}
	}

	// Generate all one vector
	template<class T, int M>
	void ONES_VEC( T V[M] ){
		for( int i=0; i<M; i++ ){
			V[i] = 1;
		}
	}

	// Generate random vector
	template<class T, int M>
	void RND_VEC( T V[M] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
			double rndnum = (double) INTEGER_SCALE * rnd ;
			V[i] = (T) rndnum;
		}
	}
	template<class T, int M>
	void RND_VEC_SCALE( T V[M], T scale ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
			double rndnum = (double) scale * INTEGER_SCALE * rnd ;
			V[i] = (T) rndnum;
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
	template<class T, int M>
	void VEC_MINUS( T V1[M], T V2[M] ){
		for( int i=0; i<M; i++ ){
			V2[i] = -V1[i];
		}
	}

	// Vector addition
	template<class T, int M>
	void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] + V2[i];
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_ADD( T1 V1[M], T2 V2[M], T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] + V2[i];
		}
	}

	// Vector subtraction
	template<class T, int M>
	void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] - V2[i];
		}
	}

	// Vector multiplication
	template<class T, int M>
	void VEC_MUL_2SCALAR( T V1[M], T V2[M], T S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V2[i];
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_MUL_2SCALAR( T1 V1[M], T2 V2[M], T3 S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V2[i];
		}
	}
	template<class T, int M, int N>
	void VEC_MUL_2MATRIX( T V1[M], T V2[N], T M3[M][N] ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				M3[i][j] = V1[i] * V2[j];
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
	template<class T, int M>
	void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] / V2[i];
		}
	}
	template<class T, int M>
	void VEC_DIV( T V1[M], T v, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] / v;
		}
	}

	// Vector norm
	template<class T, int M>
	void VEC_NORM( T V1[M], T &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
		double tmp = (double) S;
		S = (T) std::sqrt(tmp);
	}
	template<class T1, class T2, int M>
	void VEC_NORM( T1 V1[M], T2 &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
		double tmp = (double) S;
		S = (T2) std::sqrt(tmp);
	}

	// Vector norm .^2
	template<class T, int M>
	void VEC_NORM2( T V1[M], T &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
	}
	template<class T1, class T2, int M>
	void VEC_NORM2( T1 V1[M], T2 &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
		}
	}

	// Vector square norm
	template<class T, int M>
	void VEC_SQUARE_NORM( T V1[M], T &S ){
		S = 0;
		for( int i=0; i<M; i++ ){
			S += V1[i] * V1[i];
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
	template<class T, int M>
	void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] + S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_ADD( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] + S;
		}
	}

	// Vector scalar subtraction
	template<class T, int M>
	void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] - S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_SUB( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] - S;
		}
	}

	// Vector scalar multiplication
	template<class T, int M>
	void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] * S;
		}
	}
	template<class T, int M>
	void VEC_SCALAR_MUL( T S, T V1[M], T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] * S;
		}
	}

	// Vector scalar compare
	template<class T, int M>
	void VEC_SCALAR_MIN( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] < S ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MIN( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] < S ? V1[i] : S;
		}
	}
	template<class T, int M>
	void VEC_SCALAR_MAX( T V1[M], T S, T V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] > S ? V1[i] : S;
		}
	}
	template<class T1, class T2, class T3, int M>
	void VEC_SCALAR_MAX( T1 V1[M], T2 S, T3 V3[M] ){
		for( int i=0; i<M; i++ ){
			V3[i] = V1[i] > S ? V1[i] : S;
		}
	}

	// The basic matrix addition with complexity of O(NN)
	template<class T, int M, int N>
	void MAT_ADD(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] + B[i][j];
			}
		}
	}
	template<class T, int M, int N>
	void MAT_ADD(T **A,
									      T **B,
									      T **C){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] + B[i][j];
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
	template<class T, int M, int N>
	void MAT_SUB(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] - B[i][j];
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

	// Matrix Vector multiplication
	template<class T, int M, int N>
	void MAT_VEC_MUL(T A[M][N],
										  T B[N],
										  T C[M]){
		for ( int i=0; i<M; i++ ){
			C[i] = 0;
			for ( int j=0; j<N; j++ ){
				C[i] += A[i][j] * B[j];
			}
		}
	}
	template<class T, int M, int N>
	void MAT_VEC_MUL(T **A,
										  T *B,
										  T *C){
//	   	std::cout << __FILE__ << __LINE__ << std::endl;
//		printf("%lf\n",A[0][0]);
		for ( int i=0; i<M; i++ ){
			C[i] = 0;
//		   	std::cout << __FILE__ << __LINE__ << std::endl;
			for ( int j=0; j<N; j++ ){
				C[i] += A[i][j] * B[j];
//			   	std::cout << __FILE__ << __LINE__ << std::endl;
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
	template<class T, int M, int N>
	void VEC_MAT_MUL(T A[M],
										  T B[M][N],
										  T C[N]){
		for ( int i=0; i<N; i++ ){
			C[i] = 0;
			for ( int j=0; j<M; j++ ){
				C[i] += A[j] * B[j][i];
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
	template<class T, int M, int N, int P>
	void MAT_MUL(T A[M][N],
												    T B[N][P],
													T C[M][P]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<P; j++ ){
				C[i][j] = 0;
				for ( int k=0; k<N; k++ ){
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
	}
	template<class T, int M, int N, int P>
	void MAT_MUL(T **A,
												    T **B,
													T **C){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<P; j++ ){
				C[i][j] = 0;
				for ( int k=0; k<N; k++ ){
					C[i][j] += A[i][k] * B[k][j];
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
	template<class T, int M>
	void SQUARE_MAT_MUL(T A[M][M],
												    T B[M][M],
													T C[M][M]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<M; j++ ){
				C[i][j] = 0;
				for ( int k=0; k<M; k++ ){
					C[i][j] += A[i][k] * B[k][j];
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
	template<class T, int M, int N>
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
				MatR[i][j] = 0;
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				float mul = Mat[j][i] * Mat[j][i];
				norm[i] = norm[i] + mul;
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = (T) std::sqrt(norm[i]);
			for(k=0;k<M;k++){
				MatQ[k][i]=MatQ[k][i]/MatR[i][i];
			}
			for(j=i+1;j<M;j++){
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = MatR[i][j] + MatQ[k][i] * MatQ[k][j];
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = MatQ[k][j] - MatQ[k][i] * MatR[i][j];
				}
				// update norm
				norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
			}
		}
	}
	template<class T, int M, int N>
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
				MatR[i][j] = 0;
			}
		}

		// Phase 1: get the norm
		float norm[N];
		for(i=0;i<N;i++){
			norm[i] = 0;
			for(j=0;j<M;j++){
				float mul = Mat[j][i] * Mat[j][i];
				norm[i] = norm[i] + mul;
			}
		}

		// Phase 2: get the Q&R
		for(i=0;i<N;i++){
			// derive R
			MatR[i][i] = (T) std::sqrt(norm[i]);
			for(k=0;k<M;k++){
				MatQ[k][i]=MatQ[k][i]/MatR[i][i];
			}
			for(j=i+1;j<M;j++){
				float tmp;
				for(k=0;k<M;k++){
					// update R
					MatR[i][j] = MatR[i][j] + MatQ[k][i] * MatQ[k][j];
				}
				for(k=0;k<M;k++){
					// update Q
					MatQ[k][j] = MatQ[k][j] - MatQ[k][i] * MatR[i][j];
				}
			// update norm: no update for QR_MGS
				// norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
			}
		}

	}
	template<class T, int M, int N>
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
					MatQ[j][i] = (T) 1;
				else
					MatQ[j][i] = (T) 0;
			}
		}

		T g, s;
		T x[M], v[M], w[M], u[N];
		T tmp1[M][N], tmp2[M][N];
		for(k=0;k<M-1;k++){
			// x=zeros(m,1);
			for(j=0;j<M;j++){
				x[j] = (T) 0;
			}
			ZEROS_VEC<T, M>(x);
			//x(k:m,1)=R(k:m,k);
			for(j=k;j<M;j++){
				x[j] = MatR[j][k];
			}
			//g=norm(x);
			VEC_NORM<T, M>(x, g);
			// v=x;
			VEC_EQ<T, M>(x, v);
			// v(k)=x(k)+g;
			v[k] = x[k] + g;
			//s=norm(v);
			VEC_NORM<T, M>(v, s);
			if(s!=0){
				// w=v/s;
				VEC_DIV<T, M>(v, s, w);
				// u=2*R'*w;
				for(i=0;i<N;i++){
					u[i] = 0;
					for(j=0;j<M;j++){
						u[i] += 2 * MatR[j][i] * w[j];
					}
				}
				// R=R-w*u'; %Product HR
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						MatR[j][i] = MatR[j][i] - w[j] * u[i];
					}
				}
				// Q=Q-2*Q*w*w'; %Product QR
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						tmp1[j][i] = w[j] * w[i];
					}
				}
				MAT_MUL<T, M, N, N>(MatQ, tmp1, tmp2);
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						MatQ[j][i] = MatQ[j][i] - 2 * tmp2[j][i];
					}
				}
			}
		}
	}
	template<class T, int M, int N>
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
					MatQ[j][i] = (T) 1;
				else
					MatQ[j][i] = (T) 0;
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
				x[j] = (T) 0;
			}
			ZEROS_VEC<T, M>(x);
			//x(k:m,1)=R(k:m,k);
			for(j=k;j<M;j++){
				x[j] = MatR[j][k];
			}
			//g=norm(x);
			VEC_NORM<T, M>(x, g);
			// v=x;
			VEC_EQ<T, M>(x, v);
			// v(k)=x(k)+g;
			v[k] = x[k] + g;
			//s=norm(v);
			VEC_NORM<T, M>(v, s);
			if(s!=0){
				// w=v/s;
				VEC_DIV<T, M>(v, s, w);
				// u=2*R'*w;
				for(i=0;i<N;i++){
					u[i] = 0;
					for(j=0;j<M;j++){
						u[i] += 2 * MatR[j][i] * w[j];
					}
				}
				// R=R-w*u'; %Product HR
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						MatR[j][i] = MatR[j][i] - w[j] * u[i];
					}
				}
				// Q=Q-2*Q*w*w'; %Product QR
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						tmp1[j][i] = w[j] * w[i];
					}
				}
				MAT_MUL<T, M, N, N>(MatQ, tmp1, tmp2);
				for(j=0;j<M;j++){
					for(i=0;i<N;i++){
						MatQ[j][i] = MatQ[j][i] - 2 * tmp2[j][i];
					}
				}
			}
		}
	}

	template<class T, int M>
	void UPTRIANGULARMATINV(T R[M][M],T Ri[M][M]){
		int i=0,j=0,k=0;
		// R inversion
		for(i=0;i<M;i++){
			for(j=0;j<M;j++){
				Ri[i][j]=0;
			}
		}
		for(i=0;i<M;i++){
			Ri[i][i]=1/R[i][i];
			for(j=i+1;j<M;j++){
				for(k=0;k<=j-1;k++){
					Ri[i][j] = Ri[i][j] - Ri[i][k] * R[k][j] / R[j][j];
				}
			}
		}
	}
	template<class T, int M>
	void UPTRIANGULARMATINV(T **R,T **Ri){
		int i=0,j=0,k=0;
		// R inversion
		for(i=0;i<M;i++){
			for(j=0;j<M;j++){
				Ri[i][j]=0;
			}
		}
		for(i=0;i<M;i++){
			Ri[i][i]=1/R[i][i];
			for(j=i+1;j<M;j++){
				for(k=0;k<=j-1;k++){
					Ri[i][j] = Ri[i][j] - Ri[i][k] * R[k][j] / R[j][j];
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
	template<class T, int M>
	void MAT_QRINV(T A[M][M], T B[M][M]){
		T Q[M][M], R[M][M];
		T Qi[M][M], Ri[M][M];
		QRD_HH<T, M, M>(A, Q, R);
		UPTRIANGULARMATINV<T, M>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, M, M, M>(Ri, Qi, B);
	}
	template<class T, int M>
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
		QRD_HH<T, M, M>(A, Q, R);
		UPTRIANGULARMATINV<T, M>(R, Ri);
		ORTHOGONALMATINV<T, M>(Q, Qi);
		MAT_MUL<T, M, M, M>(Ri, Qi, B);
	}

	// Cholesky–Banachiewicz (row based)
	// and Cholesky–Crout (column based)
	// LU decomposition
	template<class T, unsigned int M>
	void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

		// copy the matrix
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatL[i][j] = 0;
			}
		}

		// decomposition in-place
		for(unsigned int j=0;j<M;j++){
			// compute the diagonal element
			T LL = Mat[j][j];
//			std::cout << j << "th diag update:" << std::endl;
//			std::cout << LL;
			for(unsigned int k=0;k<j;k++){
//				std::cout << "-"
//						  << MatL[j][k]
//					      << "^2";
				LL -= MatL[j][k] * MatL[j][k];
				if(LL<0){
//					std::cout << "=" LL << std::endl << std::endl;
//					std::cout << "something is wrong at ["
//							  << j << "]["
//							  << k << "]=";
					exit(-1);
				}
			}
			MatL[j][j] = (T) std::sqrt((double)LL);
//			std::cout << "=" << MatL[j][j] << std::endl << std::endl;

			// compute the non-diagonal element
//			std::cout << j << "th non-diag update:" << std::endl;
			T inv = 1 / MatL[j][j];
			for(unsigned int i=j+1;i<M;i++){
				LL = Mat[i][j];
//				std::cout << LL;
				for(unsigned int k=0;k<j;k++){
//					std::cout << "-"
//							  << MatL[i][k]
//						      << "x"
//							  << MatL[j][k];
					LL -= MatL[i][k] * MatL[j][k];
				}
				MatL[i][j] = LL * inv;
//				std::cout << "x" << inv << "=" << MatL[i][j] << std::endl;
			}
//			std::cout << std::endl;

//			std::cout << "L Matrix: " << std::endl;;
//			for(int i=0;i<4;i++){
//				for(int j=0;j<4;j++){
//					std::cout << MatL[i][j] << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
		}

		// transpose L to get U
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatU[i][j] = MatL[j][i];
			}
		}
	}
	template<class T, unsigned int M>
	void LU_CHOLBANACHROUT(T **Mat, T **MatL, T **MatU){

		// copy the matrix
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatL[i][j] = 0;
			}
		}

		// decomposition in-place
		for(unsigned int j=0;j<M;j++){
			// compute the diagonal element
			T LL = Mat[j][j];
//			std::cout << j << "th diag update:" << std::endl;
//			std::cout << LL;
			for(unsigned int k=0;k<j;k++){
//				std::cout << "-"
//						  << MatL[j][k]
//					      << "^2";
				LL -= MatL[j][k] * MatL[j][k];
				if(LL<0){
//					std::cout << "=" LL << std::endl << std::endl;
//					std::cout << "something is wrong at ["
//							  << j << "]["
//							  << k << "]=";
					exit(-1);
				}
			}
			MatL[j][j] = (T) std::sqrt((double)LL);
//			std::cout << "=" << MatL[j][j] << std::endl << std::endl;

			// compute the non-diagonal element
//			std::cout << j << "th non-diag update:" << std::endl;
			T inv = 1 / MatL[j][j];
			for(unsigned int i=j+1;i<M;i++){
				LL = Mat[i][j];
//				std::cout << LL;
				for(unsigned int k=0;k<j;k++){
//					std::cout << "-"
//							  << MatL[i][k]
//						      << "x"
//							  << MatL[j][k];
					LL -= MatL[i][k] * MatL[j][k];
				}
				MatL[i][j] = LL * inv;
//				std::cout << "x" << inv << "=" << MatL[i][j] << std::endl;
			}
//			std::cout << std::endl;

//			std::cout << "L Matrix: " << std::endl;;
//			for(int i=0;i<4;i++){
//				for(int j=0;j<4;j++){
//					std::cout << MatL[i][j] << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
		}

		// transpose L to get U
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatU[i][j] = MatL[j][i];
			}
		}
	}

	// Doolittle algorithm LU decomposition
	template<class T, unsigned int M>
	void LU_DOOLITTLE(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

		// clean the matrix
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<M;j++){
				MatL[i][j] = 0;
				MatU[i][j] = 0;
			}
		}

		// decomposition
		for(unsigned int i=0;i<M;i++){
			// upper triangular
			for(unsigned int k=0;k<M;k++){
				T tmp = Mat[i][k];
				for(unsigned int j=0;j<i;j++){
					tmp -= MatL[i][j] * MatU[j][k];
				}
				MatU[i][k] = tmp;
			}
			// lower triangular
			MatL[i][i] = 1;
			for(unsigned int k=i+1;k<M;k++){
				T tmp = Mat[k][i];
				for(unsigned int j=0;j<i;j++){
					tmp -= MatL[k][j] * MatU[j][i];
				}
				MatL[k][i] = tmp / MatU[i][i];
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
			 int dim, int order>
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
						MatTmp2[k][j] = MatTmp1[k+1][j] - MatTmp1[k][j];
						MatTmp2[k][j] *= dir;
					}
				}
				MAT_EQ<T, M, N>(MatTmp2, MatTmp1);
			}
			MAT_EQ<T, M, N>(MatTmp1, MatB);
		}else if(dim==2){
			for(int i=0;i<order;i++){
				for(int j=0;j<M;j++){
					for(int k=0;k<N-i-1;k++){
						MatTmp2[j][k] = MatTmp1[j][k+1] - MatTmp1[j][k];
						MatTmp2[j][k] *= dir;
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
	template<class T, int M, int N>
	void MAT_SIGN(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = MatA[i][j]>=0 ? 1 : -1;
	}

	// Matrix absolute value
	template<class T, int M, int N>
	void MAT_ABS(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = MatA[i][j]>=0 ? MatA[i][j] : -MatA[i][j];
	}

	// Matrix dot multiplication
	template<class T, int M, int N>
	void MAT_DOTSQUARE(T MatA[M][N],
				 T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = MatA[i][j] * MatA[i][j];
	}

	// Matrix dot multiplication
	template<class T, int M, int N>
	void MAT_DOTMUL(T MatA[M][N],
				    T MatB[M][N],
					T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] * MatB[i][j];
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
	template<class T, int M, int N>
	void MAT_COMPLEX_DOTDIV_REAL(Complex<T> MatA[M][N],
							 	 T MatB[M][N],
								 Complex<T> MatC[M][N]){
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				MatC[i][j].real = MatA[i][j].real / MatB[i][j];
				MatC[i][j].imag = MatA[i][j].imag / MatB[i][j];
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
	template<class T, int M, int N>
	void MAT_DOTINV(T MatA[M][N],
				    T MatB[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatB[i][j] = (T) 1 / MatA[i][j];
	}

	// Matrix dot multiply scalar
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
	template<class T, int M, int N>
	void MAT_SCALAR_DOTMUL(T MatA[M][N],
						   T B,
						   T MatC[M][N]){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] * B;
	}
	template<class T, int M, int N>
	void MAT_SCALAR_DOTMUL(T **MatA,
						   T B,
						   T **MatC){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i][j] = MatA[i][j] * B;
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
	template<class T, int M>
	void VEC_DIFF(T VecA,
				  T VecB,
				  int order){
		T VecTmp = VecA;
		for(int i=0;i<order;i++){
			for(int j=0;j<M-i-1;j++)
				VecB[i] = VecTmp[i+1] - VecTmp[i];
			VecTmp = VecB;
		}
	}

	// Vector absolute value
	template<class T, int M>
	void VEC_ABS(T VecA[M],
				 T VecB[M]){
		for(int i=0;i<M;i++)
			VecB[i] = VecA[i]>=0 ? VecA[i] : -VecA[i];
	}

	// Vector dot division
	template<class T, int M>
	void VEC_DOTDIV(T MatA[M],
				    T MatB[M],
					T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = MatA[i]/MatB[i];
	}

	// Vector dot division scalar
	template<class EigenT, class T, int M>
	void VEC_SCALAR_DOTDIV(T MatA[M],
						   T B,
						   T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = MatA[i]/B;
	}
	template<class EigenT, class T, int M>
	void VEC_SCALAR_DOTDIV(T A,
						   T MatB[M],
						   T MatC[M]){
		for(int i=0;i<M;i++)
			MatC[i] = A/MatB[i];
	}


private:

protected:

};



#endif /* SRC_FPT_ALGEBRA_HPP_ */
