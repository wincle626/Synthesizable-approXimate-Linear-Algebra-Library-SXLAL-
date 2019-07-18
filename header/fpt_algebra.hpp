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

	// Generate random matrix
	template<class T, int M, int N>
	void RND_MAT( T A[M][N] ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
				double rndnum = (double)( rnd ) / FLOAT_SIZE;
				A[i][j] = rndnum;
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
					int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
					double rndnum = (double)( rnd ) / FLOAT_SIZE;
					A[i][j] = rndnum;
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
			int rnd = INTEGER_SCALE * rand() % FLOAT_SIZE;
			double rndnum = (double)( rnd ) / FLOAT_SIZE;
			V[i] = rndnum;
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
		S = (T2) std::sqrt(S);
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
	void GENERAL_MAT_ADD_BASIC(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] + B[i][j];
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void GENERAL_MAT_ADD_BASIC(T1 A[M][N],
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
	void GENERAL_MAT_SUB_BASIC(T A[M][N],
									      T B[M][N],
									      T C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] - B[i][j];
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void GENERAL_MAT_SUB_BASIC(T1 A[M][N],
									      T2 B[M][N],
									      T3 C[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				C[i][j] = A[i][j] - B[i][j];
			}
		}
	}

	// Matrix Vector multiplication
	template<class T, int M, int N>
	void GENERAL_MAT_VEC_MUL_BASIC(T A[M][N],
										  T B[N],
										  T C[M]){
		for ( int i=0; i<M; i++ ){
			C[i] = 0;
			for ( int j=0; j<N; j++ ){
				C[i] += A[i][j] * B[j];
			}
		}
	}
	template<class T1, class T2, class T3, int M, int N>
	void GENERAL_MAT_VEC_MUL_BASIC(T1 A[M][N],
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
	void GENERAL_VEC_MAT_MUL_BASIC(T A[M],
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
	void GENERAL_VEC_MAT_MUL_BASIC(T1 A[M],
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
	void GENERAL_MAT_EQ_BASIC(T A[M][N],
									      T B[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				A[i][j] = B[i][j];
			}
		}
	}
	template<class T1, class T2, int M, int N>
	void GENERAL_MAT_EQ_BASIC(T1 A[M][N],
									      T2 B[M][N]){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				B[i][j] = (T2) A[i][j];
			}
		}
	}

	// The basic arbitrary matrix multiplication with complexity of O(MNP)
	template<class T, int M, int N, int P>
	void GENERAL_MAT_MUL_BASIC(T A[M][N],
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
	template<class T1, class T2, class T3, int M, int N, int P>
	void GENERAL_MAT_MUL_BASIC(T1 A[M][N],
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
	void SQUARE_MAT_MUL_BASIC(T A[M][M],
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
	void SQUARE_MAT_MUL_BASIC(T1 A[M][M],
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

	// Extend general matrix to square, the square size must be largest
	template<class T, int M, int N, int S>
	void EXTEND_GENERAL_TO_SQUARE(T A[M][N],
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
	void EXTEND_GENERAL_TO_SQUARE(T1 A[M][N],
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
	    	SQUARE_MAT_MUL_BASIC<T, 2>( A, B, C );
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
			GENERAL_MAT_ADD_BASIC<T, m/2>( A11, A22, AA );
			GENERAL_MAT_ADD_BASIC<T, m/2>( B11, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M1 );

			//Calculate M2 = (A2 + A3) × B0
			GENERAL_MAT_ADD_BASIC<T, m/2>( A21, A22, AA );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, B11, M2 );

			//Calculate M3 = A0 × (B1 - B3)
			GENERAL_MAT_SUB_BASIC<T, m/2>( B12, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, A11, BB, M3 );

			//Calculate M4 = A3 × (B2 - B0)
			GENERAL_MAT_SUB_BASIC<T, m/2>( B21, B11, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, A22, BB, M4 );

			//Calculate M5 = (A0 + A1) × B3
			GENERAL_MAT_ADD_BASIC<T, m/2>( A11, A12, AA );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, B22, M5);

			//Calculate M6 = (A2 - A0) × (B0 + B1)
			GENERAL_MAT_SUB_BASIC<T, m/2>( A21, A11, AA );
			GENERAL_MAT_ADD_BASIC<T, m/2>( B11, B12, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M6 );

			//Calculate M7 = (A1 - A3) × (B2 + B3)
			GENERAL_MAT_SUB_BASIC<T, m/2>( A12, A22, AA );
			GENERAL_MAT_ADD_BASIC<T, m/2>( B21, B22, BB );
			SQUARE_MAT_MUL_STRASSEN<T, m/2>( m/2, AA, BB, M7 );

			//Calculate C0 = M1 + M4 - M5 + M7
			GENERAL_MAT_ADD_BASIC<T, m/2>( M1, M4, AA );
			GENERAL_MAT_ADD_BASIC<T, m/2>( M5, M7, BB );
			GENERAL_MAT_SUB_BASIC<T, m/2>( AA, BB, C11 );

			//Calculate C1 = M3 + M5
			GENERAL_MAT_ADD_BASIC<T, m/2>( M3, M5, C12 );

			//Calculate C2 = M2 + M4
			GENERAL_MAT_ADD_BASIC<T, m/2>( M2, M4, C21 );

			//Calculate C3 = M1 - M2 + M3 + M6
			GENERAL_MAT_SUB_BASIC<T, m/2>( M1, M2, AA);
			GENERAL_MAT_ADD_BASIC<T, m/2>( M3, M6, BB);
			GENERAL_MAT_ADD_BASIC<T, m/2>( AA, BB, C22);

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

	// Extract submatrix
	// input matrix must larger than output matrix
	// output matrix (N-M)*(L-K) begins at row M, col K
	template<class T, int M, int N, int K, int L>
	void MAT_SUBMAT(T MatA, T &MatB){
		for(int i=M-1;i<N;i++){
			for(int j=K-1;j<L;j++){
				MatB[i,j] = MatA[i,j];
			}
		}
	}
	template<class T, int M, int N, int K, int L,
			int step_r, int step_c>
	void MAT_SUBMAT(T MatA, T &MatB){
		for(int i=M-1;i<N;i+=step_r){
			for(int j=K-1;j<L;j+=step_c){
				MatB[i,j] = MatA[i,j];
			}
		}
	}

	// Matrix merge
	// MatA: (M,N), MatB (M,K), MatC: (M, N+K)
	template<class T, int M, int N, int K>
	void MAT_MERGE_H(T MatA,
					 T MatB,
					 T &MatC){
		for(int i=0;i<M;i++)
			for(int j=0;j<N;j++)
				MatC[i,j] = MatA[i,j];
		for(int i=0;i<M;i++)
			for(int j=N,k=0;j<N+K,k<K;j++,k++)
				MatC[i,j] = MatB[i,k];
	}
	// MatA: (M,K), MatB: (N,K), MatC: (M+N,K)
	template<class T, int M, int N, int K>
	void MAT_MERGE_V(T MatA,
					 T MatB,
					 T &MatC){
		for(int i=0;i<M;i++)
			for(int j=0;j<K;j++)
				MatC[i,j] = MatA[i,j];
		for(int i=M,k=0;i<M+K,k<N;i++,k++)
			for(int j=0;j<K;j++)
				MatC[i,j] = MatB[k,j];
	}

	// Extract subvector
	template<class T, int M, int N>
	void VEC_SUBVEC(T VecA, T &VecB){
		for(int i=M-1;i<N;i++)
			VecA[i] = VecB[i];
	}
	template<class T, int M, int N, int step>
	void VEC_SUBVEC(T VecA, T &VecB){
		for(int i=M-1;i<N;i+=step)
			VecA[i] = VecB[i];
	}

	// Vector merge
	template<class T, int M, int N>
	void VEC_MERGE2VEC(T VecA,
				   T VecB,
				   T &VecC){
		for(int i=0;i<M;i++)
			VecC[i] = VecA[i];
		for(int i=M;i<N+M;i++)
			VecC[i] = VecB[i];

	}
	template<class EigenT1, class EigenT2, int M>
	void VEC_MERGE2MAT(EigenT1 VecA,
			   	   	   EigenT1 VecB,
					   EigenT2 &MatC){
		for(int i=0;i<M;i++)
			MatC[0,i] = VecA[i];
		for(int i=0;i<M;i++)
			MatC[1,i] = VecB[i];
	}

private:

protected:

};



#endif /* SRC_FPT_ALGEBRA_HPP_ */
