/*
 * matrix.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_MATH_MATRIX_HPP_
#define HEADER_MATH_MATRIX_HPP_

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
template<class T, int M, int N>
void MAT_SUB(T **A,
								      T **B,
								      T **C){
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

// The basic matrix transform
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


template<class T, int M, int N>
void ROW_OF_MATRIX(T AM[M][N], int index,
                   T V[N]){
    for(int i=0;i<N;i++){
        V[i] = AM[index][i];
    }
}

template<class T, int M, int N>
void COL_OF_MATRIX(T AM[M][N], int index,
                   T V[M]){
    for(int i=0;i<M;i++){
        V[i] = AM[i][index];
    }
}

template<class T, int M, int N>
void ROW_INTO_MATRIX(T V[M], T A[M][N], int index){
    for(int i=0;i<N;i++){
        A[index][i] = V[i];
    }
}

template<class T, int M, int N>
void COL_INTO_MATRIX(T V[M], T A[M][N], int index){
    for(int i=0;i<M;i++){
        A[i][index] = V[i];
    }
}

template<class T, int M, int N>
void SCALAR_INTO_MATRIX(T S, T A[M][N], int col, int row){
    A[col][row] = S;
}

template<class T, int M, int N>
void EYE_MATRIX(T eye[M][N], int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j)
                eye[i][j] = (T) 1;
            else
                eye[i][j] = (T) 0;
        }
    }
}

#endif /* HEADER_MATH_MATRIX_HPP_ */
