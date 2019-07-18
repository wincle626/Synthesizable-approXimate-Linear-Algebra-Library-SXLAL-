/*
 *	This is light arithmetic wrapper of Eigen
 *	library
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_EIGEN_ALGEBRA_HPP_
#define SRC_EIGEN_ALGEBRA_HPP_

#include "common.hpp"
#include "data.hpp"

class Eigen_Algebra{

public:

	// Generate all zero matrix
	template<class EigenT, int M, int N>
	void ZEROS_MAT( EigenT &Mat ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				Mat( i, j ) = 0;
			}
		}
	}
	template<class EigenT>
	void ZEROS_MAT( EigenT &Mat, int M, int N ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				Mat( i, j ) = 0;
			}
		}
	}

	// Generate all one matrix
	template<class EigenT, int M, int N>
	void ONES_MAT( EigenT &Mat ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				Mat( i, j ) = 1;
			}
		}
	}
	template<class EigenT>
	void ONES_MAT( EigenT &Mat, int M, int N ){
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				Mat( i, j ) = 1;
			}
		}
	}

	// Generate random matrix
	template<class EigenT, int M, int N>
	void RND_MAT( EigenT &Mat ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				int rnd = rand() % FLOAT_SIZE;
				double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
				Mat( i, j ) = rndnum;
			}
		}
	}
	template<class EigenT>
	void RND_MAT( EigenT &Mat, int M, int N ){
		srand (time(NULL));
		for( int i=0; i<M; i++ ){
			for( int j=0; j<N; j++ ){
				int rnd = rand() % FLOAT_SIZE;
				double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
				Mat( i, j ) = rndnum;
			}
		}
	}

	// Generate random diagonal matrix
	template<class EigenT1, class EigenT2, int M>
	void RND_DIAGMAT( EigenT1 &Mat){
		srand (time(NULL));
		int sparse_num = floor( DIAG_RATIO * M );
		EigenT2 diag_vec( M );
		for( int i=0; i<M; i++){
			if( i<sparse_num){
				int rnd = rand() % FLOAT_SIZE;
				double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
				diag_vec(i) = rndnum;
			}else{
				diag_vec(i) = 0;
			}
		}
		Mat = diag_vec.matrix().asDiagonal();
	}
	template<class EigenT1, class EigenT2>
	void RND_DIAGMAT( EigenT1 &Mat, int M){
		srand (time(NULL));
		int sparse_num = floor( DIAG_RATIO * M );
		EigenT2 diag_vec( M );
		for( int i=0; i<M; i++){
			if( i<sparse_num){
				int rnd = rand() % FLOAT_SIZE;
				double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
				diag_vec(i) = rndnum;
			}else{
				diag_vec(i) = 0;
			}
		}
		Mat = diag_vec.matrix().asDiagonal();
	}

	// Generate all zero vector
	template<class EigenT, int M>
	void ZEROS_VEC( EigenT &Vec){
		for( int i=0; i<M; i++ ){
			Vec(i) = 0;
		}
	}
	template<class EigenT>
	void ZEROS_VEC( EigenT &Vec, int M){
		for( int i=0; i<M; i++ ){
			Vec(i) = 0;
		}
	}

	// Generate all one vector
	template<class EigenT, int M>
	void ONES_VEC( EigenT &Vec){
		for( int i=0; i<M; i++ ){
			Vec(i) = 1;
		}
	}
	template<class EigenT>
	void ONES_VEC( EigenT &Vec, int M){
		for( int i=0; i<M; i++ ){
			Vec(i) = 1;
		}
	}

	// Generate random vector
	template<class EigenT, int M>
	void RND_VEC( EigenT &Vec){
		srand (time(NULL));
		for( int i=0; i<DIAG; i++ ){
			int rnd = rand() % FLOAT_SIZE;
			double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
			Vec(i) =  rndnum;
		}
	}
	template<class EigenT>
	void RND_VEC( EigenT &Vec, int M){
		srand (time(NULL));
		for( int i=0; i<DIAG; i++ ){
			int rnd = rand() % FLOAT_SIZE;
			double rndnum = INTEGER_SCALE * (double)( rnd ) / FLOAT_SIZE;
			Vec(i) =  rndnum;
		}
	}

	// Vector value given to a diagnal matrix
	template<class EigenT>
	void VEC2DIAG( EigenT Vec,
			EigenT &Mat ){
		Mat = Vec.matrix().asDiagonal();
	}

	// QR decomposition
	template<class EigenT>
	void QRD( EigenT &Mat,
			  EigenT &Q,
			  EigenT &R){
		Eigen::HouseholderQR<EigenT> qr(Mat);
		Q = qr.householderQ(); // get Q matrix
		if( ROW == COL ){
			R = qr.matrixQR().template  triangularView<Eigen::Upper>();
		}else{
			Eigen::FullPivLU<EigenT>lu_decomp(Mat);
			int Rank = lu_decomp.rank(); //retrieve rank of matrix
			R = qr.matrixQR().topLeftCorner(Rank, Rank).template
				triangularView<Eigen::Upper>(); // get R matrix
		}
	}

	// Matrix transpose
	template<class EigenT>
	void MAT_TRANS(EigenT &Mat,
			       EigenT &MatT){
		MatT = Mat.transpose();
	}

	// Matrix eigen value
	template<class EigenT1, class EigenT2>
	void MAT_EIG(EigenT1 &Mat,
				 EigenT2 &Eigvec){
		Eigvec = Mat.eigenvalues();
	}

	// Matrix multiplication
	template<class EigenT>
	void MAT_MUL( EigenT MatA,
				  EigenT MatB,
				  EigenT &MatC ){
		MatC = MatA * MatB;
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

	// Extend genera size matrix to square matrix
	template<class EigenT, int M, int N>
	void EXTEND_GENERAL_TO_SQUARE(EigenT A,
								  EigenT &B){
		B.block(0,0,M,N) = A;
	}
	template<class EigenT>
	void EXTEND_GENERAL_TO_SQUARE(EigenT A,
								  EigenT &B,
								  int M, int N){
		B.block(0,0,M,N) = A;
	}

	// Matrix multiplication with Strassen method
	template<class EigenT>
	void SQUARE_MAT_MUL_STRASSEN(int m,
			EigenT A,
			EigenT B,
			EigenT &C){

		EigenT A11(m/2,m/2), A12(m/2,m/2),
				A21(m/2,m/2), A22(m/2,m/2);
		EigenT B11(m/2,m/2), B12(m/2,m/2),
				B21(m/2,m/2), B22(m/2,m/2);
		EigenT C11(m/2,m/2), C12(m/2,m/2),
				C21(m/2,m/2), C22(m/2,m/2);
		EigenT M1(m/2,m/2), M2(m/2,m/2),
				M3(m/2,m/2), M4(m/2,m/2), M5(m/2,m/2),
				M6(m/2,m/2), M7(m/2,m/2);
		EigenT AA(m/2,m/2), BB(m/2,m/2);

		if(m == 2) {  //2-order
			C = A * B;
		} else {
			// Divide and conquer
			A11 = A.block(0,0, m/2,m/2);
			A12 = A.block(0,m/2, m/2,m/2);
			A21 = A.block(m/2,0, m/2,m/2);
			A22 = A.block(m/2,m/2, m/2,m/2);
			B11 = B.block(0,0, m/2,m/2);
			B12 = B.block(0,m/2, m/2,m/2);
			B21 = B.block(m/2,0, m/2,m/2);
			B22 = B.block(m/2,m/2, m/2,m/2);
			/*std::cout << A11 << std::endl
					  << A12 << std::endl
					  << A21 << std::endl
					  << A22 << std::endl;
			std::cout << B11 << std::endl
					  << B12 << std::endl
					  << B21 << std::endl
					  << B22 << std::endl;
			std::exit(0);*/

			//Calculate M1 = (A11 + A22) × (B11 + B22)
			AA = A11 + A22;
			BB = B11 + B22;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, AA, BB, M1);

			//Calculate M2 = (A21 + A22) × B11
			AA = A21 + A22;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, AA, B11, M2);

			//Calculate M3 = A11 × (B12 - B22)
			BB = B12 - B22;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, A11, BB, M3);

			//Calculate M4 = A22 × (B21 - B11)
			BB = B21 - B11;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, A22, BB, M4);

			//Calculate M5 = (A11 + A12) × B22
			AA = A11 + A12;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, AA, B22, M5);

			//Calculate M6 = (A21 - A11) × (B11 + B12)
			AA = A21 - A11;
			BB = B11 + B12;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, AA, BB, M6);

			//Calculate M7 = (A12 - A22) × (B21 + B22)
			AA = A12 - A22;
			BB = B21 + B22;
			SQUARE_MAT_MUL_STRASSEN<EigenT>(m/2, AA, BB, M7);

			//Calculate C11 = M1 + M4 - M5 + M7
			AA = M1 + M4;
			BB = M5 - M7;
			C11 = AA - BB;

			//Calculate C12 = M3 + M5
			C12 = M3 + M5;

			//Calculate C21 = M2 + M4
			C21 = M2 + M4;

			//Calculate C22 = M1 - M2 + M3 + M6
			AA = M1 - M2;
			BB = M3 + M6;
			C22 = AA + BB;

			//Set the result to C
			C.block(0,0,m/2,m/2) = C11;
			C.block(0,m/2,m/2,m/2) = C12;
			C.block(m/2,0,m/2,m/2) = C21;
			C.block(m/2,m/2,m/2,m/2) = C22;
		}
	}
	/*void GENERAL_MAT_MUL_DOUBLE_STRASSEN(
			int m, int n, int p,
			Eigen::MatrixXd A,
			Eigen::MatrixXd B,
			Eigen::MatrixXd &C){

		if(m<=2 || n<=2 || p<=2){
			std::cout << "normal case" << std::endl;
			MAT_MUL_DOUBLE(A, B, C);
		}else if( (m < n && n < p) || (m < n && n > p && p > m)){
			std::cout << "first case" << std::endl;
			Eigen::MatrixXd A1(m, m), A2(m, n-m);
			Eigen::MatrixXd B11(m, m), B12(m, p-m),
					B21(n-m, m), B22(n-m, p-m);
			Eigen::MatrixXd CC1(m, m), CC2(m, m), CC3(m, p-m), CC4(m, p-m);
			Eigen::MatrixXd C1(m, m), C2(m, p-m);
			A1 = A.block(0, 0, m, m);
			A2 = A.block(0, m, m, n-m);
			B11 = B.block(0, 0, m, m);
			B12 = B.block(0, m, m, p-m);
			B21 = B.block(m, 0, n-m, m);
			B22 = B.block(m, m, n-m, p-m);
			int M = FIND_NEAREST_2SQUARE(m, m, m);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, A1, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, B11, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			CC1 = C_ext.block(0, 0, m, m);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, n-m, m, A2, B21, CC2);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, m, p-m, A1, B12, CC3);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, n-m, p-m, A2, B22, CC4);
			MAT_ADD_DOUBLE(CC1, CC2, C1);
			MAT_ADD_DOUBLE(CC3, CC4, C2);
			C.block(0, 0, m, m) = C1;
			C.block(0, m, m, p-m) = C2;
		}else if( m < n && n > p && p == m ) {
			Eigen::MatrixXd A1(m, m), A2(m, n-m);
			Eigen::MatrixXd B1(m, m), B2(n-m, m);
			Eigen::MatrixXd C1(m, m), C2(m, m);
			A1 = A.block(0, 0, m, m);
			A2 = A.block(0, m, m, n-m);
			B1 = B.block(0, 0, m, m);
			B2 = B.block(m, 0, n-m, m);
			int M = FIND_NEAREST_2SQUARE(m, m, m);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, A1, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, B1, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			C1 = C_ext.block(0, 0, m, m);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, n-m, m, A2, B2, C2);
			C = C1 + C2;
		}else if( m > n && n == p  ) {
			Eigen::MatrixXd A1(n, n), A2(m-n, n);
			Eigen::MatrixXd C1(n, n), C2(m-n, n);
			A1 = A.block(0, 0, n, n);
			A2 = A.block(n, 0, m-n, n);
			int M = FIND_NEAREST_2SQUARE(n, n, n);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(n, n, A1, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(n, n, B, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			C1 = C_ext.block(0, 0, n, n);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m-n, n, n, A2, B, C2);
			C.block(0, 0, n, n) = C1;
			C.block(n, 0, m-n, n) = C2;
		}else if( m < n && n == p  ) { // wrong
			Eigen::MatrixXd A1(m, m), A2(m, n-m);
			Eigen::MatrixXd B11(m, m), B12(m, n-m), B21(n-m, m), B22(n-m, n-m);
			Eigen::MatrixXd C11(m, m), C12(m, n-m), C21(m, n-m), C22(n-m,, n-m);
			A1 = A.block(0, 0, m, m);
			A2 = A.block(0, m, m, n-m);
			int M = FIND_NEAREST_2SQUARE(m, m, m);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, A1, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, B, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			C1 = C_ext.block(0, 0, m, m);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, n-m, m, A2, B, C2);
			C.block(0, 0, m, m) = C1;
			C.block(0, m, n-m, m) = C2;
		}else if( m == n && n < p ){
			Eigen::MatrixXd B1(n, n), B2(m-n, n);
			Eigen::MatrixXd C1(n, n), C2(m-n, n);
			B1 = B.block(0, 0, m, m);
			B2 = B.block(0, m, m, p-m);
			int M = FIND_NEAREST_2SQUARE(m, m, m);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, A, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(m, m, B1, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			C1 = C_ext.block(0, 0, m, m);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m, m, p-m, A, B2, C2);
			C.block(0, 0, m, m) = C1;
			C.block(0, m, m, p-m) = C2;
		}else if( m > n && n < p ){
			std::cout << "second case" << std::endl;
			Eigen::MatrixXd A1(n, n), A2(m-n, n);
			Eigen::MatrixXd B1(n, n), B2(n, p-n);
			Eigen::MatrixXd C11(n, n), C12(n, p-n),
					C21(m-n, n), C22(m-n, p-n);
			A1 = A.block(0, 0, n, n);
			A2 = A.block(n, 0, m-n, n);
			B1 = B.block(0, 0, n, n);
			B2 = B.block(0, n, n, p-n);
			int M = FIND_NEAREST_2SQUARE(n, n, n);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(n, n, A1, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(n, n, B1, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			C11 = C_ext.block(0, 0, n, n);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(n, n, p-n, A1, B2, C12);
			std::cout << __FILE__ << "--" << __LINE__ << std::endl;
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m-n, n, n, A2, B1, C21);
			std::cout << __FILE__ << "--" << __LINE__ << std::endl;
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m-n, n, p-n, A2, B2, C22);
			std::cout << __FILE__ << "--" << __LINE__ << std::endl;
			C.block(0, 0, n, n) = C11;
			C.block(0, n, n, p-n) = C12;
			C.block(n, 0, m-n, n) = C21;
			C.block(n, n, m-n, p-n) = C22;
		}else if( (m > n && n > p) || (m < n && n > p && p < m)
				|| (m == n && n > p ) ){
			std::cout << "third case" << std::endl;
			Eigen::MatrixXd A11(p, p), A12(p, n-p),
					A21(m-p, p), A22(m-p, n-p);
			Eigen::MatrixXd B1(p, p), B2(n-p, p);
			Eigen::MatrixXd CC1(p, p), CC2(p, p),
					CC3(m-p, p), CC4(m-p, p);
			Eigen::MatrixXd C1(p, p), C2(m-p, p);
			A11 = A.block(0, 0, p, p);
			A12 = A.block(0, p, p, n-p);
			A21 = A.block(p, 0, m-p, p);
			A22 = A.block(p, p, m-p, n-p);
			B1 = B.block(0, 0, p, p);
			B2 = B.block(p, 0, n-p, p);
			int M = FIND_NEAREST_2SQUARE(p, p, p);
			Eigen::MatrixXd A_ext( M, M );
			Eigen::MatrixXd B_ext( M, M );
			Eigen::MatrixXd C_ext( M, M );
			ZEROS_MAT_DOUBLE(M, M, A_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(p, p, A11, A_ext);
			ZEROS_MAT_DOUBLE(M, M, B_ext);
			EXTEND_GENERAL_TO_SQUARE_DOUBLE(p, p, B1, B_ext);
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(M, A_ext, B_ext, C_ext);
			CC1 = C_ext.block(0, 0, p, p);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(p, n-p, p, A12, B2, CC2);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m-p, p, p, A21, B1, CC3);
			GENERAL_MAT_MUL_DOUBLE_STRASSEN(m-p, n-p, p, A22, B2, CC4);
			MAT_ADD_DOUBLE(CC1, CC2, C1);
			MAT_ADD_DOUBLE(CC3, CC4, C2);
			C.block(0, 0, p, p) = C1;
			C.block(p, 0, m-p, p) = C2;
		}else{
			SQUARE_MAT_MUL_DOUBLE_STRASSEN(m, A, B, C);
		}
	}
	*/

	// Matrix dot multiplication
	template<class EigenT, int M, int N>
	void MAT_DOT_MUL( EigenT MatA,
					  EigenT MatB,
					  EigenT &MatC ){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				MatC(i,j) = MatA(i,j) * MatB(i,j);
			}
		}
	}
	template<class EigenT>
	void MAT_DOT_MUL( EigenT MatA,
					  EigenT MatB,
					  EigenT &MatC,
					  int M, int N ){
		for ( int i=0; i<M; i++ ){
			for ( int j=0; j<N; j++ ){
				MatC(i,j) = MatA(i,j) * MatB(i,j);
			}
		}
	}
	template<class EigenT>
	void MAT_DOT_MUL( EigenT MatA,
					  EigenT MatB,
					  EigenT &MatC ){
		MatC = MatA.dot(MatB);
	}

	// Matrix addition
	template<class EigenT>
	void MAT_ADD( EigenT MatA,
				  EigenT MatB,
				  EigenT &MatC ){
		MatC = MatA + MatB;
	}

	// Matrix subtraction
	template<class EigenT>
	void MAT_SUB( EigenT MatA,
				  EigenT MatB,
				  EigenT &MatC ){
		MatC = MatA - MatB;
	}

	// Matrix Vector multiplication
	template<class EigenT1, class EigenT2>
	void MAT_VEC_MUL( EigenT1 MatA,
					  EigenT2 VecB,
					  EigenT2 &VecC){
		VecC = MatA * VecB;
	}

	// Vector matrix multiplication
	template<class EigenT1, class EigenT2>
	void VEC_MAT_MUL( EigenT1 VecA,
					  EigenT2 MatB,
					  EigenT1 &VecC){
		VecC = VecA * MatB;
	}

	// Matrix scalar addition
	template<class EigenT, class Type>
	void MAT_SCALAR_ADD( EigenT MatA,
					   Type scalar,
					   EigenT &MatC){
		MatC = MatA + scalar;
	}

	// Matrix scalar subtraction
	template<class EigenT, class Type>
	void MAT_SCALAR_SUB( EigenT MatA,
					   Type scalar,
					   EigenT &MatC){
		MatC = MatA - scalar;
	}

	// Matrix scalar multiplication
	template<class EigenT, class Type>
	void MAT_SCALAR_MUL( EigenT MatA,
					   Type scalar,
					   EigenT &MatC){
		MatC = MatA * scalar;
	}

	// Extract submatrix
	// input matrix must larger than output matrix
	// output matrix (N-M)*(L-K) begins at row M, col K
	template<class EigenT, int M, int N, int K, int L>
	void MAT_SUBMAT(EigenT MatA, EigenT &MatB){
		MatB.resize(N-M,L-K);
		MatB = MatA.block<N-M,L-K>(M,K);
	}
	template<class EigenT, int M, int N, int K, int L,
			int step_r, int step_c>
	void MAT_SUBMAT(EigenT MatA, EigenT &MatB){
		for(int i=M-1;i<N;i+=step_r){
			for(int j=K-1;j<L;j+=step_c){
				MatB(i,j) = MatA(i,j);
			}
		}
	}

	// Matrix merge
	// MatA: (M,N), MatB (M,K), MatC: (M, N+K)
	template<class EigenT, int M, int N, int K>
	void MAT_MERGE_H(EigenT MatA,
					 EigenT MatB,
					 EigenT &MatC){
		MatC.resize(M,N+K);
		MatC.block<M,N>(0,0) = MatA;
		MatC.block<M,K>(0,N) = MatB;
	}
	// MatA: (M,K), MatB: (N,K), MatC: (M+N,K)
	template<class EigenT, int M, int N, int K>
	void MAT_MERGE_V(EigenT MatA,
					 EigenT MatB,
					 EigenT &MatC){
		MatC.resize(M+N,K);
		MatC.block<M,K>(0,0) = MatA;
		MatC.block<N,K>(M,0) = MatB;
	}

	// Vector addition
	template<class EigenT>
	void VEC_ADD( EigenT VecA,
				  EigenT VecB,
				  EigenT &VecC){
		VecC = VecA + VecB;
	}

	// Vector subtraction
	template<class EigenT>
	void VEC_SUB( EigenT VecA,
				  EigenT VecB,
				  EigenT &VecC){
		VecC = VecA - VecB;
	}

	// Vector dot multiplication
	template<class EigenT, int M>
	void VEC_DOT_MUL( EigenT VecA,
					  EigenT VecB,
					  EigenT &VecC ){
		for ( int i=0; i<M; i++ ){
			VecC(i) = VecA(i) * VecB(i);
		}
	}

	// Vector multiplication
	template<class EigenT1, class EigenT2, class EigenT3>
	void VEC_MUL( EigenT1 VecA,
				  EigenT2 VecB,
				  EigenT3 &C ){
		C = VecA * VecB;
	}

	// Vector scalar addition
	template<class EigenT1, class EigenT2>
	void VEC_SCALAR_ADD( EigenT1 VecA,
						 EigenT2 scalar,
						 EigenT1 &VecC){
		VecC = VecA + scalar;
	}

	// Vector scalar subtraction
	template<class EigenT1, class EigenT2>
	void VEC_SCALAR_SUB( EigenT1 VecA,
						 EigenT2 scalar,
						 EigenT1 &VecC){
		VecC = VecA - scalar;
	}

	// Vector scalar multiplication
	template<class EigenT1, class EigenT2>
	void VEC_SCALAR_MUL( EigenT1 VecA,
						 EigenT2 scalar,
						 EigenT1 &VecC){
		VecC = VecA * scalar;
	}

	// Vector scalar comparison
	template<class EigenT1, class EigenT2, class EigenTtmp, int M>
	void VEC_SCALAR_MIN( EigenT1 V1, EigenT2 S, EigenT1 &V3 ){
		EigenT1 a_ones_vec( M );
		this->ONES_VEC<EigenT1, M>(a_ones_vec);
		a_ones_vec = a_ones_vec * S;
		EigenTtmp tmp_min = V1.array().min(a_ones_vec.array());
		V3 = tmp_min.matrix();
	}
	template<class EigenT1, class EigenT2, class EigenTtmp, int M>
	void VEC_SCALAR_MAX( EigenT1 V1, EigenT2 S, EigenT1 &V3 ){
		EigenT1 a_ones_vec( M );
		this->ONES_VEC<EigenT1, M>(a_ones_vec);
		a_ones_vec = a_ones_vec * S;
		EigenTtmp tmp_max = V1.array().max(a_ones_vec.array());
		V3 = tmp_max.matrix();
	}

	// Vector Norm
	template<class EigenT, class T>
	void VEC_NORM(EigenT V, T &norm){
		norm = V.norm();
	}

	// Extract subvector
	template<class EigenT, int M, int N>
	void VEC_SUBVEC(EigenT VecA, EigenT &VecB){
		for(int i=M-1;i<N;i++)
			VecA(i) = VecB(i);
	}
	template<class EigenT, int M, int N, int step>
	void VEC_SUBVEC(EigenT VecA, EigenT &VecB){
		for(int i=M-1;i<N;i+=step)
			VecA(i) = VecB(i);
	}

	// Vector merge
	template<class EigenT, int M, int N>
	void VEC_MERGE2VEC(EigenT VecA,
				   EigenT VecB,
				   EigenT &VecC){
		VecC.resize(VecA.size() + VecB.size());
		VecC << VecA, VecB;

	}
	template<class EigenT1, class EigenT2, int M>
	void VEC_MERGE2MAT(EigenT1 VecA,
			   	   	   EigenT1 VecB,
					   EigenT2 &MatC){
		MatC.resize(2,VecA.size());
		MatC << VecA, VecB;
	}


private:

protected:

};

#endif  /* SRC_EIGEN_ALGEBRA_HPP_ */
