/*
 *	This is a example of approximate
 *	matrix multiplication
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#include "mm.hpp"

void MATRIX_MULTIPLY_TEST_EIGEN(){

	Eigen_Algebra Eigen_Algebra_obj;

	Eigen::MatrixXd A( ROW1, COL1 );
	Eigen::MatrixXd B( ROW2, COL2 );
	int M = Eigen_Algebra_obj.FIND_NEAREST_2SQUARE(A.rows(), A.cols(), B.cols());
	Eigen::MatrixXd A_ext( M, M );
	Eigen::MatrixXd B_ext( M, M );

	Eigen::MatrixXd C1( ROW1, COL2 ), C2( M, M ), C3( ROW1, COL2 );
	Eigen::MatrixXd D1( ROW1, COL2 ), D2( ROW1, COL2 );

	Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, ROW1, COL1>( A );
	usleep(1000);
	Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, ROW2, COL2>( B );
	/*std::cout << A << std::endl << std::endl
			<< B << std::endl << std::endl;*/

	Eigen_Algebra_obj.MAT_MUL<Eigen::MatrixXd>(A, B, C1);
	//std::cout << C1 << std::endl << std::endl;

	Eigen_Algebra_obj.ZEROS_MAT<Eigen::MatrixXd>( A_ext, M, M );
	Eigen_Algebra_obj.EXTEND_GENERAL_TO_SQUARE<Eigen::MatrixXd>(
			A, A_ext, A.rows(), A.cols() );
	Eigen_Algebra_obj.ZEROS_MAT<Eigen::MatrixXd>( B_ext, M, M );
	Eigen_Algebra_obj.EXTEND_GENERAL_TO_SQUARE<Eigen::MatrixXd>(
			B, B_ext, B.rows(), B.cols() );
	/*std::cout << A_ext << std::endl << std::endl
			<< B_ext << std::endl << std::endl;*/

	Eigen_Algebra_obj.SQUARE_MAT_MUL_STRASSEN<Eigen::MatrixXd>(M, A_ext, B_ext, C2);
	//std::cout << C2 << std::endl << std::endl;

	/*Eigen_Algebra::GENERAL_MAT_MUL_DOUBLE_STRASSEN(A.rows(), A.cols(), B.cols(),
			A, B, C3);
	std::cout <<  C3 << std::endl << std::endl;*/

	D1 = C1 - C2.block(0, 0, ROW1, COL2 );
	//D2 = C1 - C3;
	std::cout << D1.sum() << std::endl;
	//std::cout << D2.sum() << std::endl;

}

void MATRIX_MULTIPLY_TEST_GENERAL(){

}

void MATRIX_MULTIPLY_TEST_XFXPT(){

}

int main(int argc, char** argv){

	MATRIX_MULTIPLY_TEST_EIGEN();
	MATRIX_MULTIPLY_TEST_GENERAL();
	MATRIX_MULTIPLY_TEST_XFXPT();
	return 0;

}
