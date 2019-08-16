/*
 *	This is a simple example of how to do
 *	proximal gradient decent algorithm in
 *	C++ with Eigen library and plot with
 *	Matplotlib C++ API
 *

 	% Solve Quadratic Program
	%
	%   minimize    0.5*x'*A*x + b'*x
	%      x
	%   subject to  -a <= x <= a
	%
	% with proximal (projected) gradient method

 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#include "pgd.hpp"
#include "pgd_eigend.hpp"
#include "pgd_eigenf.hpp"
#include "pgd_gdouble.hpp"
#include "pgd_gfloat.hpp"
#include "pgd_xfpt.hpp"
#include "pgd_xfxpt.hpp"

void PROXIMAL_GRADIENT_DECENT(){

	//Eigen::initParallel();

	//int thnum = 4;
	//omp_set_num_threads(thnum);
	//Eigen::setNbThreads(thnum);

	std::string clockname = "timeprofile.txt";
	std::string Amatrixname = "Amatrix.csv";
	std::string bvectorname = "bvector.csv";
	std::string Lname = "L.csv";

	////////////////////// generate data //////////////////////
#ifdef ALWAYS_DELETE
	std::remove(clockname.c_str());
	std::remove(Amatrixname.c_str());
	std::remove(bvectorname.c_str());
	std::remove(Lname.c_str());
	std::remove(errorhistname.c_str());
	std::remove(figurename.c_str());
#endif
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

	// check if the data file already exist
	Eigen_Algebra Eigen_Algebra_obj;
	Eigen::MatrixXd Amatrix( DIAG, DIAG );
	Eigen::MatrixXd AmatrixT( DIAG, DIAG );
	if(!checkfileexist(Amatrixname)){
		// randn(n, n)
		Eigen::MatrixXd RNDmatrix( ROW, COL );
		Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, ROW, COL>( RNDmatrix );
#ifdef DEBUG_DATA
	std::cout << RNDmatrix << std::endl << std::endl;
#endif// endif DEBUG_DATA

		// [Q, R] = qr(randn(n, n));
		Eigen::MatrixXd Qmatrix( ROW, ROW );
		Eigen::MatrixXd QmatrixT( ROW, ROW );
		Eigen::MatrixXd Rmatrix( ROW, COL );
		Eigen_Algebra_obj.QRD<Eigen::MatrixXd>(
				RNDmatrix, Qmatrix, Rmatrix, ROW, COL );
#ifdef DEBUG_DATA
	std::cout << Qmatrix << std::endl << std::endl;
	std::cout << Rmatrix << std::endl << std::endl;
#endif// endif DEBUG_DATA

		// The sampling/measurement matrix
		// A = Q*diag([50*rand(floor(0.5*n), 1); zeros(n - floor(0.5*n), 1)])*Q';
		Eigen::MatrixXd diagmatrix_tmp( DIAG, DIAG );
		Eigen::MatrixXd diagmatrix( DIAG, DIAG );
		Eigen_Algebra_obj.RND_DIAGMAT
				<Eigen::MatrixXd, Eigen::VectorXd, DIAG>(
				diagmatrix_tmp, DIAG_RATIO );
		Eigen_Algebra_obj.MAT_SCALAR_MUL<Eigen::MatrixXd, double>(
				diagmatrix_tmp, DIAG_VALUE, diagmatrix);
#ifdef DEBUG_DATA
	std::cout << diagmatrix << std::endl << std::endl;
#endif// endif DEBUG_DATA
		Eigen::MatrixXd Amatrix_tmp( DIAG, DIAG );
		Eigen_Algebra_obj.MAT_TRANS<Eigen::MatrixXd>(Qmatrix, QmatrixT);
		Eigen_Algebra_obj.MAT_MUL<Eigen::MatrixXd>(Qmatrix, diagmatrix, Amatrix_tmp);
		Eigen_Algebra_obj.MAT_MUL<Eigen::MatrixXd>(Amatrix_tmp, QmatrixT, Amatrix);
		// write to data file
		writematrix_double(Amatrix, Amatrixname, 64);
		std::cout << "created the Amatrix" << std::endl << std::endl;
	}else{
		// read from data file
		Amatrix = readmatrix_double(Amatrixname);
		std::cout << "read the Amatrix" << std::endl << std::endl;
	}
	AmatrixT = Amatrix.transpose();

	// measurement
	Eigen::VectorXd bvector( DIAG );
	if(!checkfileexist(bvectorname)){
		Eigen_Algebra_obj.RND_VEC<Eigen::VectorXd, DIAG>( bvector );
		// write to csv file
		writevector_double(bvector, bvectorname, 64);
		std::cout << "created the bvector" << std::endl << std::endl;
	}else{
		// read from data file
		bvector = readvector_double(bvectorname);
		std::cout << "read the bvector" << std::endl << std::endl;
	}

	// gradient decent step factor;
	double L;
	if(!checkfileexist(Lname)){
		Eigen::VectorXcd eigvector_complex( DIAG );
		Eigen::VectorXd eigvector_real( DIAG );
		Eigen_Algebra_obj.MAT_EIG<Eigen::MatrixXd, Eigen::VectorXcd>( Amatrix,
				eigvector_complex );
		eigvector_real = eigvector_complex.real();
		L = eigvector_real.maxCoeff();
		//L = L * std::exp(ROW);
		// write to data file
		writescalar_double(L, Lname, 64);
		std::cout << "created the L value" << std::endl << std::endl;
	}else{
		// read from data file
		L = readscalar_double(Lname);
		std::cout << "read the L value" << std::endl << std::endl;
	}
#ifdef DEBUG_DATA
	std::cout << "Amatrix: " << std::endl
			  << Amatrix
			  << std::endl << std::endl;
	std::cout << "AmatrixT: " << std::endl
			  << AmatrixT
			  << std::endl << std::endl;
	std::cout << "bvector: " << std::endl
			  << bvector
			  << std::endl << std::endl;
	std::cout << "L: " << L << std::endl << std::endl;
	std::cout << "1/L: " << 1/L << std::endl << std::endl;
#endif// endif DEBUG_DATA

//#ifdef DEBUG_DATA
//	std::cout << "Amatrix_c: " << std::endl;
//	for(int i=0;i<DIAG;i++){
//		for(int j=0;j<DIAG;j++)
//			std::cout << Amatrix_c[i][j] << " ";
//		std::cout << "\n";
//	}
//	std::cout << std::endl;
//	std::cout << "bvector_c: " << std::endl;
//	for(int i=0;i<DIAG;i++){
//		std::cout << bvector_c[i] << " ";
//	}
//	std::cout << std::endl;
//	std::cout << "L_c: " << L_c << std::endl << std::endl;
//#if defined(XILINX_FIXED_PRECISION)
//	std::cout << "factor: " << factor << std::endl << std::endl;
//#endif
//	std::exit(0);
//#endif

#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to generate and save the data"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to generate and save the data"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE
	///////////////////////////////////////////////////////////


	///////////////// Proximal gradient descent ///////////////
{
#ifdef EIGEN_DOUBLE_PRECISION
	PROXIMAL_GRADIENT_DECENT_EIGEND(Amatrix, bvector, L);
#endif
}
{
#ifdef EIGEN_FLOAT_PRECISION
	PROXIMAL_GRADIENT_DECENT_EIGENF(Amatrix, bvector, L);
#endif
}
{
#ifdef DOUBLE_PRECISION
	double Amatrix_c1[DIAG][DIAG];
	memcpy(Amatrix_c1,AmatrixT.data(),sizeof(double)*DIAG*DIAG);
	double bvector_c1[DIAG];
	memcpy(bvector_c1,bvector.data(),sizeof(double)*DIAG);
	double L_c1 = L;
	PROXIMAL_GRADIENT_DECENT_GDOUBLE(Amatrix_c1, bvector_c1, L_c1);
#endif
}
{
#ifdef FLOAT_PRECISION
	float Amatrix_c2[DIAG][DIAG];
	Eigen::MatrixXf Amatrix_f(DIAG, DIAG);
	Amatrix_f = AmatrixT.cast<float>();
	memcpy(Amatrix_c2,Amatrix_f.data(),sizeof(float)*DIAG*DIAG);
	float bvector_c2[DIAG];
	Eigen::VectorXf bvector_f(DIAG);
	bvector_f = bvector.cast<float>();
	memcpy(bvector_c2,bvector_f.data(),sizeof(float)*DIAG);
	float L_c2 = (float) L;
	PROXIMAL_GRADIENT_DECENT_GFLOAT(Amatrix_c2, bvector_c2, L_c2);
#endif// endif FLOAT_PRECISION
}
{
#if defined(COMSTOM_FLOAT_PRECISION)
	fptx2 Amatrix_c3[DIAG][DIAG];
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c3[i][j] = (fptx2) AmatrixT(i,j);
		}
	}
	fptx2 bvector_c3[DIAG];
	for(int i=0;i<DIAG;i++){
		bvector_c3[i] = (fptx2) bvector(i);
	}
	fptx2 L_c3 = (fptx2) L;
	PROXIMAL_GRADIENT_DECENT_XFLOAT2(Amatrix_c3, bvector_c3, L_c3);
#endif
}
{
#ifdef XILINX_FIXED_PRECISION
   	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebra_obj;
	double Amatrix_d[DIAG][DIAG];
	DATA_IN_T Amatrix_c3[DIAG][DIAG];
	memcpy(Amatrix_d,AmatrixT.data(),sizeof(double)*DIAG*DIAG);
	Xilinx_Fixed_Point_Algebra_obj.MAT_EQ<double, DATA_IN_T, DIAG, DIAG>(Amatrix_d, Amatrix_c3);
	double bvector_d[DIAG];
	DATA_IN_T bvector_c3[DIAG];
	memcpy(bvector_d,bvector.data(),sizeof(double)*DIAG);
	Xilinx_Fixed_Point_Algebra_obj.VEC_EQ<double, DATA_IN_T, DIAG>(bvector_d,bvector_c3);
	DATA_IN_T L_c3 = L;
	DATA_IN_T factor = 1/L;
	PROXIMAL_GRADIENT_DECENT_XFXPT(Amatrix_c3, bvector_c3, factor);
#endif // endif XILINX_FIXED_PRECISION
	///////////////////////////////////////////////////////////
}
}
int main(int argc, char** argv){
	PROXIMAL_GRADIENT_DECENT();
	return 0;
}