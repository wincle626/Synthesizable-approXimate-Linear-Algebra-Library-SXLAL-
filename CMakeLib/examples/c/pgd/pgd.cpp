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
#include "pgd_softposit32.hpp"
#include "pgd_softposit16.hpp"
#include "pgd_softposit8.hpp"
#include "pgd_softpositX.hpp"
#include "pgd_xfxpt.hpp"
#include "pgd_integer.hpp"
#include "pgd_xfloat.hpp"

//#include "MatlabEngine.hpp"
//#include "MatlabDataArray.hpp"
//
//void callFevalrandn(Eigen::MatrixXd &Mat, int M, int N){
//
//
//    // Pass vector containing MATLAB data array scalar
//    using namespace matlab::engine;
//
//    // Start MATLAB engine synchronously
//    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
//
//    // Create MATLAB data array factory
//    matlab::data::ArrayFactory factory;
//
//    // Create MATLAB input data array factory
//    std::vector<matlab::data::Array> args({
//    	factory.createScalar<int16_t>(M),
//		factory.createScalar<int16_t>(N)});
//
//    // Call MATLAB function and return result
//    matlab::data::TypedArray<double> result = matlabPtr->feval(u"randn", args);
//    std::cout << "Result: " << std::endl;
//    for(int i=0;i<M;i++){
//    	for(int j=0;j<N;j++){
//			double v = result[i][j];
//			Mat(i,j) = v;
//    	}
//    }
//}
//
//void callFevalrandn(Eigen::VectorXd &Vec, int M){
//
//
//    // Pass vector containing MATLAB data array scalar
//    using namespace matlab::engine;
//
//    // Start MATLAB engine synchronously
//    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
//
//    // Create MATLAB data array factory
//    matlab::data::ArrayFactory factory;
//
//    // Create MATLAB input data array factory
//    std::vector<matlab::data::Array> args({
//    	factory.createScalar<int16_t>(M),
//		factory.createScalar<int16_t>(1)});
//
//    // Call MATLAB function and return result
//    matlab::data::TypedArray<double> result = matlabPtr->feval(u"randn", args);
//    std::cout << "Result: " << std::endl;
//    for(int i=0;i<M;i++){
//			double v = result[i][0];
//			Vec(i) = v;
//    }
//}

void PROXIMAL_GRADIENT_DECENT(std::string path){
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;

	//Eigen::initParallel();

	//int thnum = 4;
	//omp_set_num_threads(thnum);
	//Eigen::setNbThreads(thnum);

	std::string clockname = path+"/timeprofile.txt";
	std::string Amatrixname = path+"/Amatrix.csv";
	std::string bvectorname = path+"/bvector.csv";
	std::string Lname = path+"/L.csv";

	////////////////////// generate data //////////////////////
#ifdef ALWAYS_DELETE
	std::remove(clockname.c_str());
	std::remove(Amatrixname.c_str());
	std::remove(bvectorname.c_str());
	std::remove(Lname.c_str());
//	std::remove(errorhistname.c_str());
//	std::remove(figurename.c_str());
#endif
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;

	// check if the data file already exist
	Eigen_Algebra Eigen_Algebra_obj;
	Eigen::MatrixXd Amatrix( DIAG, DIAG );
	Eigen::MatrixXd AmatrixT( DIAG, DIAG );
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;
	if(!checkfileexist(Amatrixname)){
		// randn(n, n)
		Eigen::MatrixXd RNDmatrix( ROW, COL );
		Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, ROW, COL>( RNDmatrix );
//		callFevalrandn(RNDmatrix, ROW, COL);
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
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;
	AmatrixT = Amatrix.transpose();
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;

	// measurement
	Eigen::VectorXd bvector( DIAG );
	if(!checkfileexist(bvectorname)){
		Eigen_Algebra_obj.RND_VEC<Eigen::VectorXd, DIAG>( bvector );
//		callFevalrandn(bvector, DIAG);
		// write to csv file
		writevector_double(bvector, bvectorname, 64);
		std::cout << "created the bvector" << std::endl << std::endl;
	}else{
		// read from data file
//		std::cout << __FILE__ << "," << __LINE__ << std::endl;
		bvector = readvector_double(bvectorname);
		std::cout << "read the bvector" << std::endl << std::endl;
	}
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;

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
#ifdef GENERAL_DOUBLE_PRECISION
#if DIAG < 1024
	double Amatrix_c1[DIAG][DIAG];
	double bvector_c1[DIAG];
	double L_c1 = L;
	memcpy(Amatrix_c1,AmatrixT.data(),sizeof(double)*DIAG*DIAG);
	memcpy(bvector_c1,bvector.data(),sizeof(double)*DIAG);
	PROXIMAL_GRADIENT_DECENT_GDOUBLE(Amatrix_c1, bvector_c1, L_c1);
#else
	double **Amatrix_c1;
	double *bvector_c1;
	double L_c1 = L;
	Amatrix_c1 = (double**) malloc(sizeof(double*)*DIAG);
	for(int i=0;i<DIAG;i++)
		Amatrix_c1[i] = (double*) malloc(sizeof(double)*DIAG);
	bvector_c1 = (double*) malloc(sizeof(double)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c1[i][j] = Amatrix(i,j);
		}
		bvector_c1[i] = bvector(i);
	}
//	memcpy(Amatrix_c1,AmatrixT.data(),sizeof(double)*DIAG*DIAG);
//	memcpy(bvector_c1,bvector.data(),sizeof(double)*DIAG);
	PROXIMAL_GRADIENT_DECENT_GDOUBLE(Amatrix_c1, bvector_c1, L_c1);
#endif
#endif
}

{
#ifdef GENERAL_FLOAT_PRECISION
#if DIAG < 1408
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
#else
	float **Amatrix_c2;
	float *bvector_c2;
	float L_c2 = L;
	Amatrix_c2 = (float**) malloc(sizeof(float*)*DIAG);
	for(int i=0;i<DIAG;i++)
		Amatrix_c2[i] = (float*) malloc(sizeof(float)*DIAG);
	bvector_c2 = (float*) malloc(sizeof(float)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c2[i][j] = (float)Amatrix(i,j);
		}
		bvector_c2[i] = (float)bvector(i);
	}
	PROXIMAL_GRADIENT_DECENT_GFLOAT(Amatrix_c2, bvector_c2, L_c2);
#endif
#endif// endif FLOAT_PRECISION
}

{
#if defined(GENERAL_INTEGER_PRECISION)
	int Amatrix_c4[DIAG][DIAG];
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c4[i][j] = (int) (Amatrix(i,j)*PGD_INT_SCALE);
		}
	}
	int bvector_c4[DIAG];
	for(int i=0;i<DIAG;i++){
		bvector_c4[i] = (int) (bvector(i)*PGD_INT_SCALE);
	}
	int L_c_inv = (int) ((1/L)*PGD_INT_SCALE);
	PROXIMAL_GRADIENT_DECENT_INTEGER(Amatrix_c4, bvector_c4, L_c_inv);
#endif
}

{
#if defined(COMSTOM_FLOAT_PRECISION)
#if DIAG < 1024
	fptx2 Amatrix_c3[DIAG][DIAG];
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c3[i][j] = (fptx2) Amatrix(i,j);
		}
	}
	fptx2 bvector_c3[DIAG];
	for(int i=0;i<DIAG;i++){
		bvector_c3[i] = (fptx2) bvector(i);
	}
	fptx2 L_c3 = (fptx2) L;
	PROXIMAL_GRADIENT_DECENT_XFLOAT2(Amatrix_c3, bvector_c3, L_c3);
#else
	fptx2 **Amatrix_c3;
	Amatrix_c3 = (fptx2**) malloc(sizeof(fptx2*)*DIAG);
	for(int i=0;i<DIAG;i++)
		Amatrix_c3[i] = (fptx2*) malloc(sizeof(fptx2)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Amatrix_c3[i][j] = (fptx2) Amatrix(i,j);
		}
	}
	fptx2 *bvector_c3;
	bvector_c3 = (fptx2*) malloc(sizeof(fptx2)*DIAG);
	for(int i=0;i<DIAG;i++){
		bvector_c3[i] = (fptx2) bvector(i);
	}
	fptx2 L_c3 = (fptx2) L;
	PROXIMAL_GRADIENT_DECENT_XFLOAT2(Amatrix_c3, bvector_c3, L_c3);
#endif
#endif
}

{
#ifdef XILINX_FIXED_PRECISION
#if DIAG < 1024
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
#else
	DATA_IN_T **Amatrix_c3;
	Amatrix_c3 = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	for(int i=0;i<DIAG;i++)
		Amatrix_c3[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			double tmp = Amatrix(i,j);
			Amatrix_c3[i][j] = (DATA_IN_T) tmp;
		}
	}
	DATA_IN_T *bvector_c3;
	bvector_c3 = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	for(int i=0;i<DIAG;i++){
		double tmp = bvector(i);
		bvector_c3[i] = (DATA_IN_T) tmp;
	}
	DATA_IN_T L_c3 = (DATA_IN_T) L;
	DATA_IN_T factor = 1/L;
	PROXIMAL_GRADIENT_DECENT_XFXPT(Amatrix_c3, bvector_c3, factor);
#endif
#endif // endif XILINX_FIXED_PRECISION
}

{
#if defined (SOFT_POSIT_PRECISION)
#if DIAG < 1
	double Amatrix_d[DIAG][DIAG];
	memcpy(Amatrix_d,AmatrixT.data(),sizeof(double)*DIAG*DIAG);
	double bvector_d[DIAG];
	memcpy(bvector_d,bvector.data(),sizeof(double)*DIAG);
	{
#if SOFT_POSIT_PRECISION==0
		posit32_t Amatrix_c4[DIAG][DIAG];
		posit32_t bvector_c4[DIAG];
		posit32_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP32(Amatrix_d[i][j]);

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP32(bvector_d[i]);
		L_c4 = convertDoubleToP32(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT32(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==1
		posit16_t Amatrix_c4[DIAG][DIAG];
		posit16_t bvector_c4[DIAG];
		posit16_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP16(Amatrix_d[i][j]);

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP16(bvector_d[i]);
		L_c4 = convertDoubleToP16(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT16(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==2
		posit8_t Amatrix_c4[DIAG][DIAG];
		posit8_t bvector_c4[DIAG];
		posit8_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP8(Amatrix_d[i][j]);

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP8(bvector_d[i]);
		L_c4 = convertDoubleToP8(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT8(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==3
		posit_2_t Amatrix_c4[DIAG][DIAG];
		posit_2_t bvector_c4[DIAG];
		posit_2_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = pX1_to_pX2(convertDoubleToPX1(Amatrix_d[i][j], TOTALBITS), TOTALBITS);

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = pX1_to_pX2(convertDoubleToPX1(bvector_d[i], TOTALBITS), TOTALBITS);
		L_c4 = pX1_to_pX2(convertDoubleToPX1(L, TOTALBITS), TOTALBITS);
		PROXIMAL_GRADIENT_DECENT_SPOSITX(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
#else
//	double **Amatrix_d;
//	double *bvector_d;
//	Amatrix_d = (double**) malloc(sizeof(double*)*DIAG);
//	for(int i=0;i<DIAG;i++)
//		Amatrix_d[i] = (double*) malloc(sizeof(double)*DIAG);
//	bvector_d = (double*) malloc(sizeof(double)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		for(int j=0;j<DIAG;j++){
//			Amatrix_d[i][j] = Amatrix(i,j);
//		}
//		bvector_d[i] = bvector(i);
//	}
	{
#if SOFT_POSIT_PRECISION==0
		posit32_t **Amatrix_c4;
		Amatrix_c4 = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
		for(int i=0;i<DIAG;i++)
			Amatrix_c4[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		posit32_t *bvector_c4;
		bvector_c4 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		posit32_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP32((double)Amatrix(i,j));

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP32((double)bvector(i));
		L_c4 = convertDoubleToP32(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT32(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==1
		posit16_t **Amatrix_c4;
		Amatrix_c4 = (posit16_t**) malloc(sizeof(posit16_t*)*DIAG);
		for(int i=0;i<DIAG;i++)
			Amatrix_c4[i] = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
		posit16_t *bvector_c4;
		bvector_c4 = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
		posit16_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP16((double)Amatrix(i,j));

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP16((double)bvector(i));
		L_c4 = convertDoubleToP16(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT16(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==2
		posit8_t **Amatrix_c4;
		Amatrix_c4 = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
		for(int i=0;i<DIAG;i++)
			Amatrix_c4[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		posit8_t *bvector_c4;
		bvector_c4 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		posit8_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = convertDoubleToP8((double)Amatrix(i,j));

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = convertDoubleToP8((double)bvector(i));
		L_c4 = convertDoubleToP8(L);
		PROXIMAL_GRADIENT_DECENT_SPOSIT8(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
	{
#if SOFT_POSIT_PRECISION==3
		posit_2_t **Amatrix_c4;
		Amatrix_c4 = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
		for(int i=0;i<DIAG;i++)
			Amatrix_c4[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		posit_2_t *bvector_c4;
		bvector_c4 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		posit_2_t L_c4;
		for(int i=0;i<DIAG;i++)
			for(int j=0;j<DIAG;j++)
				Amatrix_c4[i][j] = pX1_to_pX2(convertDoubleToPX1((double)Amatrix(i,j), TOTALBITS), TOTALBITS);

		for(int i=0;i<DIAG;i++)
			bvector_c4[i] = pX1_to_pX2(convertDoubleToPX1((double)bvector(i), TOTALBITS), TOTALBITS);
		L_c4 = pX1_to_pX2(convertDoubleToPX1(L, TOTALBITS), TOTALBITS);
		PROXIMAL_GRADIENT_DECENT_SPOSITX(Amatrix_c4, bvector_c4, L_c4);
#endif
	}
#endif
#endif
}
	///////////////////////////////////////////////////////////
}

int main(int argc, char** argv){
	std::string path = std::string(argv[1]);
	PROXIMAL_GRADIENT_DECENT(path)
	return 0;
}
