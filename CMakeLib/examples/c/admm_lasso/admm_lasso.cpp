/*
 * admm_lasso.cpp
 *
 *  Created on: 17 Sep 2019
 *      Author: yw106
 */

#include "admmlasso_gdouble.hpp"
#include "admmlasso_gfloat.hpp"
#include "admmlasso_xfxpt.hpp"
#include "admmlasso_xfloat.hpp"
#include "admmlasso_posit8.hpp"
#include "admmlasso_posit16.hpp"
#include "admmlasso_posit32.hpp"
#include "admmlasso_positx.hpp"

//void ADMM_LASSO_INTEGER(int A[DIAG][DIAG],
//					  int b[DIAG], int lambda,
//					  int rho, int alpha){
//
//	// parameters
//	int oneminusalpha = 1 - alpha;
//	int lambdadivrho = lambda / rho;
//
//	// variables
//	Float_Point_Algebra Float_Point_Algebraobj;
//	int Atb[DIAG];
//	int At[DIAG][DIAG];
//	int AtA[DIAG][DIAG];
//	int EYE[DIAG][DIAG];
//	int rhoEYE[DIAG][DIAG];
//	int AtAplusrhoeye[DIAG][DIAG];
//	int L[DIAG][DIAG], U[DIAG][DIAG];
//	int invL[DIAG][DIAG], invU[DIAG][DIAG];
//	int x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
//	int zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
//	int invLq[DIAG];
//	int alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
//	int x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
//	int x_hatz[DIAG];
//
//	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<int,DIAG,DIAG>(A, At);
//#if defined(DEBUG_ITER)
//	std::cout << "At:" << std::endl;
//	printmatrix(At);
//#endif
//	Float_Point_Algebraobj.MAT_VEC_MUL<int,DIAG,DIAG>(At, b, Atb);
//#if defined(DEBUG_ITER)
//	std::cout << "Atb:" << std::endl;
//	printvector(Atb);
//#endif
//	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<int,DIAG,DIAG>(At, A, AtA);
//#if defined(DEBUG_ITER)
//	std::cout << "AtA:" << std::endl;
//	printmatrix(AtA);
//#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<int,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<int,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<int,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
//#if defined(DEBUG_ITER)
//	std::cout << "AtAplusrhoeye:" << std::endl;
//	printmatrix(AtAplusrhoeye);
//#endif
//	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<int,DIAG>(
//			AtAplusrhoeye, L, U);
//#if defined(DEBUG_ITER)
//	std::cout << "L:" << std::endl;
//	printmatrix(L);
//	std::cout << "U:" << std::endl;
//	printmatrix(U);
//#endif
//	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<int,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<int,DIAG,DIAG>(invL, invU);
//	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//
//	/// iteration
//	Float_Point_Algebraobj.ZEROS_VEC<int,DIAG>(x);
//	Float_Point_Algebraobj.ZEROS_VEC<int,DIAG>(z);
//	Float_Point_Algebraobj.ZEROS_VEC<int,DIAG>(u);
//#ifdef TIME_PROFILE
//	clock_t start = clock();
//	std::ofstream TimeProfile;
//	std::string clockname = "timeprofile.txt";
//	TimeProfile.open(clockname);
//#endif
//	struct history<int> hist; // iteration record and early termination
//	for(int k=0;k<MAX_ITER;k++){
//		// q = Atb + rho*(z - u);
//		Float_Point_Algebraobj.VEC_SUB<int,DIAG>(z, u, zminusu);
//		Float_Point_Algebraobj.VEC_SCALAR_MUL<int,DIAG>(zminusu,
//				rho, rhozminusu);
//		Float_Point_Algebraobj.VEC_ADD<int,DIAG>(Atb,
//				rhozminusu, q);
//#if defined(DEBUG_ITER)
//		std::cout << "q:" << std::endl;
//		printvector(q);
//#endif
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//		// x = U \ (L \ q);
//		Float_Point_Algebraobj.MAT_VEC_MUL<int,DIAG,DIAG>(
//				invL, q, invLq);
//#if defined(DEBUG_ITER)
//		std::cout << "invL:" << std::endl;
//		printmatrix(invL);
//		std::cout << "invU:" << std::endl;
//		printmatrix(invU);
//		std::cout << "invLq:" << std::endl;
//		printvector(invLq);
//#endif
//		Float_Point_Algebraobj.MAT_VEC_MUL<int,DIAG,DIAG>(
//				invU, invLq, x);
//#if defined(DEBUG_ITER)
//		std::cout << "x:" << std::endl;
//		printvector(x);
//#endif
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//		// zold = z
//		Float_Point_Algebraobj.VEC_EQ<int,DIAG>(z, zold);
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//		//  x_hat = alpha*x + (1 - alpha)*zold;
//		Float_Point_Algebraobj.VEC_SCALAR_MUL<int,DIAG>(
//				x, alpha, alphax);
//		Float_Point_Algebraobj.VEC_SCALAR_MUL<int,DIAG>(
//				zold, oneminusalpha, oneminusalphazold);
//		Float_Point_Algebraobj.VEC_ADD<int,DIAG>(
//				alphax, oneminusalphazold, x_hat);
//#if defined(DEBUG_ITER)
//		std::cout << "x_hat:" << std::endl;
//		printvector(x_hat);
//#endif
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//		// z = shrinkage(x_hat + u, lambda/rho)
//		// 			shrinkage(x, kappa):
//		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
//		Float_Point_Algebraobj.VEC_ADD<int,DIAG>(x_hat, u, x_hatu);
//#if defined(DEBUG_ITER)
//		std::cout << "xhatu:" << std::endl;
//		printvector(x_hatu);
//#endif
//		Float_Point_Algebraobj.VEC_SCALAR_SUB<int,DIAG>(
//				x_hatu, lambdadivrho, x_hatu1);
//		Float_Point_Algebraobj.VEC_SCALAR_ADD<int,DIAG>(
//				x_hatu, lambdadivrho, x_hatu2);
//		Float_Point_Algebraobj.VEC_MINUS<int,DIAG>(
//				x_hatu2, x_hatu2);
//		Float_Point_Algebraobj.VEC_SCALAR_MAX<int,DIAG>(
//				x_hatu1, 0, x_hatu1);
//		Float_Point_Algebraobj.VEC_SCALAR_MAX<int,DIAG>(
//				x_hatu2, 0, x_hatu2);
//#if defined(DEBUG_ITER)
//		std::cout << "xhatu1:" << std::endl;
//		printvector(x_hatu1);
//		std::cout << "xhatu2:" << std::endl;
//		printvector(x_hatu2);
//#endif
//		Float_Point_Algebraobj.VEC_SUB<int,DIAG>(x_hatu1, x_hatu2, z);
//#if defined(DEBUG_ITER)
//		std::cout << "z:" << std::endl;
//		printvector(z);
//#endif
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//		// u = u + (x_hat - z);
//		Float_Point_Algebraobj.VEC_SUB<int,DIAG>(x_hat, z, x_hatz);
//#if defined(DEBUG_ITER)
//		std::cout << "x_hatz:" << std::endl;
//		printvector(x_hatz);
//#endif
//		Float_Point_Algebraobj.VEC_ADD<int,DIAG>(u, x_hatz, u);
//#if defined(DEBUG_ITER)
//		std::cout << "u:" << std::endl;
//		printvector(u);
//#endif
//#if defined(RECORD_RESULT)
//		int znorm;
//		int Ax[DIAG], Axb[DIAG];
//		int Axbnorm2;
//		int xz[DIAG], rhoxz[DIAG];
//		int xznorm, rhoxznorm;
//		int xnorm;
//		int rhou[DIAG];
//		int rhounorm;
//
//		// history.objval(k)  = objective(A, b, lambda, x, z);
//		// p = objective(A, b, lambda, x, z)
//	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
//		Float_Point_Algebraobj.VEC_NORM<int, DIAG>(z, znorm);
//		Float_Point_Algebraobj.MAT_VEC_MUL<int, DIAG, DIAG>(A, x, Ax);
//		Float_Point_Algebraobj.VEC_SUB<int, DIAG>(Ax, b, Axb);
//		Float_Point_Algebraobj.VEC_NORM2<int, DIAG>(Axb, Axbnorm2);
//		hist.objval[k] = 0.5 * Axbnorm2 + lambda * znorm;
//		// history.r_norm(k)  = norm(x - z);
//		Float_Point_Algebraobj.VEC_SUB<int, DIAG>(x, z, xz);
//		Float_Point_Algebraobj.VEC_NORM<int, DIAG>(xz, xznorm);
//		hist.r_norm[k] = xznorm;
//		// history.s_norm(k)  = norm(-rho*(z - zold));
//		Float_Point_Algebraobj.VEC_SCALAR_MUL<int, DIAG>(xz, rho, rhoxz);
//		Float_Point_Algebraobj.VEC_NORM<int, DIAG>(rhoxz, rhoxznorm);
//		hist.s_norm[k] = rhoxznorm;
//		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
//		Float_Point_Algebraobj.VEC_NORM<int, DIAG>(x, xnorm);
//		hist.eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
//		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
//		Float_Point_Algebraobj.VEC_SCALAR_MUL<int, DIAG>(rho, u, rhou);
//		Float_Point_Algebraobj.VEC_NORM<int, DIAG>(rhou, rhounorm);
//		hist.eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;
//		// record iterative solution
//		for(int i=0;i<DIAG;i++){
//			hist.u_hist[i][k] = u[i];
//			hist.x_hist[i][k] = x[i];
//			hist.z_hist[i][k] = z[i];
//		}
//#if defined(EARLY_TERMINATE)
//		if((hist.r_norm[k] < hist.eps_pri[k]) &&
//		  (hist.s_norm[k] < hist.eps_dual[k])){
//			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
//			break;
//		}
//#endif
//#endif
//		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
//
//	}
//#if defined(DEBUG_ITER)
//	std::cout << "final x:" << std::endl;
//	printvector(x);
//#endif
//#if defined(RECORD_RESULT)
//	std::string xkname = "xk_int.dat";
//	std::string ukname = "uk_int.dat";
//	std::string zkname = "zk_int.dat";
//	std::ofstream resultfile(xkname);
//	std::ofstream resultfile1(ukname);
//	std::ofstream resultfile2(zkname);
//	for(int i=0; i<MAX_ITER; i++){
//		for(int j=0;j<DIAG;j++){
//			resultfile << hist.x_hist[j][i] << ",";
//			resultfile1 << hist.u_hist[j][i] << ",";
//			resultfile2 << hist.z_hist[j][i] << ",";
//		}
//		resultfile << "\n";
//		resultfile1 << "\n";
//		resultfile2 << "\n";
//	}
//	resultfile.close();
//	resultfile1.close();
//	resultfile2.close();
//#endif
//#ifdef TIME_PROFILE
//	clock_t end = clock();
//	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
//	std::cout << "It takes "
//			  << time
//			  << " ms to finish the iteration"
//			  << std::endl << std::endl;
//	TimeProfile << "It takes "
//			  << time
//			  << " ms to finish the iteration"
//			  << std::endl << std::endl;
//#endif// endif TIME_PROFILE
//
//}

void ADMM_LASSO(std::string path){

	std::string Amatrixname = path+"A.csv";
	std::string bvectorname = path+"b.csv";
	std::string lambdaname = path+"lambda.csv";
	std::string Atmatrixname = path+"At.csv";
	std::string Umatrixname = path+"U.csv";
	std::string Lmatrixname = path+"L.csv";
	std::string Uinvmatrixname = path+"U_inv.csv";
	std::string Linvmatrixname = path+"L_inv.csv";
	std::string ULinvmatrixname = path+"UL_inv.csv";

	std::string clockname = "timeprofile.txt";

#ifdef ALWAYS_DELETE
	std::remove(clockname.c_str());
	std::remove(Amatrixname.c_str());
	std::remove(bvectorname.c_str());
	std::remove(lambdaname.c_str());
	std::remove(rhoname.c_str());
	std::remove(alphaname.c_str());
#endif
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif
#if defined(ALWAYS_DELETE)
	std::remove(clockname.c_str());
	std::remove(Amatrixname.c_str());
	std::remove(bvectorname.c_str());
	std::remove(lambdaname.c_str());
#endif
	Eigen_Algebra Eigen_Algebra_obj;
	Eigen::MatrixXd Amatrix( DIAG, DIAG );
	Eigen::MatrixXd AmatrixT( DIAG, DIAG );
	Eigen::MatrixXd bvector( DIAG, 1 );
	Eigen::MatrixXd U( DIAG, DIAG ), L( DIAG, DIAG );
	Eigen::MatrixXd U_inv( DIAG, DIAG ), L_inv( DIAG, DIAG );
	Eigen::MatrixXd U_invT( DIAG, DIAG ), L_invT( DIAG, DIAG );
	double lambda;
	if(!checkfileexist(Amatrixname)){
		Eigen::MatrixXd RNDmatrix( ROW, COL );
		// A = randn(m,n);
		Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, DIAG, DIAG>( RNDmatrix );
		// A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n);
		Eigen::VectorXd Vec = RNDmatrix.cwiseAbs2().colwise().sum().cwiseSqrt();
		Vec.resize(Vec.cols()*Vec.rows(), 1);
		for(int i=0;i<Vec.size();i++){
			Vec(i) = 1/Vec(i);
		}
		Eigen::MatrixXd Diagmatrix = Vec.matrix().asDiagonal();
		Amatrix = RNDmatrix * Diagmatrix;
		// write to data file
		writematrix_double(Amatrix, Amatrixname, 64);
		std::cout << "created the Amatrix" << std::endl << std::endl;
	}else{
		Amatrix = readmatrix_double(Amatrixname);
		std::cout << "read the Amatrix" << std::endl << std::endl;
	}
	AmatrixT = Amatrix.transpose();
	if(!checkfileexist(bvectorname)){
		// b = A*x0 + sqrt(0.001)*randn(m,1);
		Eigen::MatrixXd x0(DIAG, 1);
		Eigen::MatrixXd tmprndmat(DIAG, 1);
		Eigen_Algebra_obj.RND_SPVEC<Eigen::MatrixXd, DIAG>( x0 );
		Eigen_Algebra_obj.RND_MAT<Eigen::MatrixXd, DIAG, 1>( tmprndmat );
		tmprndmat = std::sqrt(0.001) * tmprndmat;
		bvector = Amatrix * x0 + tmprndmat;
		// write to data file
		writematrix_double(bvector, bvectorname, 64);
		std::cout << "created the bvector" << std::endl << std::endl;
	}else{
		bvector = readmatrix_double(bvectorname);
		std::cout << "read the bvector" << std::endl << std::endl;
	}
	if(!checkfileexist(lambdaname)){
		// lambda_max = norm( A'*b, 'inf' );
		// lambda = 0.1*lambda_max;
		Eigen::MatrixXd Atb = AmatrixT * bvector;
		double lambda_max = Atb.norm();
		lambda = 0.1*lambda_max;
		// write to data file
		writescalar_double(lambda, lambdaname, 64);
		std::cout << "created the lambda value" << std::endl << std::endl;
	}else{
		// read from data file
		lambda = readscalar_double(lambdaname);
		std::cout << "read the lambda value" << std::endl << std::endl;
	}
	double rho=1.0, alpha=1.0;
	Eigen::MatrixXd AtA = Amatrix * AmatrixT;
	Eigen::MatrixXd AtAplusrhoeye = AtA * rho;
	Eigen::FullPivLU<Eigen::MatrixXd> LU(AtAplusrhoeye);
	if(!checkfileexist(Umatrixname)||!checkfileexist(Lmatrixname)){
		U = LU.matrixLU().triangularView<Eigen::Upper>();
		L = LU.matrixLU().triangularView<Eigen::Lower>();
		// write to data file
		writematrix_double(U, Umatrixname, 64);
		writematrix_double(L, Lmatrixname, 64);
		std::cout << "created the Umatrix" << std::endl << std::endl;
		std::cout << "created the Lmatrix" << std::endl << std::endl;
	}else{
		U = readmatrix_double(Umatrixname);
		L = readmatrix_double(Lmatrixname);
		std::cout << "read the Umatrix" << std::endl << std::endl;
		std::cout << "read the Lmatrix" << std::endl << std::endl;
	}
	if(!checkfileexist(Uinvmatrixname)||!checkfileexist(Linvmatrixname)){
		U_inv = U.inverse();
		L_inv = L.inverse();
		// write to data file
		writematrix_double(U_inv, Uinvmatrixname, 64);
		writematrix_double(L_inv, Linvmatrixname, 64);
		std::cout << "created the Uinvmatrix" << std::endl << std::endl;
		std::cout << "created the Linvmatrix" << std::endl << std::endl;
	}else{
		U_inv = readmatrix_double(Uinvmatrixname);
		L_inv = readmatrix_double(Linvmatrixname);
		std::cout << "read the Uinvmatrix" << std::endl << std::endl;
		std::cout << "read the Linvmatrix" << std::endl << std::endl;
	}
	U_invT = U_inv.transpose();
	L_invT = L_inv.transpose();
//	std::cout << Amatrix << std::endl << std::endl;
//	std::cout << AmatrixT << std::endl << std::endl;
//	std::cout << U << std::endl << std::endl;
//	std::cout << L << std::endl << std::endl;
//	std::cout << U_inv << std::endl << std::endl;
//	std::cout << L_inv << std::endl << std::endl;
//	std::cout << U_invT << std::endl << std::endl;
//	std::cout << L_invT << std::endl << std::endl;

#if defined(DEBUG_DATA)
	std::cout << "A:" << std::endl << Amatrix << std::endl << std::endl;
	std::cout << "b:" << std::endl << bvector << std::endl << std::endl;
	std::cout << "lambda:" << std::endl << lambda << std::endl << std::endl;
#endif

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

#if defined(GENERAL_DOUBLE_PRECISION)
#if DIAG < 1024
	double A_c[DIAG][DIAG], b_c[DIAG];
	double Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	double lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (double) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (double) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (double) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (double) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (double) bvector(i);
	}
	lambda_c = (double) lambda;
	rho_c = (double) rho;
	alpha_c = (double) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_DOUBLE(A_c, At_c, Linv_c, Uinv_c, b_c,
			lambda_c, rho_c, alpha_c);
#else
	double **A_c, *b_c;
	double **Uinv_c, **Linv_c, **At_c;
	double lambda_c, rho_c, alpha_c;
	A_c = (double**) malloc(sizeof(double*)*DIAG);
	Uinv_c = (double**) malloc(sizeof(double*)*DIAG);
	Linv_c = (double**) malloc(sizeof(double*)*DIAG);
	At_c = (double**) malloc(sizeof(double*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (double*) malloc(sizeof(double)*DIAG);
		Uinv_c[i] = (double*) malloc(sizeof(double)*DIAG);
		Linv_c[i] = (double*) malloc(sizeof(double)*DIAG);
		At_c[i] = (double*) malloc(sizeof(double)*DIAG);
	}
	b_c = (double*) malloc(sizeof(double)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (double) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (double) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (double) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (double) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (double) bvector(i);
	}
	lambda_c = (double) lambda;
	rho_c = (double) rho;
	alpha_c = (double) alpha;
//	std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_DOUBLE(A_c, At_c, Linv_c, Uinv_c, b_c,
			lambda_c, rho_c, alpha_c);
//	ADMM_LASSO_DOUBLE(A_c, b_c, lambda_c, rho_c, alpha_c);
//	ADMM_LASSO_DOUBLE(A_c, b_c, lambda_c);
#endif
#endif
#if defined(GENERAL_FLOAT_PRECISION)
#if DIAG  < 1024
	Eigen::MatrixXf Amatrix_f(DIAG, DIAG);
	Eigen::MatrixXf AmatrixT_f(DIAG, DIAG);
	Eigen::MatrixXf U_invT_f(DIAG, DIAG);
	Eigen::MatrixXf L_invT_f(DIAG, DIAG);
	Eigen::VectorXf bvector_f(DIAG);
	Amatrix_f = Amatrix.cast<float>();
	AmatrixT_f = AmatrixT.cast<float>();
	U_invT_f = U_invT.cast<float>();
	L_invT_f = L_invT.cast<float>();
	bvector_f = bvector.cast<float>();
	float A_c[DIAG][DIAG], b_c[DIAG];
	float Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	float lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (float) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (float) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (float) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (float) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (float) bvector(i);
	}
	lambda_c = (float) lambda;
	rho_c = (float) rho;
	alpha_c = (float) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_FLOAT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	float **A_c, *b_c;
	float **Uinv_c, **Linv_c, **At_c;
	float lambda_c, rho_c, alpha_c;
	A_c = (float**) malloc(sizeof(float*)*DIAG);
	Uinv_c = (float**) malloc(sizeof(float*)*DIAG);
	Linv_c = (float**) malloc(sizeof(float*)*DIAG);
	At_c = (float**) malloc(sizeof(float*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (float*) malloc(sizeof(float)*DIAG);
		Uinv_c[i] = (float*) malloc(sizeof(float)*DIAG);
		Linv_c[i] = (float*) malloc(sizeof(float)*DIAG);
		At_c[i] = (float*) malloc(sizeof(float)*DIAG);
	}
	b_c = (float*) malloc(sizeof(float)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (float) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (float) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (float) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (float) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (float) bvector(i);
	}
	lambda_c = (float) lambda;
	rho_c = (float) rho;
	alpha_c = (float) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_FLOAT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if defined(COMSTOM_FLOAT_PRECISION)
#if DIAG < 1024
	fptx_admmlasso A_c[DIAG][DIAG], b_c[DIAG];
	fptx_admmlasso Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	fptx_admmlasso lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (fptx_admmlasso) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (fptx_admmlasso) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (fptx_admmlasso) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (fptx_admmlasso) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (fptx_admmlasso) bvector(i);
	}
	lambda_c = (fptx_admmlasso) lambda;
	rho_c = (fptx_admmlasso) rho;
	alpha_c = (fptx_admmlasso) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_XFPT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	fptx_admmlasso **A_c, *b_c;
	fptx_admmlasso **Uinv_c, **Linv_c, **At_c;
	fptx_admmlasso lambda_c, rho_c, alpha_c;
	A_c = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	Uinv_c = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	Linv_c = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	At_c = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		Uinv_c[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		Linv_c[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		At_c[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	}
	b_c = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (fptx_admmlasso) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (fptx_admmlasso) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (fptx_admmlasso) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (fptx_admmlasso) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (fptx_admmlasso) bvector(i);
	}
	lambda_c = (fptx_admmlasso) lambda;
	rho_c = (fptx_admmlasso) rho;
	alpha_c = (fptx_admmlasso) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_XFPT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if defined(XILINX_FIXED_PRECISION)
#if DIAG < 1024
	DATA_IN_T A_c[DIAG][DIAG], b_c[DIAG];
	DATA_IN_T Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	DATA_IN_T lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (DATA_IN_T) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (DATA_IN_T) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (DATA_IN_T) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (DATA_IN_T) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (DATA_IN_T) bvector(i);
	}
	lambda_c = (DATA_IN_T) lambda;
	rho_c = (DATA_IN_T) rho;
	alpha_c = (DATA_IN_T) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_FXPT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	DATA_IN_T **A_c, *b_c;
	DATA_IN_T **Uinv_c, **Linv_c, **At_c;
	DATA_IN_T lambda_c, rho_c, alpha_c;
	A_c = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	Uinv_c = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	Linv_c = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	At_c = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		Uinv_c[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		Linv_c[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		At_c[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	}
	b_c = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = (DATA_IN_T) Amatrix(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = (DATA_IN_T) AmatrixT(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = (DATA_IN_T) U_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = (DATA_IN_T) L_inv(i,j);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = (DATA_IN_T) bvector(i);
	}
	lambda_c = (DATA_IN_T) lambda;
	rho_c = (DATA_IN_T) rho;
	alpha_c = (DATA_IN_T) alpha;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_FXPT(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if defined(SOFT_POSIT_PRECISION)
#if SOFT_POSIT_PRECISION==0
#if DIAG < 1024
	posit32_t A_c[DIAG][DIAG], b_c[DIAG];
	posit32_t Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	posit32_t lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP32((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP32((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP32((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP32((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP32((double)bvector(i));
	}
	lambda_c = convertDoubleToP32(lambda);
	rho_c = convertDoubleToP32(rho);
	alpha_c = convertDoubleToP32(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT32(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	posit32_t **A_c, *b_c;
	posit32_t **Uinv_c, **Linv_c, **At_c;
	posit32_t lambda_c, rho_c, alpha_c;
	A_c = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	Uinv_c = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	Linv_c = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	At_c = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		Uinv_c[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		Linv_c[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		At_c[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	}
	b_c = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP32((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP32((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP32((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP32((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP32((double)bvector(i));
	}
	lambda_c = convertDoubleToP32(lambda);
	rho_c = convertDoubleToP32(rho);
	alpha_c = convertDoubleToP32(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT32(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if SOFT_POSIT_PRECISION==1
#if DIAG < 1024
	posit16_t A_c[DIAG][DIAG], b_c[DIAG];
	posit16_t Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	posit16_t lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP16((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP16((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP16((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP16((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP16((double)bvector(i));
	}
	lambda_c = convertDoubleToP16(lambda);
	rho_c = convertDoubleToP16(rho);
	alpha_c = convertDoubleToP16(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT16(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	posit16_t **A_c, *b_c;
	posit16_t **Uinv_c, **Linv_c, **At_c;
	posit16_t lambda_c, rho_c, alpha_c;
	A_c = (posit16_t**) malloc(sizeof(posit16_t*)*DIAG);
	Uinv_c = (posit16_t**) malloc(sizeof(posit16_t*)*DIAG);
	Linv_c = (posit16_t**) malloc(sizeof(posit16_t*)*DIAG);
	At_c = (posit16_t**) malloc(sizeof(posit16_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
		Uinv_c[i] = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
		Linv_c[i] = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
		At_c[i] = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
	}
	b_c = (posit16_t*) malloc(sizeof(posit16_t)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP16((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP16((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP16((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP16((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP16((double)bvector(i));
	}
	lambda_c = convertDoubleToP16(lambda);
	rho_c = convertDoubleToP16(rho);
	alpha_c = convertDoubleToP16(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT16(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if SOFT_POSIT_PRECISION==2
#if DIAG < 1024
	posit8_t A_c[DIAG][DIAG], b_c[DIAG];
	posit8_t Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	posit8_t lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP8((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP8((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP8((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP8((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP8((double)bvector(i));
	}
	lambda_c = convertDoubleToP8(lambda);
	rho_c = convertDoubleToP8(rho);
	alpha_c = convertDoubleToP8(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT8(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	posit8_t **A_c, *b_c;
	posit8_t **Uinv_c, **Linv_c, **At_c;
	posit8_t lambda_c, rho_c, alpha_c;
	A_c = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	Uinv_c = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	Linv_c = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	At_c = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		Uinv_c[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		Linv_c[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		At_c[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	}
	b_c = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = convertDoubleToP8((double)Amatrix(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = convertDoubleToP8((double)AmatrixT(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = convertDoubleToP8((double)U_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = convertDoubleToP8((double)L_inv(i,j));
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = convertDoubleToP8((double)bvector(i));
	}
	lambda_c = convertDoubleToP8(lambda);
	rho_c = convertDoubleToP8(rho);
	alpha_c = convertDoubleToP8(alpha);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSIT8(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#if SOFT_POSIT_PRECISION==3
#if DIAG < 1024
	posit_2_t A_c[DIAG][DIAG], b_c[DIAG];
	posit_2_t Uinv_c[DIAG][DIAG], Linv_c[DIAG][DIAG], At_c[DIAG][DIAG];
	posit_2_t lambda_c, rho_c, alpha_c;
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)Amatrix(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)AmatrixT(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)U_inv(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)L_inv(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = pX1_to_pX2(convertDoubleToPX1((double)bvector(i),TOTALBITS),TOTALBITS);
	}
	lambda_c = pX1_to_pX2(convertDoubleToPX1(lambda,TOTALBITS),TOTALBITS);
	rho_c = pX1_to_pX2(convertDoubleToPX1(rho,TOTALBITS),TOTALBITS);
	alpha_c = pX1_to_pX2(convertDoubleToPX1(alpha,TOTALBITS),TOTALBITS);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSITX(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#else
	posit_2_t **A_c, *b_c;
	posit_2_t **Uinv_c, **Linv_c, **At_c;
	posit_2_t lambda_c, rho_c, alpha_c;
	A_c = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	Uinv_c = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	Linv_c = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	At_c = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		A_c[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		Uinv_c[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		Linv_c[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		At_c[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	}
	b_c = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			A_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)Amatrix(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			At_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)AmatrixT(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Uinv_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)U_inv(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			Linv_c[i][j] = pX1_to_pX2(convertDoubleToPX1((double)L_inv(i,j),TOTALBITS),TOTALBITS);
		}
	}
	for(int i=0;i<DIAG;i++){
		b_c[i] = pX1_to_pX2(convertDoubleToPX1((double)bvector(i),TOTALBITS),TOTALBITS);
	}
	lambda_c = pX1_to_pX2(convertDoubleToPX1(lambda,TOTALBITS),TOTALBITS);
	rho_c = pX1_to_pX2(convertDoubleToPX1(rho,TOTALBITS),TOTALBITS);
	alpha_c = pX1_to_pX2(convertDoubleToPX1(alpha,TOTALBITS),TOTALBITS);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_POSITX(A_c, At_c, Linv_c, Uinv_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#endif
#endif
}

int main(int argc, char** argv){
	std::string path = std::string(argv[1]);
	ADMM_LASSO(path);
	return 0;
}
	
