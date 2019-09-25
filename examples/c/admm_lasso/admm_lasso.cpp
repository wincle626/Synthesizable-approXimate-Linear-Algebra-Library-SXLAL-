/*
 * admm_lasso.cpp
 *
 *  Created on: 17 Sep 2019
 *      Author: yw106
 */

#include "admm_lasso.hpp"

#define QUIET      0
#define ABSTOL     1e-4
#define RELTOL     1e-2
#define MAX_ITER   10

template<typename dtype>
struct history{
	dtype u_hist[DIAG][MAX_ITER];
	dtype x_hist[DIAG][MAX_ITER];
	dtype z_hist[DIAG][MAX_ITER];
	dtype objval[MAX_ITER];
	dtype r_norm[MAX_ITER];
	dtype s_norm[MAX_ITER];
	dtype eps_pri[MAX_ITER];
	dtype eps_dual[MAX_ITER];
};

struct histry{
	double u_hist[DIAG][MAX_ITER];
	double x_hist[DIAG][MAX_ITER];
	double z_hist[DIAG][MAX_ITER];
	double objval[MAX_ITER];
	double r_norm[MAX_ITER];
	double s_norm[MAX_ITER];
	double eps_pri[MAX_ITER];
	double eps_dual[MAX_ITER];
};

void printmatrix(double M[DIAG][DIAG]){
	for(int i=0;i<DIAG;i++){
		for(int j=0;j<DIAG;j++){
			std::cout << M[i][j] << ",";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printvector(double V[DIAG]){
	for(int j=0;j<DIAG;j++){
		std::cout << V[j] << ",";
	}
	std::cout << std::endl << std::endl;
}

void ADMM_LASSO_DOUBLE(double A[DIAG][DIAG],
					   double b[DIAG], double lambda,
					   double rho, double alpha){

	// parameters
	double oneminusalpha = 1 - alpha;
	double lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	double Atb[DIAG];
	double At[DIAG][DIAG];
	double AtA[DIAG][DIAG];
	double EYE[DIAG][DIAG];
	double rhoEYE[DIAG][DIAG];
	double AtAplusrhoeye[DIAG][DIAG];
	double L[DIAG][DIAG], U[DIAG][DIAG];
	double invL[DIAG][DIAG], invU[DIAG][DIAG];
	double x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	double zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	double invLq[DIAG];
	double alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	double x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	double x_hatz[DIAG];

	// A'*b
	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Float_Point_Algebraobj.MAT_MUL<double,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	Float_Point_Algebraobj.IDENDTITY_MAT<double,DIAG,DIAG>(EYE);
	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<double,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Float_Point_Algebraobj.MAT_ADD<double,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Float_Point_Algebraobj.LU_CHOLBANACHROUT<double,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	Float_Point_Algebraobj.MAT_QRINV<double,DIAG>(L, invL);
	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	Float_Point_Algebraobj.ZEROS_VEC<double,DIAG>(x);
	Float_Point_Algebraobj.ZEROS_VEC<double,DIAG>(z);
	Float_Point_Algebraobj.ZEROS_VEC<double,DIAG>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct histry hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		Float_Point_Algebraobj.VEC_EQ<double,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(
				x, alpha, alphax);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(
				zold, oneminusalpha, oneminusalphazold);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		Float_Point_Algebraobj.VEC_SCALAR_SUB<double,DIAG>(
				x_hatu, lambdadivrho, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_ADD<double,DIAG>(
				x_hatu, lambdadivrho, x_hatu2);
		Float_Point_Algebraobj.VEC_MINUS<double,DIAG>(
				x_hatu2, x_hatu2);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<double,DIAG>(
				x_hatu1, 0, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<double,DIAG>(
				x_hatu2, 0, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		double znorm;
		double Ax[DIAG], Axb[DIAG];
		double Axbnorm2;
		double xz[DIAG], rhoxz[DIAG];
		double xznorm, rhoxznorm;
		double xnorm;
		double rhou[DIAG];
		double rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<double, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<double, DIAG>(Axb, Axbnorm2);
		hist.objval[k] = 0.5 * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(xz, rho, rhoxz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(x, xnorm);
		hist.eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhou, rhounorm);
		hist.eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if((hist.r_norm[k] < hist.eps_pri[k]) &&
		  (hist.s_norm[k] < hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gdouble.dat";
	std::string ukname = "uk_gdouble.dat";
	std::string zkname = "zk_gdouble.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << hist.x_hist[j][i] << ",";
			resultfile1 << hist.u_hist[j][i] << ",";
			resultfile2 << hist.z_hist[j][i] << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE
}

void ADMM_LASSO_FLOAT(float A[DIAG][DIAG],
					  float b[DIAG], float lambda,
					  float rho, float alpha){

	// parameters
	float oneminusalpha = 1 - alpha;
	float lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	float Atb[DIAG];
	float At[DIAG][DIAG];
	float AtA[DIAG][DIAG];
	float EYE[DIAG][DIAG];
	float rhoEYE[DIAG][DIAG];
	float AtAplusrhoeye[DIAG][DIAG];
	float L[DIAG][DIAG], U[DIAG][DIAG];
	float invL[DIAG][DIAG], invU[DIAG][DIAG];
	float x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	float zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	float invLq[DIAG];
	float alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	float x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	float x_hatz[DIAG];

	// A'*b
	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Float_Point_Algebraobj.MAT_MUL<float,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	Float_Point_Algebraobj.IDENDTITY_MAT<float,DIAG,DIAG>(EYE);
	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<float,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Float_Point_Algebraobj.MAT_ADD<float,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Float_Point_Algebraobj.LU_CHOLBANACHROUT<float,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	Float_Point_Algebraobj.MAT_QRINV<float,DIAG>(L, invL);
	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	Float_Point_Algebraobj.ZEROS_VEC<float,DIAG>(x);
	Float_Point_Algebraobj.ZEROS_VEC<float,DIAG>(z);
	Float_Point_Algebraobj.ZEROS_VEC<float,DIAG>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<float> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		Float_Point_Algebraobj.VEC_EQ<float,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float,DIAG>(
				x, alpha, alphax);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float,DIAG>(
				zold, oneminusalpha, oneminusalphazold);
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		Float_Point_Algebraobj.VEC_SCALAR_SUB<float,DIAG>(
				x_hatu, lambdadivrho, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_ADD<float,DIAG>(
				x_hatu, lambdadivrho, x_hatu2);
		Float_Point_Algebraobj.VEC_MINUS<float,DIAG>(
				x_hatu2, x_hatu2);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<float,DIAG>(
				x_hatu1, 0, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<float,DIAG>(
				x_hatu2, 0, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		float znorm;
		float Ax[DIAG], Axb[DIAG];
		float Axbnorm2;
		float xz[DIAG], rhoxz[DIAG];
		float xznorm, rhoxznorm;
		float xnorm;
		float rhou[DIAG];
		float rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<float, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<float, DIAG>(Axb, Axbnorm2);
		hist.objval[k] = 0.5 * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float, DIAG>(xz, rho, rhoxz);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(x, xnorm);
		hist.eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(rhou, rhounorm);
		hist.eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if((hist.r_norm[k] < hist.eps_pri[k]) &&
		  (hist.s_norm[k] < hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gfloat.dat";
	std::string ukname = "uk_gfloat.dat";
	std::string zkname = "zk_gfloat.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << hist.x_hist[j][i] << ",";
			resultfile1 << hist.u_hist[j][i] << ",";
			resultfile2 << hist.z_hist[j][i] << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

void ADMM_LASSO_FXPT(DATA_IN_T A[DIAG][DIAG],
					 DATA_IN_T b[DIAG], DATA_IN_T lambda,
					 DATA_IN_T rho, float alpha){

	// parameters
	DATA_IN_T oneminusalpha = 1 - alpha;
	DATA_IN_T lambdadivrho = lambda / rho;

	// variables
	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebraobj;
	DATA_IN_T Atb[DIAG];
	DATA_IN_T At[DIAG][DIAG];
	DATA_IN_T AtA[DIAG][DIAG];
	DATA_IN_T EYE[DIAG][DIAG];
	DATA_IN_T rhoEYE[DIAG][DIAG];
	DATA_IN_T AtAplusrhoeye[DIAG][DIAG];
	DATA_IN_T L[DIAG][DIAG], U[DIAG][DIAG];
	DATA_IN_T invL[DIAG][DIAG], invU[DIAG][DIAG];
	DATA_IN_T x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	DATA_IN_T zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	DATA_IN_T invLq[DIAG];
	DATA_IN_T alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	DATA_IN_T x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	DATA_IN_T x_hatz[DIAG];

	// A'*b
	Xilinx_Fixed_Point_Algebraobj.MAT_TRANS<DATA_IN_T,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Xilinx_Fixed_Point_Algebraobj.MAT_MUL<DATA_IN_T,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	Xilinx_Fixed_Point_Algebraobj.IDENDTITY_MAT<DATA_IN_T,DIAG,DIAG>(EYE);
	Xilinx_Fixed_Point_Algebraobj.MAT_SCALAR_DOTMUL<DATA_IN_T,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Xilinx_Fixed_Point_Algebraobj.MAT_ADD<DATA_IN_T,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Xilinx_Fixed_Point_Algebraobj.LU_CHOLBANACHROUT<DATA_IN_T,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	Xilinx_Fixed_Point_Algebraobj.MAT_QRINV<DATA_IN_T,DIAG>(L, invL);
	Xilinx_Fixed_Point_Algebraobj.MAT_TRANS<DATA_IN_T,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	Xilinx_Fixed_Point_Algebraobj.ZEROS_VEC<DATA_IN_T,DIAG>(x);
	Xilinx_Fixed_Point_Algebraobj.ZEROS_VEC<DATA_IN_T,DIAG>(z);
	Xilinx_Fixed_Point_Algebraobj.ZEROS_VEC<DATA_IN_T,DIAG>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<DATA_IN_T> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(z, u, zminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T,DIAG>(zminusu,
				rho, rhozminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		Xilinx_Fixed_Point_Algebraobj.VEC_EQ<DATA_IN_T,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T,DIAG>(
				x, alpha, alphax);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T,DIAG>(
				zold, oneminusalpha, oneminusalphazold);
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_SUB<DATA_IN_T,DIAG>(
				x_hatu, lambdadivrho, x_hatu1);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_ADD<DATA_IN_T,DIAG>(
				x_hatu, lambdadivrho, x_hatu2);
		Xilinx_Fixed_Point_Algebraobj.VEC_MINUS<DATA_IN_T,DIAG>(
				x_hatu2, x_hatu2);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MAX<DATA_IN_T,DIAG>(
				x_hatu1, 0, x_hatu1);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MAX<DATA_IN_T,DIAG>(
				x_hatu2, 0, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		DATA_IN_T znorm;
		DATA_IN_T Ax[DIAG], Axb[DIAG];
		DATA_IN_T Axbnorm2;
		DATA_IN_T xz[DIAG], rhoxz[DIAG];
		DATA_IN_T xznorm, rhoxznorm;
		DATA_IN_T xnorm;
		DATA_IN_T rhou[DIAG];
		DATA_IN_T rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(z, znorm);
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T, DIAG, DIAG>(A, x, Ax);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(Ax, b, Axb);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM2<DATA_IN_T, DIAG>(Axb, Axbnorm2);
		hist.objval[k] = ((DATA_IN_T)0.5) * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(x, z, xz);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>(xz, rho, rhoxz);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(x, xnorm);
		hist.eps_pri[k] = ((DATA_IN_T)std::sqrt(DIAG))*((DATA_IN_T)ABSTOL)+((DATA_IN_T)RELTOL)*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>(rho, u, rhou);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(rhou, rhounorm);
		hist.eps_dual[k] = ((DATA_IN_T)std::sqrt(DIAG))*((DATA_IN_T)ABSTOL)+((DATA_IN_T)RELTOL)*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if((hist.r_norm[k] < hist.eps_pri[k]) &&
		  (hist.s_norm[k] < hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gDATA_IN_T.dat";
	std::string ukname = "uk_gDATA_IN_T.dat";
	std::string zkname = "zk_gDATA_IN_T.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << (float)hist.x_hist[j][i] << ",";
			resultfile1 << (float)hist.u_hist[j][i] << ",";
			resultfile2 << (float)hist.z_hist[j][i] << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

void ADMM_LASSO_XFPT(fptx_admmlasso A[DIAG][DIAG],
					 fptx_admmlasso b[DIAG], fptx_admmlasso lambda,
				 	 fptx_admmlasso rho, float alpha){

	// parameters
	fptx_admmlasso oneminusalpha = 1 - alpha;
	fptx_admmlasso lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	fptx_admmlasso Atb[DIAG];
	fptx_admmlasso At[DIAG][DIAG];
	fptx_admmlasso AtA[DIAG][DIAG];
	fptx_admmlasso EYE[DIAG][DIAG];
	fptx_admmlasso rhoEYE[DIAG][DIAG];
	fptx_admmlasso AtAplusrhoeye[DIAG][DIAG];
	fptx_admmlasso L[DIAG][DIAG], U[DIAG][DIAG];
	fptx_admmlasso invL[DIAG][DIAG], invU[DIAG][DIAG];
	fptx_admmlasso x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	fptx_admmlasso zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	fptx_admmlasso invLq[DIAG];
	fptx_admmlasso alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	fptx_admmlasso x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	fptx_admmlasso x_hatz[DIAG];

	// A'*b
	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Float_Point_Algebraobj.MAT_MUL<fptx_admmlasso,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	Float_Point_Algebraobj.IDENDTITY_MAT<fptx_admmlasso,DIAG,DIAG>(EYE);
	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<fptx_admmlasso,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Float_Point_Algebraobj.MAT_ADD<fptx_admmlasso,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Float_Point_Algebraobj.LU_CHOLBANACHROUT<fptx_admmlasso,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	Float_Point_Algebraobj.MAT_QRINV<fptx_admmlasso,DIAG>(L, invL);
	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	Float_Point_Algebraobj.ZEROS_VEC<fptx_admmlasso,DIAG>(x);
	Float_Point_Algebraobj.ZEROS_VEC<fptx_admmlasso,DIAG>(z);
	Float_Point_Algebraobj.ZEROS_VEC<fptx_admmlasso,DIAG>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<fptx_admmlasso> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		Float_Point_Algebraobj.VEC_EQ<fptx_admmlasso,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso,DIAG>(
				x, alpha, alphax);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso,DIAG>(
				zold, oneminusalpha, oneminusalphazold);
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		Float_Point_Algebraobj.VEC_SCALAR_SUB<fptx_admmlasso,DIAG>(
				x_hatu, lambdadivrho, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_ADD<fptx_admmlasso,DIAG>(
				x_hatu, lambdadivrho, x_hatu2);
		Float_Point_Algebraobj.VEC_MINUS<fptx_admmlasso,DIAG>(
				x_hatu2, x_hatu2);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<fptx_admmlasso,DIAG>(
				x_hatu1, 0, x_hatu1);
		Float_Point_Algebraobj.VEC_SCALAR_MAX<fptx_admmlasso,DIAG>(
				x_hatu2, 0, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		fptx_admmlasso znorm;
		fptx_admmlasso Ax[DIAG], Axb[DIAG];
		fptx_admmlasso Axbnorm2;
		fptx_admmlasso xz[DIAG], rhoxz[DIAG];
		fptx_admmlasso xznorm, rhoxznorm;
		fptx_admmlasso xnorm;
		fptx_admmlasso rhou[DIAG];
		fptx_admmlasso rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<fptx_admmlasso, DIAG>(Axb, Axbnorm2);
		hist.objval[k] = 0.5 * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso, DIAG>(xz, rho, rhoxz);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(x, xnorm);
		hist.eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(rhou, rhounorm);
		hist.eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if((hist.r_norm[k] < hist.eps_pri[k]) &&
		  (hist.s_norm[k] < hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gfptx_admmlasso.dat";
	std::string ukname = "uk_gfptx_admmlasso.dat";
	std::string zkname = "zk_gfptx_admmlasso.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << (float)hist.x_hist[j][i] << ",";
			resultfile1 << (float)hist.u_hist[j][i] << ",";
			resultfile2 << (float)hist.z_hist[j][i] << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

void ADMM_LASSO_POSIT8(posit8_t A[DIAG][DIAG],
					   posit8_t b[DIAG], posit8_t lambda,
					   posit8_t rho, posit8_t alpha){

	// parameters
	posit8_t oneminusalpha = p8_sub(convertDoubleToP8(1), alpha);
	posit8_t lambdadivrho = p8_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	posit8_t Atb[DIAG];
	posit8_t At[DIAG][DIAG];
	posit8_t AtA[DIAG][DIAG];
	posit8_t EYE[DIAG][DIAG];
	posit8_t rhoEYE[DIAG][DIAG];
	posit8_t AtAplusrhoeye[DIAG][DIAG];
	posit8_t L[DIAG][DIAG], U[DIAG][DIAG];
	posit8_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit8_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit8_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit8_t invLq[DIAG];
	posit8_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit8_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit8_t x_hatz[DIAG];

	// A'*b
	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	SoftPosit_Algebraobj.MAT_MUL<posit8_t,DIAG,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	SoftPosit_Algebraobj.IDENDTITY_MAT<posit8_t,DIAG,DIAG,convertDoubleToP8>(EYE);
	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit8_t,DIAG,DIAG,p8_mul>(
			EYE, rho, rhoEYE);
	SoftPosit_Algebraobj.MAT_ADD<posit8_t,DIAG,DIAG,p8_add>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit8_t,DIAG, convertDoubleToP8, p8_add, p8_sub, p8_mul, p8_div, p8_sqrt, p8_lt>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	SoftPosit_Algebraobj.MAT_QRINV<posit8_t,DIAG, convertDoubleToP8, p8_mul,
			p8_add, p8_sub, p8_sqrt, p8_div, p8_eq>(L, invL);
	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	SoftPosit_Algebraobj.ZEROS_VEC<posit8_t,DIAG,convertDoubleToP8>(x);
	SoftPosit_Algebraobj.ZEROS_VEC<posit8_t,DIAG,convertDoubleToP8>(z);
	SoftPosit_Algebraobj.ZEROS_VEC<posit8_t,DIAG,convertDoubleToP8>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<posit8_t> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,
				convertDoubleToP8,p8_mul,p8_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		SoftPosit_Algebraobj.VEC_EQ<posit8_t,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(
				x, alpha, alphax);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(
				zold, oneminusalpha, oneminusalphazold);
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		SoftPosit_Algebraobj.VEC_SCALAR_SUB<posit8_t,DIAG,p8_sub>(
				x_hatu, lambdadivrho, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_ADD<posit8_t,DIAG,p8_add>(
				x_hatu, lambdadivrho, x_hatu2);
		SoftPosit_Algebraobj.VEC_MINUS<posit8_t,DIAG, convertDoubleToP8, p8_mul>(
				x_hatu2, x_hatu2);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit8_t,DIAG,p8_lt>(
				x_hatu1, convertDoubleToP8(0), x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit8_t,DIAG,p8_lt>(
				x_hatu2, convertDoubleToP8(0), x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		posit8_t znorm;
		posit8_t Ax[DIAG], Axb[DIAG];
		posit8_t Axbnorm2;
		posit8_t xz[DIAG], rhoxz[DIAG];
		posit8_t xznorm, rhoxznorm;
		posit8_t xnorm;
		posit8_t rhou[DIAG];
		posit8_t rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t, DIAG, DIAG,convertDoubleToP8,p8_mul,p8_add>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add>(Axb, Axbnorm2);
		hist.objval[k] = p8_add(p8_mul(convertDoubleToP8(0.5), Axbnorm2), p8_mul(lambda, znorm));
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t, DIAG,p8_mul>(xz, rho, rhoxz);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(x, xnorm);
		hist.eps_pri[k] = p8_add(p8_mul(p8_sqrt(convertDoubleToP8(DIAG)),convertDoubleToP8(ABSTOL)),p8_mul(convertDoubleToP8(RELTOL),(p8_lt(znorm,xnorm)?xnorm:znorm)));
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t,DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(rhou, rhounorm);
		hist.eps_dual[k] = p8_add(p8_mul(p8_sqrt(convertDoubleToP8(DIAG)),convertDoubleToP8(ABSTOL)),p8_mul(convertDoubleToP8(RELTOL),rhounorm));
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if(p8_lt(hist.r_norm[k], hist.eps_pri[k]) &&
		  p8_lt(hist.s_norm[k], hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gposit8_t.dat";
	std::string ukname = "uk_gposit8_t.dat";
	std::string zkname = "zk_gposit8_t.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << convertP8ToDouble(hist.x_hist[j][i]) << ",";
			resultfile1 << convertP8ToDouble(hist.u_hist[j][i]) << ",";
			resultfile2 << convertP8ToDouble(hist.z_hist[j][i]) << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

void ADMM_LASSO_POSIT16(posit16_t A[DIAG][DIAG],
						posit16_t b[DIAG], posit16_t lambda,
						posit16_t rho, posit16_t alpha){

	// parameters
	posit16_t oneminusalpha = p16_sub(convertDoubleToP16(1), alpha);
	posit16_t lambdadivrho = p16_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	posit16_t Atb[DIAG];
	posit16_t At[DIAG][DIAG];
	posit16_t AtA[DIAG][DIAG];
	posit16_t EYE[DIAG][DIAG];
	posit16_t rhoEYE[DIAG][DIAG];
	posit16_t AtAplusrhoeye[DIAG][DIAG];
	posit16_t L[DIAG][DIAG], U[DIAG][DIAG];
	posit16_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit16_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit16_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit16_t invLq[DIAG];
	posit16_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit16_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit16_t x_hatz[DIAG];

	// A'*b
	SoftPosit_Algebraobj.MAT_TRANS<posit16_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit16_t,DIAG,DIAG,convertDoubleToP16,p16_mul,p16_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	SoftPosit_Algebraobj.MAT_MUL<posit16_t,DIAG,DIAG,DIAG,convertDoubleToP16,p16_mul,p16_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	SoftPosit_Algebraobj.IDENDTITY_MAT<posit16_t,DIAG,DIAG,convertDoubleToP16>(EYE);
	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit16_t,DIAG,DIAG,p16_mul>(
			EYE, rho, rhoEYE);
	SoftPosit_Algebraobj.MAT_ADD<posit16_t,DIAG,DIAG,p16_add>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit16_t,DIAG, convertDoubleToP16, p16_add, p16_sub, p16_mul, p16_div, p16_sqrt, p16_lt>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	SoftPosit_Algebraobj.MAT_QRINV<posit16_t,DIAG, convertDoubleToP16, p16_mul, p16_add,
			p16_sub, p16_sqrt, p16_div,p16_eq>(L, invL);
	SoftPosit_Algebraobj.MAT_TRANS<posit16_t,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	SoftPosit_Algebraobj.ZEROS_VEC<posit16_t,DIAG,convertDoubleToP16>(x);
	SoftPosit_Algebraobj.ZEROS_VEC<posit16_t,DIAG,convertDoubleToP16>(z);
	SoftPosit_Algebraobj.ZEROS_VEC<posit16_t,DIAG,convertDoubleToP16>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<posit16_t> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit16_t,DIAG,p16_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit16_t,DIAG,p16_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit16_t,DIAG,DIAG,
				convertDoubleToP16,p16_mul,p16_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit16_t,DIAG,DIAG,convertDoubleToP16,p16_mul,p16_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		SoftPosit_Algebraobj.VEC_EQ<posit16_t,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(
				x, alpha, alphax);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(
				zold, oneminusalpha, oneminusalphazold);
		SoftPosit_Algebraobj.VEC_ADD<posit16_t,DIAG,p16_add>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit16_t,DIAG,p16_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		SoftPosit_Algebraobj.VEC_SCALAR_SUB<posit16_t,DIAG,p16_sub>(
				x_hatu, lambdadivrho, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_ADD<posit16_t,DIAG,p16_add>(
				x_hatu, lambdadivrho, x_hatu2);
		SoftPosit_Algebraobj.VEC_MINUS<posit16_t,DIAG, convertDoubleToP16, p16_mul>(
				x_hatu2, x_hatu2);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit16_t,DIAG,p16_lt>(
				x_hatu1, convertDoubleToP16(0), x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit16_t,DIAG,p16_lt>(
				x_hatu2, convertDoubleToP16(0), x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit16_t,DIAG,p16_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit16_t,DIAG,p16_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit16_t,DIAG,p16_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		posit16_t znorm;
		posit16_t Ax[DIAG], Axb[DIAG];
		posit16_t Axbnorm2;
		posit16_t xz[DIAG], rhoxz[DIAG];
		posit16_t xznorm, rhoxznorm;
		posit16_t xnorm;
		posit16_t rhou[DIAG];
		posit16_t rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit16_t, DIAG,convertDoubleToP16,p16_mul,p16_add,p16_sqrt>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit16_t, DIAG, DIAG,convertDoubleToP16,p16_mul,p16_add>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit16_t, DIAG,p16_sub>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit16_t, DIAG,convertDoubleToP16,p16_mul,p16_add>(Axb, Axbnorm2);
		hist.objval[k] = p16_add(p16_mul(convertDoubleToP16(0.5), Axbnorm2), p16_mul(lambda, znorm));
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit16_t, DIAG,p16_sub>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit16_t, DIAG,convertDoubleToP16,p16_mul,p16_add,p16_sqrt>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit16_t, DIAG,p16_mul>(xz, rho, rhoxz);
		SoftPosit_Algebraobj.VEC_NORM<posit16_t, DIAG,convertDoubleToP16,p16_mul,p16_add,p16_sqrt>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit16_t, DIAG,convertDoubleToP16,p16_mul,p16_add,p16_sqrt>(x, xnorm);
		hist.eps_pri[k] = p16_add(p16_mul(p16_sqrt(convertDoubleToP16(DIAG)),convertDoubleToP16(ABSTOL)),p16_mul(convertDoubleToP16(RELTOL),(p16_lt(znorm,xnorm)?xnorm:znorm)));
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit16_t,DIAG,convertDoubleToP16,p16_mul,p16_add,p16_sqrt>(rhou, rhounorm);
		hist.eps_dual[k] = p16_add(p16_mul(p16_sqrt(convertDoubleToP16(DIAG)),convertDoubleToP16(ABSTOL)),p16_mul(convertDoubleToP16(RELTOL),rhounorm));
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if(p16_lt(hist.r_norm[k], hist.eps_pri[k]) &&
		  p16_lt(hist.s_norm[k], hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gposit16_t.dat";
	std::string ukname = "uk_gposit16_t.dat";
	std::string zkname = "zk_gposit16_t.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << convertP16ToDouble(hist.x_hist[j][i]) << ",";
			resultfile1 << convertP16ToDouble(hist.u_hist[j][i]) << ",";
			resultfile2 << convertP16ToDouble(hist.z_hist[j][i]) << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

void ADMM_LASSO_POSIT32(posit32_t A[DIAG][DIAG],
						posit32_t b[DIAG], posit32_t lambda,
						posit32_t rho, posit32_t alpha){

	// parameters
	posit32_t oneminusalpha = p32_sub(convertDoubleToP32(1), alpha);
	posit32_t lambdadivrho = p32_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	posit32_t Atb[DIAG];
	posit32_t At[DIAG][DIAG];
	posit32_t AtA[DIAG][DIAG];
	posit32_t EYE[DIAG][DIAG];
	posit32_t rhoEYE[DIAG][DIAG];
	posit32_t AtAplusrhoeye[DIAG][DIAG];
	posit32_t L[DIAG][DIAG], U[DIAG][DIAG];
	posit32_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit32_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit32_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit32_t invLq[DIAG];
	posit32_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit32_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit32_t x_hatz[DIAG];

	// A'*b
	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	printmatrix(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	printvector(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	SoftPosit_Algebraobj.MAT_MUL<posit32_t,DIAG,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	printmatrix(AtA);
#endif
	SoftPosit_Algebraobj.IDENDTITY_MAT<posit32_t,DIAG,DIAG,convertDoubleToP32>(EYE);
	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit32_t,DIAG,DIAG,p32_mul>(
			EYE, rho, rhoEYE);
	SoftPosit_Algebraobj.MAT_ADD<posit32_t,DIAG,DIAG,p32_add>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	printmatrix(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit32_t,DIAG, convertDoubleToP32, p32_add, p32_sub, p32_mul, p32_div, p32_sqrt, p32_lt>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	printmatrix(L);
	std::cout << "U:" << std::endl;
	printmatrix(U);
#endif
	// invers L and U;
	SoftPosit_Algebraobj.MAT_QRINV<posit32_t,DIAG, convertDoubleToP32, p32_mul,
			p32_add, p32_sub, p32_sqrt, p32_div, p32_eq>(L, invL);
	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	SoftPosit_Algebraobj.ZEROS_VEC<posit32_t,DIAG,convertDoubleToP32>(x);
	SoftPosit_Algebraobj.ZEROS_VEC<posit32_t,DIAG,convertDoubleToP32>(z);
	SoftPosit_Algebraobj.ZEROS_VEC<posit32_t,DIAG,convertDoubleToP32>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<posit32_t> hist; // iteration record and early termination
	for(int k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		printvector(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,
				convertDoubleToP32,p32_mul,p32_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		printmatrix(invL);
		std::cout << "invU:" << std::endl;
		printmatrix(invU);
		std::cout << "invLq:" << std::endl;
		printvector(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		printvector(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		SoftPosit_Algebraobj.VEC_EQ<posit32_t,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(
				x, alpha, alphax);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(
				zold, oneminusalpha, oneminusalphazold);
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		printvector(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		printvector(x_hatu);
#endif
		SoftPosit_Algebraobj.VEC_SCALAR_SUB<posit32_t,DIAG,p32_sub>(
				x_hatu, lambdadivrho, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_ADD<posit32_t,DIAG,p32_add>(
				x_hatu, lambdadivrho, x_hatu2);
		SoftPosit_Algebraobj.VEC_MINUS<posit32_t,DIAG, convertDoubleToP32, p32_mul>(
				x_hatu2, x_hatu2);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit32_t,DIAG,p32_lt>(
				x_hatu1, convertDoubleToP32(0), x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit32_t,DIAG,p32_lt>(
				x_hatu2, convertDoubleToP32(0), x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		printvector(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		printvector(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		printvector(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		printvector(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		printvector(u);
#endif
#if defined(RECORD_RESULT)
		posit32_t znorm;
		posit32_t Ax[DIAG], Axb[DIAG];
		posit32_t Axbnorm2;
		posit32_t xz[DIAG], rhoxz[DIAG];
		posit32_t xznorm, rhoxznorm;
		posit32_t xnorm;
		posit32_t rhou[DIAG];
		posit32_t rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t, DIAG, DIAG,convertDoubleToP32,p32_mul,p32_add>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add>(Axb, Axbnorm2);
		hist.objval[k] = p32_add(p32_mul(convertDoubleToP32(0.5), Axbnorm2), p32_mul(lambda, znorm));
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t, DIAG,p32_mul>(xz, rho, rhoxz);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(rhoxz, rhoxznorm);
		hist.s_norm[k] = rhoxznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(x, xnorm);
		hist.eps_pri[k] = p32_add(p32_mul(p32_sqrt(convertDoubleToP32(DIAG)),convertDoubleToP32(ABSTOL)),p32_mul(convertDoubleToP32(RELTOL),(p32_lt(znorm,xnorm)?xnorm:znorm)));
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t,DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(rhou, rhounorm);
		hist.eps_dual[k] = p32_add(p32_mul(p32_sqrt(convertDoubleToP32(DIAG)),convertDoubleToP32(ABSTOL)),p32_mul(convertDoubleToP32(RELTOL),rhounorm));
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
#if defined(EARLY_TERMINATE)
		if(p32_lt(hist.r_norm[k], hist.eps_pri[k]) &&
		  p32_lt(hist.s_norm[k], hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	printvector(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk_gposit32_t.dat";
	std::string ukname = "uk_gposit32_t.dat";
	std::string zkname = "zk_gposit32_t.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	for(int i=0; i<MAX_ITER; i++){
		for(int j=0;j<DIAG;j++){
			resultfile << convertP32ToDouble(hist.x_hist[j][i]) << ",";
			resultfile1 << convertP32ToDouble(hist.u_hist[j][i]) << ",";
			resultfile2 << convertP32ToDouble(hist.z_hist[j][i]) << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
#endif
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to finish the iteration"
			  << std::endl << std::endl;
#endif// endif TIME_PROFILE

}

int main(int argc, char** argv){

	std::string clockname = "timeprofile.txt";
	std::string Amatrixname = "A.csv";
	std::string bvectorname = "b.csv";
	std::string lambdaname = "lambda.csv";
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
	double A_c[DIAG][DIAG], b_c[DIAG];
	double lambda_c, rho_c, alpha_c;
	memcpy(A_c, AmatrixT.data(), sizeof(double)*DIAG*DIAG);
	memcpy(b_c, bvector.data(), sizeof(double)*DIAG);
	lambda_c = lambda;
	rho_c = 1;
	alpha_c = 1;
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	ADMM_LASSO_DOUBLE(A_c, b_c, lambda_c, rho_c, alpha_c);
#endif
#if defined(GENERAL_FLOAT_PRECISION)
	float A[DIAG][DIAG], b[DIAG], lambda, rho, alpha;
	ADMM_LASSO_FLOAT(A, b, lambda, rho, alpha);
#endif
#if defined(COMSTOM_FLOAT_PRECISION)
	fptx_admmlasso A[DIAG][DIAG], b[DIAG], lambda, rho, alpha;
	ADMM_LASSO_FXPT(A, b, lambda, rho, alpha);
#endif
#if defined(XILINX_FIXED_PRECISION)
#endif
#if defined(SOFT_POSIT_PRECISION)
#endif

	return 0;
}

