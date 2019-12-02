/*
 * admmlasso_gdouble.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_gdouble.hpp"

namespace plt = matplotlibcpp;

void ADMM_LASSO_DOUBLE(double A[DIAG][DIAG], double At[DIAG][DIAG],
					   double invL[DIAG][DIAG],double invU[DIAG][DIAG],
					   double b[DIAG], double lambda,
					   double rho, double alpha){

	// parameters
	double oneminusalpha = 1 - alpha;
	double lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
#if defined(DEBUG_ITER)
	DISPLAY DISPLAYobj;
#endif

	double Atb[DIAG];
//	double At[DIAG][DIAG];
//	double AtA[DIAG][DIAG];
//	double EYE[DIAG][DIAG];
//	double rhoEYE[DIAG][DIAG];
//	double AtAplusrhoeye[DIAG][DIAG];
//	double L[DIAG][DIAG], U[DIAG][DIAG];
//	double invL[DIAG][DIAG], invU[DIAG][DIAG];
	double x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	double zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	double invLq[DIAG];
	double alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	double x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	double x_hatz[DIAG];
#if defined(DEBUG_ITER)
	std::cout << "A:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(A);
	std::cout << "U_inv:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(invU);
	std::cout << "L_inv:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(invL);
#endif

	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<double,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<double,DIAG,DIAG>(At, A, AtA);
//#if defined(DEBUG_ITER)
//	std::cout << "AtA:" << std::endl;
//	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtA);
//#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<double,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<double,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<double,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
//#if defined(DEBUG_ITER)
//	std::cout << "AtAplusrhoeye:" << std::endl;
//	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtAplusrhoeye);
//#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<double,DIAG>(
//			AtAplusrhoeye, L, U);
//#if defined(DEBUG_ITER)
//	std::cout << "L:" << std::endl;
//	DISPLAYobj.printmatrix<double,DIAG,DIAG>(L);
//	std::cout << "U:" << std::endl;
//	DISPLAYobj.printmatrix<double,DIAG,DIAG>(U);
//#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<double,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	struct histry hist; // iteration record and early termination
    int k = 0;
	for(k=0;k<MAX_ITER;k++){
//	for(int k=0;k<1;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x);
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
		DISPLAYobj.printvector<double,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<double,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(u);
#endif
#if defined(RECORD_RESULT)
		double znorm;
		double Ax[DIAG], Axb[DIAG];
		double Axbnorm2;
		double xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		double xznorm, rhozoldznorm;
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
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			hist.r_norm[k], hist.eps_pri[k],
			hist.s_norm[k], hist.eps_dual[k], hist.objval[k]);
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
	DISPLAYobj.printvector<double,DIAG>(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk.dat";
	std::string ukname = "uk.dat";
	std::string zkname = "zk.dat";
	std::string objvalname = "objvalk.dat";
	std::string r_normname = "r_normk.dat";
	std::string s_normname = "s_normk.dat";
	std::string eps_priname = "eps_prik.dat";
	std::string eps_dualname = "eps_dualk.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	std::ofstream resultfile3(objvalname);
	std::ofstream resultfile4(r_normname);
	std::ofstream resultfile5(s_normname);
	std::ofstream resultfile6(eps_priname);
	std::ofstream resultfile7(eps_dualname);
	for(int i=0; i<k; i++){
		resultfile3 << hist.objval[i] << "\n";
		resultfile4 << hist.r_norm[i] << "\n";
		resultfile5 << hist.s_norm[i] << "\n";
		resultfile6 << hist.eps_pri[i] << "\n";
		resultfile7 << hist.eps_dual[i] << "\n";
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
	resultfile3.close();
	resultfile4.close();
	resultfile5.close();
	resultfile6.close();
	resultfile7.close();
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


#ifdef PLOT_FIGURE
	std::string figurename1 = "ADMM_Obj_gfloat.png";
	// Matplotlib plotting
	std::vector<double> xplot, yplot;
	for(int i=0;i<k;i++){
		xplot.push_back(i+1);
		yplot.push_back(hist.objval[i]);
	}
	plt::named_plot( "Objective function", xplot, yplot);
	plt::title("ADMM Objective");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Objective Value");
	plt::legend();
	plt::save(figurename1);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
	std::string figurename2 = "ADMM_NORM_gfloat.png";
	// Matplotlib plotting
	std::vector<double> yplot1, yplot2, yplot3, yplot4;
	for(int i=0;i<k;i++){
		yplot1.push_back(hist.r_norm[i]);
		yplot2.push_back(hist.eps_pri[i]);
		yplot3.push_back(hist.s_norm[i]);
		yplot4.push_back(hist.eps_dual[i]);
	}
	plt::title("ADMM NORM");
	plt::subplot(2,1,1);
	plt::named_semilogy("r NORM", xplot, yplot1, "k");
	plt::named_semilogy("eps pri", xplot, yplot2, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||r||_2");
	plt::legend();
	plt::subplot(2,1,2);
	plt::named_semilogy("s NORM", xplot, yplot3, "k");
	plt::named_semilogy("eps dual", xplot, yplot4, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||s||_2");
	plt::legend();
	plt::save(figurename2);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
#endif// endif PLOT_FIGURE
}


void ADMM_LASSO_DOUBLE(double **A, double **At,
					   double **invL,double **invU,
					   double *b, double lambda,
					   double rho, double alpha){

//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// parameters
	double oneminusalpha = 1 - alpha;
	double lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
#if defined(DEBUG_ITER)
	DISPLAY DISPLAYobj;
#endif
//	double **At;
//	At = (double**) malloc(sizeof(double*)*DIAG);
//	for(int i=0;i<DIAG;i++)
//		At[i] = (double*) malloc(sizeof(double)*DIAG);
//	double **AtA;
//	double **EYE;
//	double **rhoEYE;
//	double **AtAplusrhoeye;
//	double **L, **U;
//	double **invL, **invU;
	double *Atb;
	double *x, *z, *u;
	double  *zold;
	Atb = (double*) malloc(sizeof(double)*DIAG);
	x = (double*) malloc(sizeof(double)*DIAG);
	z = (double*) malloc(sizeof(double)*DIAG);
	u = (double*) malloc(sizeof(double)*DIAG);
	zold = (double*) malloc(sizeof(double)*DIAG);
//	At = (double**) malloc(sizeof(double*)*DIAG);
//	AtA = (double**) malloc(sizeof(double*)*DIAG);
//	EYE = (double**) malloc(sizeof(double*)*DIAG);
//	rhoEYE = (double**) malloc(sizeof(double*)*DIAG);
//	AtAplusrhoeye = (double**) malloc(sizeof(double*)*DIAG);
//	L = (double**) malloc(sizeof(double*)*DIAG);
//	U = (double**) malloc(sizeof(double*)*DIAG);
//	invL = (double**) malloc(sizeof(double*)*DIAG);
//	invU = (double**) malloc(sizeof(double*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (double*) malloc(sizeof(double)*DIAG);
//		AtA[i] = (double*) malloc(sizeof(double)*DIAG);
//		EYE[i] = (double*) malloc(sizeof(double)*DIAG);
//		rhoEYE[i] = (double*) malloc(sizeof(double)*DIAG);
//		AtAplusrhoeye[i] = (double*) malloc(sizeof(double)*DIAG);
//		L[i] = (double*) malloc(sizeof(double)*DIAG);
//		U[i] = (double*) malloc(sizeof(double)*DIAG);
//		invL[i] = (double*) malloc(sizeof(double)*DIAG);
//		invU[i] = (double*) malloc(sizeof(double)*DIAG);
//	}

//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<double,DIAG>(Atb);
#endif
//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<double,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtA);
#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<double,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<double,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<double,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtAplusrhoeye);
#endif
//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<double,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<double,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(invL, invU);
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
	// iteration record and early termination
	double** u_hist = NULL;
	double** x_hist = NULL;
	double** z_hist = NULL;
	double* objval = NULL;
	double* r_norm = NULL;
	double* s_norm = NULL;
	double* eps_pri = NULL;
	double* eps_dual = NULL;
	u_hist = (double**) malloc(sizeof(double*)*DIAG);
	x_hist = (double**) malloc(sizeof(double*)*DIAG);
	z_hist = (double**) malloc(sizeof(double*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (double*) malloc(sizeof(double)*MAX_ITER);
		x_hist[i] = (double*) malloc(sizeof(double)*MAX_ITER);
		z_hist[i] = (double*) malloc(sizeof(double)*MAX_ITER);
	}
	objval = (double*) malloc(sizeof(double)*MAX_ITER);
	r_norm = (double*) malloc(sizeof(double)*MAX_ITER);
	s_norm = (double*) malloc(sizeof(double)*MAX_ITER);
	eps_pri = (double*) malloc(sizeof(double)*MAX_ITER);
	eps_dual = (double*) malloc(sizeof(double)*MAX_ITER);

	int k = 0;
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	for(k=0;k<MAX_ITER;k++){
		double *zminusu=NULL, *rhozminusu=NULL, *q=NULL;
		double *invLq=NULL;
		double *alphax=NULL,*oneminusalphazold=NULL,*x_hat=NULL;
		double *x_hatu=NULL,*x_hatu1=NULL,*x_hatu2=NULL;
		double *x_hatz=NULL;
		zminusu = (double*) malloc(sizeof(double)*DIAG);
		rhozminusu = (double*) malloc(sizeof(double)*DIAG);
		q = (double*) malloc(sizeof(double)*DIAG);
		invLq = (double*) malloc(sizeof(double)*DIAG);
		alphax = (double*) malloc(sizeof(double)*DIAG);
		oneminusalphazold = (double*) malloc(sizeof(double)*DIAG);
		x_hat = (double*) malloc(sizeof(double)*DIAG);
		x_hatu = (double*) malloc(sizeof(double)*DIAG);
		x_hatu1 = (double*) malloc(sizeof(double)*DIAG);
		x_hatu2 = (double*) malloc(sizeof(double)*DIAG);
		x_hatz = (double*) malloc(sizeof(double)*DIAG);
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(q);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
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
		DISPLAYobj.printvector<double,DIAG>(x_hat);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<double,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(z);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(u);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
#if defined(RECORD_RESULT)
		double znorm;
		double Axbnorm2;
		double xznorm, rhozoldznorm;
		double xnorm;
		double rhounorm;
		double *Ax=NULL, *Axb=NULL;
		double *xz=NULL, *zoldz=NULL, *rhozoldz=NULL;
		double *rhou=NULL;
		Ax = (double*) realloc(Ax,sizeof(double)*DIAG);
		Axb = (double*) realloc(Axb,sizeof(double)*DIAG);
		xz = (double*) realloc(xz,sizeof(double)*DIAG);
		zoldz = (double*) realloc(zoldz,sizeof(double)*DIAG);
		rhozoldz = (double*) realloc(rhozoldz,sizeof(double)*DIAG);
		rhou = (double*) realloc(rhou,sizeof(double)*DIAG);

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<double, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<double, DIAG>(Axb, Axbnorm2);
		objval[k] = 0.5 * Axbnorm2 + lambda * znorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(xz, xznorm);
		r_norm[k] = xznorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(x, xnorm);
		eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhou, rhounorm);
		eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			r_norm[k], eps_pri[k],
			s_norm[k], eps_dual[k], objval[k]);

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

//		free(Ax);
//		free(Axb);
//		free(xz);
//		free(zoldz);
//		free(rhozoldz);
//		free(rhou);
#if defined(EARLY_TERMINATE)
		if((r_norm[k] < eps_pri[k]) &&
		  (s_norm[k] < eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<double,DIAG>(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk.dat";
	std::string ukname = "uk.dat";
	std::string zkname = "zk.dat";
	std::string objvalname = "objvalk.dat";
	std::string r_normname = "r_normk.dat";
	std::string s_normname = "s_normk.dat";
	std::string eps_priname = "eps_prik.dat";
	std::string eps_dualname = "eps_dualk.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	std::ofstream resultfile3(objvalname);
	std::ofstream resultfile4(r_normname);
	std::ofstream resultfile5(s_normname);
	std::ofstream resultfile6(eps_priname);
	std::ofstream resultfile7(eps_dualname);
	for(int i=0; i<k; i++){
		resultfile3 << objval[i] << "\n";
		resultfile4 << r_norm[i] << "\n";
		resultfile5 << s_norm[i] << "\n";
		resultfile6 << eps_pri[i] << "\n";
		resultfile7 << eps_dual[i] << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << x_hist[j][i] << ",";
			resultfile1 << u_hist[j][i] << ",";
			resultfile2 << z_hist[j][i] << ",";
		}
		resultfile << "\n";
		resultfile1 << "\n";
		resultfile2 << "\n";
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();
	resultfile3.close();
	resultfile4.close();
	resultfile5.close();
	resultfile6.close();
	resultfile7.close();
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


#ifdef PLOT_FIGURE
	std::string figurename1 = "ADMM_Obj_gfloat.png";
	// Matplotlib plotting
	std::vector<double> xplot, yplot;
	for(int i=0;i<k;i++){
		xplot.push_back(i+1);
		yplot.push_back(objval[i]);
	}
	plt::named_plot( "Objective function", xplot, yplot);
	plt::title("ADMM Objective");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Objective Value");
	plt::legend();
	plt::save(figurename1);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
	std::string figurename2 = "ADMM_NORM_gfloat.png";
	// Matplotlib plotting
	std::vector<double> yplot1, yplot2, yplot3, yplot4;
	for(int i=0;i<k;i++){
		yplot1.push_back(r_norm[i]);
		yplot2.push_back(eps_pri[i]);
		yplot3.push_back(s_norm[i]);
		yplot4.push_back(eps_dual[i]);
	}
	plt::title("ADMM NORM");
	plt::subplot(2,1,1);
	plt::named_semilogy("r NORM", xplot, yplot1, "k");
	plt::named_semilogy("eps pri", xplot, yplot2, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||r||_2");
	plt::legend();
	plt::subplot(2,1,2);
	plt::named_semilogy("s NORM", xplot, yplot3, "k");
	plt::named_semilogy("eps dual", xplot, yplot4, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||s||_2");
	plt::legend();
	plt::save(figurename2);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
#endif// endif PLOT_FIGURE
}



void ADMM_LASSO_DOUBLE(double **A, double *b, double lambda,
					   double rho, double alpha){

	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// parameters
	double oneminusalpha = 1 - alpha;
	double lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
#if defined(DEBUG_ITER)
	DISPLAY DISPLAYobj;
#endif
	double *Atb;
	double **At;
	double **AtA;
	double **EYE;
	double **rhoEYE;
	double **AtAplusrhoeye;
	double **L, **U;
	double **invL, **invU;
	double *x, *zold, *z, *u;
	double *zminusu, *rhozminusu, *q;
	double *invLq;
	double *alphax,*oneminusalphazold,*x_hat;
	double *x_hatu,*x_hatu1,*x_hatu2;
	double *x_hatz;
	Atb = (double*) malloc(sizeof(double)*DIAG);
	x = (double*) malloc(sizeof(double)*DIAG);
	zold = (double*) malloc(sizeof(double)*DIAG);
	z = (double*) malloc(sizeof(double)*DIAG);
	u = (double*) malloc(sizeof(double)*DIAG);
	zminusu = (double*) malloc(sizeof(double)*DIAG);
	rhozminusu = (double*) malloc(sizeof(double)*DIAG);
	q = (double*) malloc(sizeof(double)*DIAG);
	invLq = (double*) malloc(sizeof(double)*DIAG);
	alphax = (double*) malloc(sizeof(double)*DIAG);
	oneminusalphazold = (double*) malloc(sizeof(double)*DIAG);
	x_hat = (double*) malloc(sizeof(double)*DIAG);
	x_hatu = (double*) malloc(sizeof(double)*DIAG);
	x_hatu1 = (double*) malloc(sizeof(double)*DIAG);
	x_hatu2 = (double*) malloc(sizeof(double)*DIAG);
	x_hatz = (double*) malloc(sizeof(double)*DIAG);
	At = (double**) malloc(sizeof(double*)*DIAG);
	AtA = (double**) malloc(sizeof(double*)*DIAG);
	EYE = (double**) malloc(sizeof(double*)*DIAG);
	rhoEYE = (double**) malloc(sizeof(double*)*DIAG);
	AtAplusrhoeye = (double**) malloc(sizeof(double*)*DIAG);
	L = (double**) malloc(sizeof(double*)*DIAG);
	U = (double**) malloc(sizeof(double*)*DIAG);
	invL = (double**) malloc(sizeof(double*)*DIAG);
	invU = (double**) malloc(sizeof(double*)*DIAG);
	for(int i=0;i<DIAG;i++){
		At[i] = (double*) malloc(sizeof(double)*DIAG);
		AtA[i] = (double*) malloc(sizeof(double)*DIAG);
		EYE[i] = (double*) malloc(sizeof(double)*DIAG);
		rhoEYE[i] = (double*) malloc(sizeof(double)*DIAG);
		AtAplusrhoeye[i] = (double*) malloc(sizeof(double)*DIAG);
		L[i] = (double*) malloc(sizeof(double)*DIAG);
		U[i] = (double*) malloc(sizeof(double)*DIAG);
		invL[i] = (double*) malloc(sizeof(double)*DIAG);
		invU[i] = (double*) malloc(sizeof(double)*DIAG);
	}

	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*b
	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<double,DIAG>(Atb);
#endif
	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Float_Point_Algebraobj.MAT_MUL<double,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtA);
#endif
	Float_Point_Algebraobj.IDENDTITY_MAT<double,DIAG,DIAG>(EYE);
	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<double,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Float_Point_Algebraobj.MAT_ADD<double,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtAplusrhoeye);
#endif
	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Float_Point_Algebraobj.LU_CHOLBANACHROUT<double,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(U);
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
	int k = 0;
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(q);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
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
		DISPLAYobj.printvector<double,DIAG>(x_hat);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<double,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(z);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(u);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
#if defined(RECORD_RESULT)
		double znorm;
		double *Ax, *Axb;
		double Axbnorm2;
		double *xz, *zoldz, *rhozoldz;
		double xznorm, rhozoldznorm;
		double xnorm;
		double rhou[DIAG];
		double rhounorm;
		Ax = (double*) malloc(sizeof(double)*DIAG);
		Axb = (double*) malloc(sizeof(double)*DIAG);
		xz = (double*) malloc(sizeof(double)*DIAG);
		zoldz = (double*) malloc(sizeof(double)*DIAG);
		rhozoldz = (double*) malloc(sizeof(double)*DIAG);

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
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			hist.r_norm[k], hist.eps_pri[k],
			hist.s_norm[k], hist.eps_dual[k], hist.objval[k]);
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
	DISPLAYobj.printvector<double,DIAG>(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk.dat";
	std::string ukname = "uk.dat";
	std::string zkname = "zk.dat";
	std::string objvalname = "objvalk.dat";
	std::string r_normname = "r_normk.dat";
	std::string s_normname = "s_normk.dat";
	std::string eps_priname = "eps_prik.dat";
	std::string eps_dualname = "eps_dualk.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	std::ofstream resultfile3(objvalname);
	std::ofstream resultfile4(r_normname);
	std::ofstream resultfile5(s_normname);
	std::ofstream resultfile6(eps_priname);
	std::ofstream resultfile7(eps_dualname);
	for(int i=0; i<k; i++){
		resultfile3 << hist.objval[i] << "\n";
		resultfile4 << hist.r_norm[i] << "\n";
		resultfile5 << hist.s_norm[i] << "\n";
		resultfile6 << hist.eps_pri[i] << "\n";
		resultfile7 << hist.eps_dual[i] << "\n";
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
	resultfile3.close();
	resultfile4.close();
	resultfile5.close();
	resultfile6.close();
	resultfile7.close();
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

#ifdef PLOT_FIGURE
	std::string figurename1 = "ADMM_Obj_gfloat.png";
	// Matplotlib plotting
	std::vector<double> xplot, yplot;
	for(int i=0;i<k;i++){
		xplot.push_back(i+1);
		yplot.push_back(hist.objval[i]);
	}
	plt::named_plot( "Objective function", xplot, yplot);
	plt::title("ADMM Objective");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Objective Value");
	plt::legend();
	plt::save(figurename1);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
	std::string figurename2 = "ADMM_NORM_gfloat.png";
	// Matplotlib plotting
	std::vector<double> yplot1, yplot2, yplot3, yplot4;
	for(int i=0;i<k;i++){
		yplot1.push_back(hist.r_norm[i]);
		yplot2.push_back(hist.eps_pri[i]);
		yplot3.push_back(hist.s_norm[i]);
		yplot4.push_back(hist.eps_dual[i]);
	}
	plt::title("ADMM NORM");
	plt::subplot(2,1,1);
	plt::named_semilogy("r NORM", xplot, yplot1, "k");
	plt::named_semilogy("eps pri", xplot, yplot2, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||r||_2");
	plt::legend();
	plt::subplot(2,1,2);
	plt::named_semilogy("s NORM", xplot, yplot3, "k");
	plt::named_semilogy("eps dual", xplot, yplot4, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||s||_2");
	plt::legend();
	plt::save(figurename2);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
#endif// endif PLOT_FIGURE
}


void ADMM_LASSO_DOUBLE(double **A, double *b, double lambda){


	double rho = 1;
	double alpha =1;
	std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// parameters
	double oneminusalpha = 1 - alpha;
	double lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
#if defined(DEBUG_ITER)
	DISPLAY DISPLAYobj;
#endif
	double *Atb;
	double **At;
	double **AtA;
	double **EYE;
	double **rhoEYE;
	double **AtAplusrhoeye;
	double **L, **U;
	double **invL, **invU;
	double *x, *zold, *z, *u;
	double *zminusu, *rhozminusu, *q;
	double *invLq;
	double *alphax,*oneminusalphazold,*x_hat;
	double *x_hatu,*x_hatu1,*x_hatu2;
	double *x_hatz;
	Atb = (double*) malloc(sizeof(double)*DIAG);
	x = (double*) malloc(sizeof(double)*DIAG);
	zold = (double*) malloc(sizeof(double)*DIAG);
	z = (double*) malloc(sizeof(double)*DIAG);
	u = (double*) malloc(sizeof(double)*DIAG);
	zminusu = (double*) malloc(sizeof(double)*DIAG);
	rhozminusu = (double*) malloc(sizeof(double)*DIAG);
	q = (double*) malloc(sizeof(double)*DIAG);
	invLq = (double*) malloc(sizeof(double)*DIAG);
	alphax = (double*) malloc(sizeof(double)*DIAG);
	oneminusalphazold = (double*) malloc(sizeof(double)*DIAG);
	x_hat = (double*) malloc(sizeof(double)*DIAG);
	x_hatu = (double*) malloc(sizeof(double)*DIAG);
	x_hatu1 = (double*) malloc(sizeof(double)*DIAG);
	x_hatu2 = (double*) malloc(sizeof(double)*DIAG);
	x_hatz = (double*) malloc(sizeof(double)*DIAG);
	At = (double**) malloc(sizeof(double*)*DIAG);
	AtA = (double**) malloc(sizeof(double*)*DIAG);
	EYE = (double**) malloc(sizeof(double*)*DIAG);
	rhoEYE = (double**) malloc(sizeof(double*)*DIAG);
	AtAplusrhoeye = (double**) malloc(sizeof(double*)*DIAG);
	L = (double**) malloc(sizeof(double*)*DIAG);
	U = (double**) malloc(sizeof(double*)*DIAG);
	invL = (double**) malloc(sizeof(double*)*DIAG);
	invU = (double**) malloc(sizeof(double*)*DIAG);
	for(int i=0;i<DIAG;i++){
		At[i] = (double*) malloc(sizeof(double)*DIAG);
		AtA[i] = (double*) malloc(sizeof(double)*DIAG);
		EYE[i] = (double*) malloc(sizeof(double)*DIAG);
		rhoEYE[i] = (double*) malloc(sizeof(double)*DIAG);
		AtAplusrhoeye[i] = (double*) malloc(sizeof(double)*DIAG);
		L[i] = (double*) malloc(sizeof(double)*DIAG);
		U[i] = (double*) malloc(sizeof(double)*DIAG);
		invL[i] = (double*) malloc(sizeof(double)*DIAG);
		invU[i] = (double*) malloc(sizeof(double)*DIAG);
	}

	std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*b
	Float_Point_Algebraobj.MAT_TRANS<double,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<double,DIAG>(Atb);
#endif
	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
	Float_Point_Algebraobj.MAT_MUL<double,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtA);
#endif
	Float_Point_Algebraobj.IDENDTITY_MAT<double,DIAG,DIAG>(EYE);
	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<double,DIAG,DIAG>(
			EYE, rho, rhoEYE);
	Float_Point_Algebraobj.MAT_ADD<double,DIAG,DIAG>(AtA,
			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(AtAplusrhoeye);
#endif
	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
	Float_Point_Algebraobj.LU_CHOLBANACHROUT<double,DIAG>(
			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<double,DIAG,DIAG>(U);
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
	int k = 0;
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(q);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<double,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<double,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
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
		DISPLAYobj.printvector<double,DIAG>(x_hat);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<double,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(z);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<double,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<double,DIAG>(u);
#endif
		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
#if defined(RECORD_RESULT)
		double znorm;
		double *Ax, *Axb;
		double Axbnorm2;
		double *xz, *zoldz, *rhozoldz;
		double xznorm, rhozoldznorm;
		double xnorm;
		double rhou[DIAG];
		double rhounorm;
		Ax = (double*) malloc(sizeof(double)*DIAG);
		Axb = (double*) malloc(sizeof(double)*DIAG);
		xz = (double*) malloc(sizeof(double)*DIAG);
		zoldz = (double*) malloc(sizeof(double)*DIAG);
		rhozoldz = (double*) malloc(sizeof(double)*DIAG);

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
		Float_Point_Algebraobj.VEC_SUB<double, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<double, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<double, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			hist.r_norm[k], hist.eps_pri[k],
			hist.s_norm[k], hist.eps_dual[k], hist.objval[k]);
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
	DISPLAYobj.printvector<double,DIAG>(x);
#endif
#if defined(RECORD_RESULT)
	std::string xkname = "xk.dat";
	std::string ukname = "uk.dat";
	std::string zkname = "zk.dat";
	std::string objvalname = "objvalk.dat";
	std::string r_normname = "r_normk.dat";
	std::string s_normname = "s_normk.dat";
	std::string eps_priname = "eps_prik.dat";
	std::string eps_dualname = "eps_dualk.dat";
	std::ofstream resultfile(xkname);
	std::ofstream resultfile1(ukname);
	std::ofstream resultfile2(zkname);
	std::ofstream resultfile3(objvalname);
	std::ofstream resultfile4(r_normname);
	std::ofstream resultfile5(s_normname);
	std::ofstream resultfile6(eps_priname);
	std::ofstream resultfile7(eps_dualname);
	for(int i=0; i<k; i++){
		resultfile3 << hist.objval[i] << "\n";
		resultfile4 << hist.r_norm[i] << "\n";
		resultfile5 << hist.s_norm[i] << "\n";
		resultfile6 << hist.eps_pri[i] << "\n";
		resultfile7 << hist.eps_dual[i] << "\n";
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
	resultfile3.close();
	resultfile4.close();
	resultfile5.close();
	resultfile6.close();
	resultfile7.close();
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


#ifdef PLOT_FIGURE
	std::string figurename1 = "ADMM_Obj_gfloat.png";
	// Matplotlib plotting
	std::vector<double> xplot, yplot;
	for(int i=0;i<k;i++){
		xplot.push_back(i+1);
		yplot.push_back(hist.objval[i]);
	}
	plt::named_plot( "Objective function", xplot, yplot);
	plt::title("ADMM Objective");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Objective Value");
	plt::legend();
	plt::save(figurename1);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
	std::string figurename2 = "ADMM_NORM_gfloat.png";
	// Matplotlib plotting
	std::vector<double> yplot1, yplot2, yplot3, yplot4;
	for(int i=0;i<k;i++){
		yplot1.push_back(hist.r_norm[i]);
		yplot2.push_back(hist.eps_pri[i]);
		yplot3.push_back(hist.s_norm[i]);
		yplot4.push_back(hist.eps_dual[i]);
	}
	plt::title("ADMM NORM");
	plt::subplot(2,1,1);
	plt::named_semilogy("r NORM", xplot, yplot1, "k");
	plt::named_semilogy("eps pri", xplot, yplot2, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||r||_2");
	plt::legend();
	plt::subplot(2,1,2);
	plt::named_semilogy("s NORM", xplot, yplot3, "k");
	plt::named_semilogy("eps dual", xplot, yplot4, "k--");
	plt::xlabel("Number of Iteration");
	plt::ylabel("||s||_2");
	plt::legend();
	plt::save(figurename2);
#ifdef SHOW_FIGURE
   	plt::show();
   	plt::close();
#endif// endif SHOW_FIGURE
#endif// endif PLOT_FIGURE
}

