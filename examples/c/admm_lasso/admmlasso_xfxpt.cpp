/*
 * admmlasso_xfxpt.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */


#include "admmlasso_xfxpt.hpp"

namespace plt = matplotlibcpp;


void ADMM_LASSO_FXPT(DATA_IN_T A[DIAG][DIAG], DATA_IN_T At[DIAG][DIAG],
					   DATA_IN_T invL[DIAG][DIAG],DATA_IN_T invU[DIAG][DIAG],
					   DATA_IN_T b[DIAG], DATA_IN_T lambda,
					   DATA_IN_T rho, DATA_IN_T alpha){

	// parameters
	DATA_IN_T oneminusalpha = 1 - alpha;
	DATA_IN_T lambdadivrho = lambda / rho;

	// variables
	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebraobj;
	DISPLAY DISPLAYobj;

	DATA_IN_T Atb[DIAG];
//	DATA_IN_T At[DIAG][DIAG];
//	DATA_IN_T AtA[DIAG][DIAG];
//	DATA_IN_T EYE[DIAG][DIAG];
//	DATA_IN_T rhoEYE[DIAG][DIAG];
//	DATA_IN_T AtAplusrhoeye[DIAG][DIAG];
//	DATA_IN_T L[DIAG][DIAG], U[DIAG][DIAG];
//	DATA_IN_T invL[DIAG][DIAG], invU[DIAG][DIAG];
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
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(At);
#endif
	Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<DATA_IN_T,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Xilinx_Fixed_Point_Algebraobj.MAT_MUL<DATA_IN_T,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(AtA);
#endif
//	Xilinx_Fixed_Point_Algebraobj.IDENDTITY_MAT<DATA_IN_T,DIAG,DIAG>(EYE);
//	Xilinx_Fixed_Point_Algebraobj.MAT_SCALAR_DOTMUL<DATA_IN_T,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Xilinx_Fixed_Point_Algebraobj.MAT_ADD<DATA_IN_T,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Xilinx_Fixed_Point_Algebraobj.LU_CHOLBANACHROUT<DATA_IN_T,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Xilinx_Fixed_Point_Algebraobj.MAT_QRINV<DATA_IN_T,DIAG>(L, invL);
//	Xilinx_Fixed_Point_Algebraobj.MAT_TRANS<DATA_IN_T,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	struct history<DATA_IN_T> hist; // iteration record and early termination
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(z, u, zminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T,DIAG>(zminusu,
				rho, rhozminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(invLq);
#endif
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x);
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
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu2);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatz);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(u);
#endif
#if defined(RECORD_RESULT)
		DATA_IN_T znorm;
		DATA_IN_T Ax[DIAG], Axb[DIAG];
		DATA_IN_T Axbnorm2;
		DATA_IN_T xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		DATA_IN_T xznorm, rhozoldznorm;
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
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(zold, z, zoldz);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>(zoldz, rho, rhozoldz);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			(float)hist.r_norm[k], (float)hist.eps_pri[k],
			(float)hist.s_norm[k], (float)hist.eps_dual[k],
			(float)hist.objval[k]);
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
	DISPLAYobj.printvector<DATA_IN_T,DIAG>(x);
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
		resultfile3 << (float)(hist.objval[i]) << "\n";
		resultfile4 << (float)(hist.r_norm[i]) << "\n";
		resultfile5 << (float)(hist.s_norm[i]) << "\n";
		resultfile6 << (float)(hist.eps_pri[i]) << "\n";
		resultfile7 << (float)(hist.eps_dual[i]) << "\n";
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
		yplot.push_back((float)hist.objval[i]);
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
		yplot1.push_back((float)hist.r_norm[i]);
		yplot2.push_back((float)hist.eps_pri[i]);
		yplot3.push_back((float)hist.s_norm[i]);
		yplot4.push_back((float)hist.eps_dual[i]);
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

void ADMM_LASSO_FXPT(DATA_IN_T **A, DATA_IN_T **At,
					   DATA_IN_T **invL,DATA_IN_T **invU,
					   DATA_IN_T *b, DATA_IN_T lambda,
					   DATA_IN_T rho, DATA_IN_T alpha){

	// parameters
	DATA_IN_T oneminusalpha = 1 - alpha;
	DATA_IN_T lambdadivrho = lambda / rho;

	// variables
	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebraobj;
	DISPLAY DISPLAYobj;

	DATA_IN_T *Atb;
//	DATA_IN_T **At;
//	DATA_IN_T **AtA;
//	DATA_IN_T **EYE;
//	DATA_IN_T **rhoEYE;
//	DATA_IN_T **AtAplusrhoeye;
//	DATA_IN_T **L, **U;
//	DATA_IN_T **invL, **invU;
	DATA_IN_T *x, *zold, *z, *u;
	DATA_IN_T *zminusu, *rhozminusu, *q;
	DATA_IN_T *invLq;
	DATA_IN_T *alphax,*oneminusalphazold,*x_hat;
	DATA_IN_T *x_hatu,*x_hatu1,*x_hatu2;
	DATA_IN_T *x_hatz;
	Atb = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	zold = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	z = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	u = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	zminusu = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	rhozminusu = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	q = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	invLq = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	alphax = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	oneminusalphazold = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x_hat = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x_hatu = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x_hatu1 = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x_hatu2 = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
	x_hatz = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//	At = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	AtA = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	EYE = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	rhoEYE = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	AtAplusrhoeye = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	L = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	U = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	invL = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	invU = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		AtA[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		EYE[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		rhoEYE[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		AtAplusrhoeye[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		L[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		U[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		invL[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//		invU[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
//	}

	// A'*b
//	Xilinx_Fixed_Point_Algebraobj.MAT_TRANS<DATA_IN_T,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(At);
#endif
	Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<DATA_IN_T,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Xilinx_Fixed_Point_Algebraobj.MAT_MUL<DATA_IN_T,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(AtA);
#endif
//	Xilinx_Fixed_Point_Algebraobj.IDENDTITY_MAT<DATA_IN_T,DIAG,DIAG>(EYE);
//	Xilinx_Fixed_Point_Algebraobj.MAT_SCALAR_DOTMUL<DATA_IN_T,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Xilinx_Fixed_Point_Algebraobj.MAT_ADD<DATA_IN_T,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Xilinx_Fixed_Point_Algebraobj.LU_CHOLBANACHROUT<DATA_IN_T,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Xilinx_Fixed_Point_Algebraobj.MAT_QRINV<DATA_IN_T,DIAG>(L, invL);
//	Xilinx_Fixed_Point_Algebraobj.MAT_TRANS<DATA_IN_T,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	DATA_IN_T** u_hist = NULL;
	DATA_IN_T** x_hist = NULL;
	DATA_IN_T** z_hist = NULL;
	DATA_IN_T* objval = NULL;
	DATA_IN_T* r_norm = NULL;
	DATA_IN_T* s_norm = NULL;
	DATA_IN_T* eps_pri = NULL;
	DATA_IN_T* eps_dual = NULL;
	u_hist = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	x_hist = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	z_hist = (DATA_IN_T**) malloc(sizeof(DATA_IN_T*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
		x_hist[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
		z_hist[i] = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	}
	objval = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	r_norm = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	s_norm = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	eps_pri = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	eps_dual = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*MAX_ITER);
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(z, u, zminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T,DIAG>(zminusu,
				rho, rhozminusu);
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<DATA_IN_T,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(invLq);
#endif
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x);
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
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatu2);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(x_hatz);
#endif
		Xilinx_Fixed_Point_Algebraobj.VEC_ADD<DATA_IN_T,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<DATA_IN_T,DIAG>(u);
#endif
#if defined(RECORD_RESULT)

		DATA_IN_T znorm;
		DATA_IN_T *Ax, *Axb;
		DATA_IN_T Axbnorm2;
		DATA_IN_T *xz, *zoldz, *rhozoldz;
		DATA_IN_T xznorm, rhozoldznorm;
		DATA_IN_T xnorm;
		DATA_IN_T rhou[DIAG];
		DATA_IN_T rhounorm;
		Ax = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		Axb = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		xz = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		zoldz = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);
		rhozoldz = (DATA_IN_T*) malloc(sizeof(DATA_IN_T)*DIAG);

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(z, znorm);
		Xilinx_Fixed_Point_Algebraobj.MAT_VEC_MUL<DATA_IN_T, DIAG, DIAG>(A, x, Ax);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(Ax, b, Axb);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM2<DATA_IN_T, DIAG>(Axb, Axbnorm2);
		objval[k] = ((DATA_IN_T)0.5) * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(x, z, xz);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(xz, xznorm);
		r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Xilinx_Fixed_Point_Algebraobj.VEC_SUB<DATA_IN_T, DIAG>(zold, z, zoldz);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>(zoldz, rho, rhozoldz);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(x, xnorm);
		eps_pri[k] = ((DATA_IN_T)std::sqrt(DIAG))*((DATA_IN_T)ABSTOL)+((DATA_IN_T)RELTOL)*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Xilinx_Fixed_Point_Algebraobj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>(rho, u, rhou);
		Xilinx_Fixed_Point_Algebraobj.VEC_NORM<DATA_IN_T, DIAG>(rhou, rhounorm);
		eps_dual[k] = ((DATA_IN_T)std::sqrt(DIAG))*((DATA_IN_T)ABSTOL)+((DATA_IN_T)RELTOL)*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			(float)r_norm[k], (float)eps_pri[k],
			(float)s_norm[k], (float)eps_dual[k],
			(float)objval[k]);
#if defined(EARLY_TERMINATE)
		if((r_norm[k] < eps_pri[k]) &&
		  (s_norm[k] < eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<DATA_IN_T,DIAG>(x);
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
		resultfile3 << (float)(objval[i]) << "\n";
		resultfile4 << (float)(r_norm[i]) << "\n";
		resultfile5 << (float)(s_norm[i]) << "\n";
		resultfile6 << (float)(eps_pri[i]) << "\n";
		resultfile7 << (float)(eps_dual[i]) << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << (float)x_hist[j][i] << ",";
			resultfile1 << (float)u_hist[j][i] << ",";
			resultfile2 << (float)z_hist[j][i] << ",";
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
		yplot.push_back((float)objval[i]);
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
		yplot1.push_back((float)r_norm[i]);
		yplot2.push_back((float)eps_pri[i]);
		yplot3.push_back((float)s_norm[i]);
		yplot4.push_back((float)eps_dual[i]);
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
