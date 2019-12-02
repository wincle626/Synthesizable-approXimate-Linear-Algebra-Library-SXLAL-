/*
 * admmlasso_xfloat.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_xfloat.hpp"

namespace plt = matplotlibcpp;

void ADMM_LASSO_XFPT(fptx_admmlasso A[DIAG][DIAG], fptx_admmlasso At[DIAG][DIAG],
					   fptx_admmlasso invL[DIAG][DIAG],fptx_admmlasso invU[DIAG][DIAG],
					   fptx_admmlasso b[DIAG], fptx_admmlasso lambda,
					   fptx_admmlasso rho, fptx_admmlasso alpha){

	// parameters
	fptx_admmlasso oneminusalpha = 1 - alpha;
	fptx_admmlasso lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	DISPLAY DISPLAYobj;

	fptx_admmlasso Atb[DIAG];
//	fptx_admmlasso At[DIAG][DIAG];
//	fptx_admmlasso AtA[DIAG][DIAG];
//	fptx_admmlasso EYE[DIAG][DIAG];
//	fptx_admmlasso rhoEYE[DIAG][DIAG];
//	fptx_admmlasso AtAplusrhoeye[DIAG][DIAG];
//	fptx_admmlasso L[DIAG][DIAG], U[DIAG][DIAG];
//	fptx_admmlasso invL[DIAG][DIAG], invU[DIAG][DIAG];
	fptx_admmlasso x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	fptx_admmlasso zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	fptx_admmlasso invLq[DIAG];
	fptx_admmlasso alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	fptx_admmlasso x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	fptx_admmlasso x_hatz[DIAG];

	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<fptx_admmlasso,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<fptx_admmlasso,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(AtA);
#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<fptx_admmlasso,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<fptx_admmlasso,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<fptx_admmlasso,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<fptx_admmlasso,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<fptx_admmlasso,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
    int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x);
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
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(u);
#endif
#if defined(RECORD_RESULT)
		fptx_admmlasso znorm;
		fptx_admmlasso Ax[DIAG], Axb[DIAG];
		fptx_admmlasso Axbnorm2;
		fptx_admmlasso xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		fptx_admmlasso xznorm, rhozoldznorm;
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
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			(float)hist.r_norm[k], (float)hist.eps_pri[k],
			(float)hist.s_norm[k], (float)hist.eps_dual[k], (float)hist.objval[k]);
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
	DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x);
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
	for(int i=0; i<MAX_ITER; i++){
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



void ADMM_LASSO_XFPT(fptx_admmlasso **A, fptx_admmlasso **At,
					   fptx_admmlasso **invL,fptx_admmlasso **invU,
					   fptx_admmlasso *b, fptx_admmlasso lambda,
					   fptx_admmlasso rho, fptx_admmlasso alpha){

	// parameters
	fptx_admmlasso oneminusalpha = 1 - alpha;
	fptx_admmlasso lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	DISPLAY DISPLAYobj;

	fptx_admmlasso *Atb;
//	fptx_admmlasso **At;
//	fptx_admmlasso **AtA;
//	fptx_admmlasso **EYE;
//	fptx_admmlasso **rhoEYE;
//	fptx_admmlasso **AtAplusrhoeye;
//	fptx_admmlasso **L, **U;
//	fptx_admmlasso **invL, **invU;
	fptx_admmlasso *x, *zold, *z, *u;
	fptx_admmlasso *zminusu, *rhozminusu, *q;
	fptx_admmlasso *invLq;
	fptx_admmlasso *alphax,*oneminusalphazold,*x_hat;
	fptx_admmlasso *x_hatu,*x_hatu1,*x_hatu2;
	fptx_admmlasso *x_hatz;
	Atb = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	zold = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	z = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	u = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	zminusu = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	rhozminusu = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	q = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	invLq = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	alphax = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	oneminusalphazold = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x_hat = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x_hatu = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x_hatu1 = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x_hatu2 = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
	x_hatz = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//	At = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	AtA = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	EYE = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	rhoEYE = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	AtAplusrhoeye = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	L = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	U = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	invL = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	invU = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		AtA[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		EYE[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		rhoEYE[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		AtAplusrhoeye[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		L[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		U[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		invL[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//		invU[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
//	}

	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<fptx_admmlasso,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<fptx_admmlasso,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(AtA);
#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<fptx_admmlasso,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<fptx_admmlasso,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<fptx_admmlasso,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<fptx_admmlasso,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<fptx_admmlasso,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<fptx_admmlasso,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");

	fptx_admmlasso** u_hist = NULL;
	fptx_admmlasso** x_hist = NULL;
	fptx_admmlasso** z_hist = NULL;
	fptx_admmlasso* objval = NULL;
	fptx_admmlasso* r_norm = NULL;
	fptx_admmlasso* s_norm = NULL;
	fptx_admmlasso* eps_pri = NULL;
	fptx_admmlasso* eps_dual = NULL;
	u_hist = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	x_hist = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	z_hist = (fptx_admmlasso**) malloc(sizeof(fptx_admmlasso*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
		x_hist[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
		z_hist[i] = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
	}
	objval = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
	r_norm = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
	s_norm = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
	eps_pri = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);
	eps_dual = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*MAX_ITER);

	int k=0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(invL);

		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<fptx_admmlasso,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x);
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
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<fptx_admmlasso,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<fptx_admmlasso,DIAG>(u);
#endif
#if defined(RECORD_RESULT)

		fptx_admmlasso znorm;
		fptx_admmlasso *Ax, *Axb;
		fptx_admmlasso Axbnorm2;
		fptx_admmlasso *xz, *zoldz, *rhozoldz;
		fptx_admmlasso xznorm, rhozoldznorm;
		fptx_admmlasso xnorm;
		fptx_admmlasso rhou[DIAG];
		fptx_admmlasso rhounorm;
		Ax = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		Axb = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		xz = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		zoldz = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);
		rhozoldz = (fptx_admmlasso*) malloc(sizeof(fptx_admmlasso)*DIAG);

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<fptx_admmlasso, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<fptx_admmlasso, DIAG>(Axb, Axbnorm2);
		objval[k] = 0.5 * Axbnorm2 + lambda * znorm;
		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(xz, xznorm);
		r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SUB<fptx_admmlasso, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(x, xnorm);
		eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<fptx_admmlasso, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<fptx_admmlasso, DIAG>(rhou, rhounorm);
		eps_dual[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*rhounorm;
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
			(float)r_norm[k], (float)eps_pri[k],
			(float)s_norm[k], (float)eps_dual[k], (float)objval[k]);
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
	DISPLAYobj.printvector<fptx_admmlasso,DIAG>(x);
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
	for(int i=0; i<MAX_ITER; i++){
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


