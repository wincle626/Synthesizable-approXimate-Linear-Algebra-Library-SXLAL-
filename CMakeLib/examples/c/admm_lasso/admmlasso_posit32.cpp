/*
 * admmlasso_posit32.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_posit32.hpp"

namespace plt = matplotlibcpp;


void ADMM_LASSO_POSIT32(posit32_t A[DIAG][DIAG], posit32_t At[DIAG][DIAG],
					   posit32_t invL[DIAG][DIAG],posit32_t invU[DIAG][DIAG],
					   posit32_t b[DIAG], posit32_t lambda,
					   posit32_t rho, posit32_t alpha){

	// parameters
	posit32_t oneminusalpha = p32_sub(convertDoubleToP32(1), alpha);
	posit32_t lambdadivrho = p32_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit32_t Atb[DIAG];
//	posit32_t At[DIAG][DIAG];
//	posit32_t AtA[DIAG][DIAG];
//	posit32_t EYE[DIAG][DIAG];
//	posit32_t rhoEYE[DIAG][DIAG];
//	posit32_t AtAplusrhoeye[DIAG][DIAG];
//	posit32_t L[DIAG][DIAG], U[DIAG][DIAG];
//	posit32_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit32_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit32_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit32_t invLq[DIAG];
	posit32_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit32_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit32_t x_hatz[DIAG];

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit32_t,DIAG,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit32_t,DIAG,DIAG,convertDoubleToP32>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit32_t,DIAG,DIAG,p32_mul>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit32_t,DIAG,DIAG,p32_add>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit32_t,DIAG, convertDoubleToP32, p32_add, p32_sub, p32_mul, p32_div, p32_sqrt, p32_lt>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit32_t,DIAG, convertDoubleToP32, p32_mul, p32_add,
//			p32_sub, p32_sqrt, p32_div,p32_eq>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	struct history<posit32_t> hist; // iteration record and early termination
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,
				convertDoubleToP32,p32_mul,p32_add>(
				invL, q, invLq);
#if defined(DEBUG_ITERADMM_LASSO_POSIT32)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x);
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
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu);
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
		DISPLAYobj.priADMM_LASSO_POSIT32ntvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(u);
#endif
#if defined(RECORD_RESULT)
		posit32_t znorm;
		posit32_t Ax[DIAG], Axb[DIAG];
		posit32_t Axbnorm2;
		posit32_t xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		posit32_t xznorm, rhozoldznorm;
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
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t, DIAG,p32_mul>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertP32ToDouble(hist.r_norm[k]), convertP32ToDouble(hist.eps_pri[k]),
				convertP32ToDouble(hist.s_norm[k]), convertP32ToDouble(hist.eps_dual[k]),
				convertP32ToDouble(hist.objval[k]));
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
	DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x);
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
		resultfile3 << convertP32ToDouble(hist.objval[i]) << "\n";
		resultfile4 << convertP32ToDouble(hist.r_norm[i]) << "\n";
		resultfile5 << convertP32ToDouble(hist.s_norm[i]) << "\n";
		resultfile6 << convertP32ToDouble(hist.eps_pri[i]) << "\n";
		resultfile7 << convertP32ToDouble(hist.eps_dual[i]) << "\n";
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
		yplot.push_back(convertP32ToDouble(hist.objval[i]));
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
		yplot1.push_back(convertP32ToDouble(hist.r_norm[i]));
		yplot2.push_back(convertP32ToDouble(hist.eps_pri[i]));
		yplot3.push_back(convertP32ToDouble(hist.s_norm[i]));
		yplot4.push_back(convertP32ToDouble(hist.eps_dual[i]));
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

void ADMM_LASSO_POSIT32(posit32_t **A, posit32_t **At,
						   posit32_t **invL,posit32_t **invU,
						   posit32_t *b, posit32_t lambda,
						   posit32_t rho, posit32_t alpha){

	// parameters
	posit32_t oneminusalpha = p32_sub(convertDoubleToP32(1), alpha);
	posit32_t lambdadivrho = p32_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit32_t *Atb;
//	posit32_t **At;
//	posit32_t **AtA;
//	posit32_t **EYE;
//	posit32_t **rhoEYE;
//	posit32_t **AtAplusrhoeye;
//	posit32_t **L, **U;
//	posit32_t **invL, **invU;
	posit32_t *x, *zold, *z, *u;
	posit32_t *zminusu, *rhozminusu, *q;
	posit32_t *invLq;
	posit32_t *alphax,*oneminusalphazold,*x_hat;
	posit32_t *x_hatu,*x_hatu1,*x_hatu2;
	posit32_t *x_hatz;
	Atb = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	zold = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	z = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	u = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	zminusu = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	rhozminusu = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	q = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	invLq = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	alphax = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	oneminusalphazold = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x_hat = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x_hatu = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x_hatu1 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x_hatu2 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
	x_hatz = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//	At = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	AtA = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	EYE = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	rhoEYE = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	AtAplusrhoeye = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	L = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	U = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	invL = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	invU = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		AtA[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		EYE[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		rhoEYE[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		AtAplusrhoeye[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		L[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		U[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		invL[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//		invU[i] = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
//	}

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit32_t,DIAG,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit32_t,DIAG,DIAG,convertDoubleToP32>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit32_t,DIAG,DIAG,p32_mul>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit32_t,DIAG,DIAG,p32_add>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit32_t,DIAG, convertDoubleToP32, p32_add, p32_sub, p32_mul, p32_div, p32_sqrt, p32_lt>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit32_t,DIAG, convertDoubleToP32, p32_mul,
//			p32_add, p32_sub, p32_sqrt, p32_div, p32_eq>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit32_t,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	posit32_t** u_hist = NULL;
	posit32_t** x_hist = NULL;
	posit32_t** z_hist = NULL;
	posit32_t* objval = NULL;
	posit32_t* r_norm = NULL;
	posit32_t* s_norm = NULL;
	posit32_t* eps_pri = NULL;
	posit32_t* eps_dual = NULL;
	u_hist = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	x_hist = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	z_hist = (posit32_t**) malloc(sizeof(posit32_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
		x_hist[i] = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
		z_hist[i] = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
	}
	objval = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
	r_norm = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
	s_norm = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
	eps_pri = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);
	eps_dual = (posit32_t*) malloc(sizeof(posit32_t)*MAX_ITER);

	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,
				convertDoubleToP32,p32_mul,p32_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit32_t,DIAG,DIAG,convertP32ToDouble>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x);
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
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu);
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
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t,DIAG,p32_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit32_t,DIAG,p32_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(u);
#endif
#if defined(RECORD_RESULT)

		posit32_t znorm;
		posit32_t *Ax, *Axb;
		posit32_t Axbnorm2;
		posit32_t *xz, *zoldz, *rhozoldz;
		posit32_t xznorm, rhozoldznorm;
		posit32_t xnorm;
		posit32_t rhou[DIAG];
		posit32_t rhounorm;
		Ax = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		Axb = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		xz = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		zoldz = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		rhozoldz = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit32_t, DIAG, DIAG,convertDoubleToP32,p32_mul,p32_add>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add>(Axb, Axbnorm2);
		objval[k] = p32_add(p32_mul(convertDoubleToP32(0.5), Axbnorm2), p32_mul(lambda, znorm));
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(xz, xznorm);
		r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SUB<posit32_t, DIAG,p32_sub>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t, DIAG,p32_mul>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit32_t, DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(x, xnorm);
		eps_pri[k] = p32_add(p32_mul(p32_sqrt(convertDoubleToP32(DIAG)),convertDoubleToP32(ABSTOL)),p32_mul(convertDoubleToP32(RELTOL),(p32_lt(znorm,xnorm)?xnorm:znorm)));
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit32_t,DIAG,convertDoubleToP32,p32_mul,p32_add,p32_sqrt>(rhou, rhounorm);
		eps_dual[k] = p32_add(p32_mul(p32_sqrt(convertDoubleToP32(DIAG)),convertDoubleToP32(ABSTOL)),p32_mul(convertDoubleToP32(RELTOL),rhounorm));
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertP32ToDouble(r_norm[k]), convertP32ToDouble(eps_pri[k]),
				convertP32ToDouble(s_norm[k]), convertP32ToDouble(eps_dual[k]),
				convertP32ToDouble(objval[k]));
#if defined(EARLY_TERMINATE)
		if(p32_lt(r_norm[k], eps_pri[k]) &&
		  p32_lt(s_norm[k], eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<posit32_t,DIAG,convertP32ToDouble>(x);
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
		resultfile3 << convertP32ToDouble(objval[i]) << "\n";
		resultfile4 << convertP32ToDouble(r_norm[i]) << "\n";
		resultfile5 << convertP32ToDouble(s_norm[i]) << "\n";
		resultfile6 << convertP32ToDouble(eps_pri[i]) << "\n";
		resultfile7 << convertP32ToDouble(eps_dual[i]) << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << convertP32ToDouble(x_hist[j][i]) << ",";
			resultfile1 << convertP32ToDouble(u_hist[j][i]) << ",";
			resultfile2 << convertP32ToDouble(z_hist[j][i]) << ",";
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
		yplot.push_back(convertP32ToDouble(objval[i]));
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
		yplot1.push_back(convertP32ToDouble(r_norm[i]));
		yplot2.push_back(convertP32ToDouble(eps_pri[i]));
		yplot3.push_back(convertP32ToDouble(s_norm[i]));
		yplot4.push_back(convertP32ToDouble(eps_dual[i]));
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


