/*
 * admmlasso_posit8.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_posit8.hpp"

namespace plt = matplotlibcpp;



void ADMM_LASSO_POSIT8(posit8_t A[DIAG][DIAG], posit8_t At[DIAG][DIAG],
					   posit8_t invL[DIAG][DIAG],posit8_t invU[DIAG][DIAG],
					   posit8_t b[DIAG], posit8_t lambda,
					   posit8_t rho, posit8_t alpha){

	// parameters
	posit8_t oneminusalpha = p8_sub(convertDoubleToP8(1), alpha);
	posit8_t lambdadivrho = p8_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit8_t Atb[DIAG];
//	posit8_t At[DIAG][DIAG];
//	posit8_t AtA[DIAG][DIAG];
//	posit8_t EYE[DIAG][DIAG];
//	posit8_t rhoEYE[DIAG][DIAG];
//	posit8_t AtAplusrhoeye[DIAG][DIAG];
//	posit8_t L[DIAG][DIAG], U[DIAG][DIAG];
//	posit8_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit8_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit8_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit8_t invLq[DIAG];
	posit8_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit8_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit8_t x_hatz[DIAG];

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit8_t,DIAG,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit8_t,DIAG,DIAG,convertDoubleToP8>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit8_t,DIAG,DIAG,p8_mul>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit8_t,DIAG,DIAG,p8_add>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit8_t,DIAG, convertDoubleToP8, p8_add, p8_sub, p8_mul, p8_div, p8_sqrt, p8_lt>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit8_t,DIAG, convertDoubleToP8, p8_mul, p8_add,
//			p8_sub, p8_sqrt, p8_div,p8_eq>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	struct history<posit8_t> hist; // iteration record and early termination
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,
				convertDoubleToP8,p8_mul,p8_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x);
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
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu);
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
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(u);
#endif
#if defined(RECORD_RESULT)
		posit8_t znorm;
		posit8_t Ax[DIAG], Axb[DIAG];
		posit8_t Axbnorm2;
		posit8_t xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		posit8_t xznorm, rhozoldznorm;
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
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t, DIAG,p8_mul>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertP8ToDouble(hist.r_norm[k]), convertP8ToDouble(hist.eps_pri[k]),
				convertP8ToDouble(hist.s_norm[k]), convertP8ToDouble(hist.eps_dual[k]),
				convertP8ToDouble(hist.objval[k]));
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
	DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x);
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
		resultfile3 << convertP8ToDouble(hist.objval[i]) << "\n";
		resultfile4 << convertP8ToDouble(hist.r_norm[i]) << "\n";
		resultfile5 << convertP8ToDouble(hist.s_norm[i]) << "\n";
		resultfile6 << convertP8ToDouble(hist.eps_pri[i]) << "\n";
		resultfile7 << convertP8ToDouble(hist.eps_dual[i]) << "\n";
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
		yplot.push_back(convertP8ToDouble(hist.objval[i]));
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
		yplot1.push_back(convertP8ToDouble(hist.r_norm[i]));
		yplot2.push_back(convertP8ToDouble(hist.eps_pri[i]));
		yplot3.push_back(convertP8ToDouble(hist.s_norm[i]));
		yplot4.push_back(convertP8ToDouble(hist.eps_dual[i]));
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

void ADMM_LASSO_POSIT8(posit8_t **A, posit8_t **At,
						   posit8_t **invL,posit8_t **invU,
						   posit8_t *b, posit8_t lambda,
						   posit8_t rho, posit8_t alpha){

	// parameters
	posit8_t oneminusalpha = p8_sub(convertDoubleToP8(1), alpha);
	posit8_t lambdadivrho = p8_div(lambda, rho);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit8_t *Atb;
//	posit8_t **At;
//	posit8_t **AtA;
//	posit8_t **EYE;
//	posit8_t **rhoEYE;
//	posit8_t **AtAplusrhoeye;
//	posit8_t **L, **U;
//	posit8_t **invL, **invU;
	posit8_t *x, *zold, *z, *u;
	posit8_t *zminusu, *rhozminusu, *q;
	posit8_t *invLq;
	posit8_t *alphax,*oneminusalphazold,*x_hat;
	posit8_t *x_hatu,*x_hatu1,*x_hatu2;
	posit8_t *x_hatz;
	Atb = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	zold = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	z = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	u = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	zminusu = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	rhozminusu = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	q = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	invLq = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	alphax = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	oneminusalphazold = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x_hat = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x_hatu = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x_hatu1 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x_hatu2 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
	x_hatz = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//	At = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	AtA = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	EYE = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	rhoEYE = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	AtAplusrhoeye = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	L = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	U = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	invL = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	invU = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		AtA[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		EYE[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		rhoEYE[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		AtAplusrhoeye[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		L[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		U[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		invL[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//		invU[i] = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
//	}

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit8_t,DIAG,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit8_t,DIAG,DIAG,convertDoubleToP8>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit8_t,DIAG,DIAG,p8_mul>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit8_t,DIAG,DIAG,p8_add>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit8_t,DIAG, convertDoubleToP8, p8_add, p8_sub, p8_mul, p8_div, p8_sqrt, p8_lt>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit8_t,DIAG, convertDoubleToP8, p8_mul,
//			p8_add, p8_sub, p8_sqrt, p8_div, p8_eq>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit8_t,DIAG,DIAG>(invL, invU);
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
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	posit8_t** u_hist = NULL;
	posit8_t** x_hist = NULL;
	posit8_t** z_hist = NULL;
	posit8_t* objval = NULL;
	posit8_t* r_norm = NULL;
	posit8_t* s_norm = NULL;
	posit8_t* eps_pri = NULL;
	posit8_t* eps_dual = NULL;
	u_hist = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	x_hist = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	z_hist = (posit8_t**) malloc(sizeof(posit8_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
		x_hist[i] = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
		z_hist[i] = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	}
	objval = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	r_norm = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	s_norm = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	eps_pri = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	eps_dual = (posit8_t*) malloc(sizeof(posit8_t)*MAX_ITER);
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,
				convertDoubleToP8,p8_mul,p8_add>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit8_t,DIAG,DIAG,convertP8ToDouble>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x);
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
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu);
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
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t,DIAG,p8_sub>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit8_t,DIAG,p8_add>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(u);
#endif
#if defined(RECORD_RESULT)

		posit8_t znorm;
		posit8_t *Ax, *Axb;
		posit8_t Axbnorm2;
		posit8_t *xz, *zoldz, *rhozoldz;
		posit8_t xznorm, rhozoldznorm;
		posit8_t xnorm;
		posit8_t rhou[DIAG];
		posit8_t rhounorm;
		Ax = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		Axb = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		xz = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		zoldz = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		rhozoldz = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit8_t, DIAG, DIAG,convertDoubleToP8,p8_mul,p8_add>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add>(Axb, Axbnorm2);
		objval[k] = p8_add(p8_mul(convertDoubleToP8(0.5), Axbnorm2), p8_mul(lambda, znorm));
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(xz, xznorm);
		r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SUB<posit8_t, DIAG,p8_sub>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t, DIAG,p8_mul>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit8_t, DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(x, xnorm);
		eps_pri[k] = p8_add(p8_mul(p8_sqrt(convertDoubleToP8(DIAG)),convertDoubleToP8(ABSTOL)),p8_mul(convertDoubleToP8(RELTOL),(p8_lt(znorm,xnorm)?xnorm:znorm)));
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit8_t,DIAG,p8_mul>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit8_t,DIAG,convertDoubleToP8,p8_mul,p8_add,p8_sqrt>(rhou, rhounorm);
		eps_dual[k] = p8_add(p8_mul(p8_sqrt(convertDoubleToP8(DIAG)),convertDoubleToP8(ABSTOL)),p8_mul(convertDoubleToP8(RELTOL),rhounorm));
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertP8ToDouble(r_norm[k]), convertP8ToDouble(eps_pri[k]),
				convertP8ToDouble(s_norm[k]), convertP8ToDouble(eps_dual[k]),
				convertP8ToDouble(objval[k]));
#if defined(EARLY_TERMINATE)
		if(p8_lt(r_norm[k], eps_pri[k]) &&
		  p8_lt(s_norm[k], eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<posit8_t,DIAG,convertP8ToDouble>(x);
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
		resultfile3 << convertP8ToDouble(objval[i]) << "\n";
		resultfile4 << convertP8ToDouble(r_norm[i]) << "\n";
		resultfile5 << convertP8ToDouble(s_norm[i]) << "\n";
		resultfile6 << convertP8ToDouble(eps_pri[i]) << "\n";
		resultfile7 << convertP8ToDouble(eps_dual[i]) << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << convertP8ToDouble(x_hist[j][i]) << ",";
			resultfile1 << convertP8ToDouble(u_hist[j][i]) << ",";
			resultfile2 << convertP8ToDouble(z_hist[j][i]) << ",";
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
		yplot.push_back(convertP8ToDouble(objval[i]));
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
		yplot1.push_back(convertP8ToDouble(r_norm[i]));
		yplot2.push_back(convertP8ToDouble(eps_pri[i]));
		yplot3.push_back(convertP8ToDouble(s_norm[i]));
		yplot4.push_back(convertP8ToDouble(eps_dual[i]));
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


