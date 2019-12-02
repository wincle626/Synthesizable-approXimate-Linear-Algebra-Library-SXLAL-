/*
 * admmlasso_positx.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_positx.hpp"


namespace plt = matplotlibcpp;


void ADMM_LASSO_POSITX(posit_2_t A[DIAG][DIAG], posit_2_t At[DIAG][DIAG],
					   posit_2_t invL[DIAG][DIAG],posit_2_t invU[DIAG][DIAG],
					   posit_2_t b[DIAG], posit_2_t lambda,
					   posit_2_t rho, posit_2_t alpha){

	// parameters
	posit_2_t oneminusalpha = pX2_sub(pX1_to_pX2(convertDoubleToPX1(1, TOTALBITS), TOTALBITS), alpha, TOTALBITS);
	posit_2_t lambdadivrho = pX2_div(lambda, rho, TOTALBITS);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit_2_t Atb[DIAG];
//	posit_2_t At[DIAG][DIAG];
//	posit_2_t AtA[DIAG][DIAG];
//	posit_2_t EYE[DIAG][DIAG];
//	posit_2_t rhoEYE[DIAG][DIAG];
//	posit_2_t AtAplusrhoeye[DIAG][DIAG];
//	posit_2_t L[DIAG][DIAG], U[DIAG][DIAG];
//	posit_2_t invL[DIAG][DIAG], invU[DIAG][DIAG];
	posit_2_t x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	posit_2_t zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	posit_2_t invLq[DIAG];
	posit_2_t alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	posit_2_t x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	posit_2_t x_hatz[DIAG];

	posit_2_t one = pX1_to_pX2(convertDoubleToPX1(1, TOTALBITS), TOTALBITS);
	posit_2_t zero = pX1_to_pX2(convertDoubleToPX1(0, TOTALBITS), TOTALBITS);

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit_2_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit_2_t,posit_1_t,DIAG,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit_2_t,DIAG,DIAG,pX2_mul,TOTALBITS>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit_2_t,DIAG,DIAG,pX2_add,TOTALBITS>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2,
//			pX2_add, pX2_sub, pX2_mul, pX2_div, pX2_sqrt, pX2_lt,TOTALBITS>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2,pX2_mul,
//			pX2_add, pX2_sub, pX2_sqrt, pX2_div, pX2_eq,TOTALBITS>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit_2_t,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(x);
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(z);
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
	struct history<posit_2_t> hist; // iteration record and early termination
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,
				convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		SoftPosit_Algebraobj.VEC_EQ<posit_2_t,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
				x, alpha, alphax);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
				zold, oneminusalpha, oneminusalphazold);
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu);
#endif
		SoftPosit_Algebraobj.VEC_SCALAR_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(
				x_hatu, lambdadivrho, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
				x_hatu, lambdadivrho, x_hatu2);
		SoftPosit_Algebraobj.VEC_MINUS<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2, pX2_mul,TOTALBITS>(
				x_hatu2, x_hatu2);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
				x_hatu1, zero, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
				x_hatu2, zero, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(u);
#endif
#if defined(RECORD_RESULT)
		posit_2_t znorm;
		posit_2_t Ax[DIAG], Axb[DIAG];
		posit_2_t Axbnorm2;
		posit_2_t xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		posit_2_t xznorm, rhozoldznorm;
		posit_2_t xnorm;
		posit_2_t rhou[DIAG];
		posit_2_t rhounorm;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t, posit_1_t, DIAG, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(Axb, Axbnorm2);
		hist.objval[k] = pX2_add(
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(0.5,TOTALBITS),TOTALBITS), Axbnorm2,TOTALBITS),
				pX2_mul(lambda, znorm,TOTALBITS),TOTALBITS);
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(xz, xznorm);
		hist.r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t, DIAG,pX2_mul,TOTALBITS>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(x, xnorm);
		hist.eps_pri[k] = pX2_add(
				pX2_mul(pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(DIAG,TOTALBITS),TOTALBITS),TOTALBITS),
						pX1_to_pX2(convertDoubleToPX1(ABSTOL,TOTALBITS),TOTALBITS),TOTALBITS),
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(RELTOL,TOTALBITS),TOTALBITS),(pX2_lt(znorm,xnorm)?xnorm:znorm),TOTALBITS),TOTALBITS);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul, TOTALBITS>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(rhou, rhounorm);
		hist.eps_dual[k] = pX2_add(
				pX2_mul(pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(DIAG,TOTALBITS),TOTALBITS),TOTALBITS),
						pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(ABSTOL,TOTALBITS),TOTALBITS),TOTALBITS),TOTALBITS),
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(RELTOL,TOTALBITS),TOTALBITS),rhounorm,TOTALBITS),TOTALBITS);
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			hist.u_hist[i][k] = u[i];
			hist.x_hist[i][k] = x[i];
			hist.z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertPX1ToDouble(pX2_to_pX1(hist.r_norm[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(hist.eps_pri[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(hist.s_norm[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(hist.eps_dual[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(hist.objval[k],TOTALBITS)));
#if defined(EARLY_TERMINATE)
		if(pX2_lt(hist.r_norm[k], hist.eps_pri[k]) &&
		  pX2_lt(hist.s_norm[k], hist.eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x);
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
		resultfile3 << convertPX1ToDouble(pX2_to_pX1(hist.objval[i],TOTALBITS)) << "\n";
		resultfile4 << convertPX1ToDouble(pX2_to_pX1(hist.r_norm[i],TOTALBITS)) << "\n";
		resultfile5 << convertPX1ToDouble(pX2_to_pX1(hist.s_norm[i],TOTALBITS)) << "\n";
		resultfile6 << convertPX1ToDouble(pX2_to_pX1(hist.eps_pri[i],TOTALBITS)) << "\n";
		resultfile7 << convertPX1ToDouble(pX2_to_pX1(hist.eps_dual[i],TOTALBITS)) << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << convertPX1ToDouble(pX2_to_pX1(hist.x_hist[j][i],TOTALBITS)) << ",";
			resultfile1 << convertPX1ToDouble(pX2_to_pX1(hist.u_hist[j][i],TOTALBITS)) << ",";
			resultfile2 << convertPX1ToDouble(pX2_to_pX1(hist.z_hist[j][i],TOTALBITS)) << ",";
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
		yplot.push_back(convertPX1ToDouble(pX2_to_pX1(hist.objval[i],TOTALBITS)));
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
		yplot1.push_back(convertPX1ToDouble(pX2_to_pX1(hist.r_norm[i],TOTALBITS)));
		yplot2.push_back(convertPX1ToDouble(pX2_to_pX1(hist.eps_pri[i],TOTALBITS)));
		yplot3.push_back(convertPX1ToDouble(pX2_to_pX1(hist.s_norm[i],TOTALBITS)));
		yplot4.push_back(convertPX1ToDouble(pX2_to_pX1(hist.eps_dual[i],TOTALBITS)));
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


void ADMM_LASSO_POSITX(posit_2_t **A, posit_2_t **At,
						   posit_2_t **invL,posit_2_t **invU,
						   posit_2_t *b, posit_2_t lambda,
						   posit_2_t rho, posit_2_t alpha){

	// parameters
	posit_2_t oneminusalpha = pX2_sub(pX1_to_pX2(convertDoubleToPX1(1, TOTALBITS), TOTALBITS), alpha, TOTALBITS);
	posit_2_t lambdadivrho = pX2_div(lambda, rho, TOTALBITS);

	// variables
	SoftPosit_Algebra SoftPosit_Algebraobj;
	DISPLAY DISPLAYobj;

	posit_2_t *Atb;
//	posit_2_t **At;
//	posit_2_t **AtA;
//	posit_2_t **EYE;
//	posit_2_t **rhoEYE;
//	posit_2_t **AtAplusrhoeye;
//	posit_2_t **L, **U;
//	posit_2_t **invL, **invU;
	posit_2_t *x, *zold, *z, *u;
	posit_2_t *zminusu, *rhozminusu, *q;
	posit_2_t *invLq;
	posit_2_t *alphax,*oneminusalphazold,*x_hat;
	posit_2_t *x_hatu,*x_hatu1,*x_hatu2;
	posit_2_t *x_hatz;
	Atb = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	zold = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	z = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	u = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	zminusu = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	rhozminusu = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	q = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	invLq = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	alphax = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	oneminusalphazold = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x_hat = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x_hatu = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x_hatu1 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x_hatu2 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
	x_hatz = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//	At = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	AtA = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	EYE = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	rhoEYE = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	AtAplusrhoeye = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	L = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	U = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	invL = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	invU = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		AtA[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		EYE[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		rhoEYE[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		AtAplusrhoeye[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		L[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		U[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		invL[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//		invU[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
//	}

	posit_2_t one = pX1_to_pX2(convertDoubleToPX1(1, TOTALBITS), TOTALBITS);
	posit_2_t zero = pX1_to_pX2(convertDoubleToPX1(0, TOTALBITS), TOTALBITS);

	// A'*b
//	SoftPosit_Algebraobj.MAT_TRANS<posit_2_t,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(At);
#endif
	SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	SoftPosit_Algebraobj.MAT_MUL<posit_2_t,posit_1_t,DIAG,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(AtA);
#endif
//	SoftPosit_Algebraobj.IDENDTITY_MAT<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(EYE);
//	SoftPosit_Algebraobj.MAT_SCALAR_DOTMUL<posit_2_t,DIAG,DIAG,pX2_mul,TOTALBITS>(
//			EYE, rho, rhoEYE);
//	SoftPosit_Algebraobj.MAT_ADD<posit_2_t,DIAG,DIAG,pX2_add,TOTALBITS>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(AtAplusrhoeye);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	SoftPosit_Algebraobj.LU_CHOLBANACHROUT<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2,
//			pX2_add, pX2_sub, pX2_mul, pX2_div, pX2_sqrt, pX2_lt,TOTALBITS>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(U);
#endif
	// invers L and U;
//	SoftPosit_Algebraobj.MAT_QRINV<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2,pX2_mul,
//			pX2_add, pX2_sub, pX2_sqrt, pX2_div, pX2_eq,TOTALBITS>(L, invL);
//	SoftPosit_Algebraobj.MAT_TRANS<posit_2_t,DIAG,DIAG>(invL, invU);
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	/// iteration
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(x);
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(z);
	SoftPosit_Algebraobj.ZEROS_VEC<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,TOTALBITS>(u);
#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	std::string clockname = "timeprofile.txt";
	TimeProfile.open(clockname);
#endif
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	posit_2_t** u_hist = NULL;
	posit_2_t** x_hist = NULL;
	posit_2_t** z_hist = NULL;
	posit_2_t* objval = NULL;
	posit_2_t* r_norm = NULL;
	posit_2_t* s_norm = NULL;
	posit_2_t* eps_pri = NULL;
	posit_2_t* eps_dual = NULL;
	u_hist = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	x_hist = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	z_hist = (posit_2_t**) malloc(sizeof(posit_2_t*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
		x_hist[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
		z_hist[i] = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	}
	objval = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	r_norm = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	s_norm = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	eps_pri = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	eps_dual = (posit_2_t*) malloc(sizeof(posit_2_t)*MAX_ITER);
	int k = 0;
	for(k=0;k<MAX_ITER;k++){
		// q = Atb + rho*(z - u);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(z, u, zminusu);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(zminusu,
				rho, rhozminusu);
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,
				convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<posit_2_t,posit_1_t,DIAG,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(invLq);
#endif
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// zold = z
		SoftPosit_Algebraobj.VEC_EQ<posit_2_t,DIAG>(z, zold);
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//  x_hat = alpha*x + (1 - alpha)*zold;
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
				x, alpha, alphax);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
				zold, oneminusalpha, oneminusalphazold);
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
				alphax, oneminusalphazold, x_hat);
#if defined(DEBUG_ITER)
		std::cout << "x_hat:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu);
#endif
		SoftPosit_Algebraobj.VEC_SCALAR_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(
				x_hatu, lambdadivrho, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
				x_hatu, lambdadivrho, x_hatu2);
		SoftPosit_Algebraobj.VEC_MINUS<posit_2_t,posit_1_t,DIAG, convertDoubleToPX1,pX1_to_pX2, pX2_mul,TOTALBITS>(
				x_hatu2, x_hatu2);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
				x_hatu1, zero, x_hatu1);
		SoftPosit_Algebraobj.VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
				x_hatu2, zero, x_hatu2);
#if defined(DEBUG_ITER)
		std::cout << "xhatu1:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatu2);
#endif
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x_hatz);
#endif
		SoftPosit_Algebraobj.VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(u);
#endif
#if defined(RECORD_RESULT)

		posit_2_t znorm;
		posit_2_t *Ax, *Axb;
		posit_2_t Axbnorm2;
		posit_2_t *xz, *zoldz, *rhozoldz;
		posit_2_t xznorm, rhozoldznorm;
		posit_2_t xnorm;
		posit_2_t rhou[DIAG];
		posit_2_t rhounorm;
		Ax = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		Axb = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		xz = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		zoldz = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		rhozoldz = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(z, znorm);
		SoftPosit_Algebraobj.MAT_VEC_MUL<posit_2_t, posit_1_t, DIAG, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(A, x, Ax);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(Ax, b, Axb);
		SoftPosit_Algebraobj.VEC_NORM2<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(Axb, Axbnorm2);
		objval[k] = pX2_add(
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(0.5,TOTALBITS),TOTALBITS), Axbnorm2,TOTALBITS),
				pX2_mul(lambda, znorm,TOTALBITS),TOTALBITS);
		// history.r_norm(k)  = norm(x - z);
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(x, z, xz);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(xz, xznorm);
		r_norm[k] = xznorm;
		// history.s_norm(k)  = norm(-rho*(z - zold));
		SoftPosit_Algebraobj.VEC_SUB<posit_2_t, DIAG,pX2_sub,TOTALBITS>(zold, z, zoldz);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t, DIAG,pX2_mul,TOTALBITS>(zoldz, rho, rhozoldz);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;
		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(x, xnorm);
		eps_pri[k] = pX2_add(
				pX2_mul(pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(DIAG,TOTALBITS),TOTALBITS),TOTALBITS),
						pX1_to_pX2(convertDoubleToPX1(ABSTOL,TOTALBITS),TOTALBITS),TOTALBITS),
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(RELTOL,TOTALBITS),TOTALBITS),(pX2_lt(znorm,xnorm)?xnorm:znorm),TOTALBITS),TOTALBITS);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		SoftPosit_Algebraobj.VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul, TOTALBITS>(rho, u, rhou);
		SoftPosit_Algebraobj.VEC_NORM<posit_2_t, posit_1_t, DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,pX2_sqrt,TOTALBITS>(rhou, rhounorm);
		eps_dual[k] = pX2_add(
				pX2_mul(pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(DIAG,TOTALBITS),TOTALBITS),TOTALBITS),
						pX2_sqrt(pX1_to_pX2(convertDoubleToPX1(ABSTOL,TOTALBITS),TOTALBITS),TOTALBITS),TOTALBITS),
				pX2_mul(pX1_to_pX2(convertDoubleToPX1(RELTOL,TOTALBITS),TOTALBITS),rhounorm,TOTALBITS),TOTALBITS);
		// record iterative solution
		for(int i=0;i<DIAG;i++){
			u_hist[i][k] = u[i];
			x_hist[i][k] = x[i];
			z_hist[i][k] = z[i];
		}
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", k,
				convertPX1ToDouble(pX2_to_pX1(r_norm[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(eps_pri[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(s_norm[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(eps_dual[k],TOTALBITS)),
				convertPX1ToDouble(pX2_to_pX1(objval[k],TOTALBITS)));
#if defined(EARLY_TERMINATE)
		if(pX2_lt(r_norm[k], eps_pri[k]) &&
		  pX2_lt(s_norm[k], eps_dual[k])){
			std::cout << k << "th iteration. Oho! Terminated! " << std::endl;
			break;
		}
#endif
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;

	}
#if defined(DEBUG_ITER)
	std::cout << "final x:" << std::endl;
	DISPLAYobj.printvector<posit_2_t,posit_1_t,DIAG,pX2_to_pX1,convertPX1ToDouble,TOTALBITS>(x);
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
		resultfile3 << convertPX1ToDouble(pX2_to_pX1(objval[i],TOTALBITS)) << "\n";
		resultfile4 << convertPX1ToDouble(pX2_to_pX1(r_norm[i],TOTALBITS)) << "\n";
		resultfile5 << convertPX1ToDouble(pX2_to_pX1(s_norm[i],TOTALBITS)) << "\n";
		resultfile6 << convertPX1ToDouble(pX2_to_pX1(eps_pri[i],TOTALBITS)) << "\n";
		resultfile7 << convertPX1ToDouble(pX2_to_pX1(eps_dual[i],TOTALBITS)) << "\n";
		for(int j=0;j<DIAG;j++){
			resultfile << convertPX1ToDouble(pX2_to_pX1(x_hist[j][i],TOTALBITS)) << ",";
			resultfile1 << convertPX1ToDouble(pX2_to_pX1(u_hist[j][i],TOTALBITS)) << ",";
			resultfile2 << convertPX1ToDouble(pX2_to_pX1(z_hist[j][i],TOTALBITS)) << ",";
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
		yplot.push_back(convertPX1ToDouble(pX2_to_pX1(objval[i],TOTALBITS)));
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
		yplot1.push_back(convertPX1ToDouble(pX2_to_pX1(r_norm[i],TOTALBITS)));
		yplot2.push_back(convertPX1ToDouble(pX2_to_pX1(eps_pri[i],TOTALBITS)));
		yplot3.push_back(convertPX1ToDouble(pX2_to_pX1(s_norm[i],TOTALBITS)));
		yplot4.push_back(convertPX1ToDouble(pX2_to_pX1(eps_dual[i],TOTALBITS)));
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

