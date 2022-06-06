/*
 * admmlasso_gfloat.cpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#include "admmlasso_gfloat.hpp"

namespace plt = matplotlibcpp;

void ADMM_LASSO_FLOAT(float A[DIAG][DIAG], float At[DIAG][DIAG],
					   float invL[DIAG][DIAG],float invU[DIAG][DIAG],
					   float b[DIAG], float lambda,
					   float rho, float alpha){

	// parameters
	float oneminusalpha = 1 - alpha;
	float lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
	DISPLAY DISPLAYobj;

	float Atb[DIAG];
//	float At[DIAG][DIAG];
//	float AtA[DIAG][DIAG];
//	float EYE[DIAG][DIAG];
//	float rhoEYE[DIAG][DIAG];
//	float AtAplusrhoeye[DIAG][DIAG];
//	float L[DIAG][DIAG], U[DIAG][DIAG];
//	float invL[DIAG][DIAG], invU[DIAG][DIAG];
	float x[DIAG], zold[DIAG], z[DIAG], u[DIAG];
	float zminusu[DIAG], rhozminusu[DIAG], q[DIAG];
	float invLq[DIAG];
	float alphax[DIAG],oneminusalphazold[DIAG],x_hat[DIAG];
	float x_hatu[DIAG],x_hatu1[DIAG],x_hatu2[DIAG];
	float x_hatz[DIAG];
#if defined(DEBUG_ITER)
	std::cout << "A:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(A);
	std::cout << "U_inv:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(invU);
	std::cout << "L_inv:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(invL);
#endif

	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<float,DIAG>(Atb);
#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<float,DIAG,DIAG>(At, A, AtA);
//#if defined(DEBUG_ITER)
//	std::cout << "AtA:" << std::endl;
//	DISPLAYobj.printmatrix<float,DIAG,DIAG>(AtA);
//#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<float,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<float,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<float,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
//#if defined(DEBUG_ITER)
//	std::cout << "AtAplusrhoeye:" << std::endl;
//	DISPLAYobj.printmatrix<float,DIAG,DIAG>(AtAplusrhoeye);
//#endif
	// std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<float,DIAG>(
//			AtAplusrhoeye, L, U);
//#if defined(DEBUG_ITER)
//	std::cout << "L:" << std::endl;
//	DISPLAYobj.printmatrix<float,DIAG,DIAG>(L);
//	std::cout << "U:" << std::endl;
//	DISPLAYobj.printmatrix<float,DIAG,DIAG>(U);
//#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<float,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(invL, invU);
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
	struct histry hist; // iteration record and early termination
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
    int k = 0;
	for(k=0;k<MAX_ITER;k++){
//	for(int k=0;k<1;k++){
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(q);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<float,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<float,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x);
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
		DISPLAYobj.printvector<float,DIAG>(x_hat);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<float,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(z);
#endif
		// std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(u);
#endif
#if defined(RECORD_RESULT)
		float znorm;
		float Ax[DIAG], Axb[DIAG];
		float Axbnorm2;
		float xz[DIAG], zoldz[DIAG], rhozoldz[DIAG];
		float xznorm, rhozoldznorm;
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
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(rhozoldz, rhozoldznorm);
		hist.s_norm[k] = rhozoldznorm;
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
	DISPLAYobj.printvector<float,DIAG>(x);
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



void ADMM_LASSO_FLOAT(float **A, float **At,
					   float **invL,float **invU,
					   float *b, float lambda,
					   float rho, float alpha){

//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// parameters
	float oneminusalpha = 1 - alpha;
	float lambdadivrho = lambda / rho;

	// variables
	Float_Point_Algebra Float_Point_Algebraobj;
#if defined(DEBUG_ITER)
	DISPLAY DISPLAYobj;
#endif
//	float **At;
//	At = (float**) malloc(sizeof(float*)*DIAG);
//	for(int i=0;i<DIAG;i++)
//		At[i] = (float*) malloc(sizeof(float)*DIAG);
//	float **AtA;
//	float **EYE;
//	float **rhoEYE;
//	float **AtAplusrhoeye;
//	float **L, **U;
//	float **invL, **invU;
	float *Atb;
	float *x, *z, *u;
	float  *zold;
	Atb = (float*) malloc(sizeof(float)*DIAG);
	x = (float*) malloc(sizeof(float)*DIAG);
	z = (float*) malloc(sizeof(float)*DIAG);
	u = (float*) malloc(sizeof(float)*DIAG);
	zold = (float*) malloc(sizeof(float)*DIAG);
//	At = (float**) malloc(sizeof(float*)*DIAG);
//	AtA = (float**) malloc(sizeof(float*)*DIAG);
//	EYE = (float**) malloc(sizeof(float*)*DIAG);
//	rhoEYE = (float**) malloc(sizeof(float*)*DIAG);
//	AtAplusrhoeye = (float**) malloc(sizeof(float*)*DIAG);
//	L = (float**) malloc(sizeof(float*)*DIAG);
//	U = (float**) malloc(sizeof(float*)*DIAG);
//	invL = (float**) malloc(sizeof(float*)*DIAG);
//	invU = (float**) malloc(sizeof(float*)*DIAG);
//	for(int i=0;i<DIAG;i++){
//		At[i] = (float*) malloc(sizeof(float)*DIAG);
//		AtA[i] = (float*) malloc(sizeof(float)*DIAG);
//		EYE[i] = (float*) malloc(sizeof(float)*DIAG);
//		rhoEYE[i] = (float*) malloc(sizeof(float)*DIAG);
//		AtAplusrhoeye[i] = (float*) malloc(sizeof(float)*DIAG);
//		L[i] = (float*) malloc(sizeof(float)*DIAG);
//		U[i] = (float*) malloc(sizeof(float)*DIAG);
//		invL[i] = (float*) malloc(sizeof(float)*DIAG);
//		invU[i] = (float*) malloc(sizeof(float)*DIAG);
//	}

//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*b
//	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(A, At);
#if defined(DEBUG_ITER)
	std::cout << "At:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(At);
#endif
	Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(At, b, Atb);
#if defined(DEBUG_ITER)
	std::cout << "Atb:" << std::endl;
	DISPLAYobj.printvector<float,DIAG>(Atb);
#endif
//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// A'*A + rho*speye(n)
//	Float_Point_Algebraobj.MAT_MUL<float,DIAG,DIAG,DIAG>(At, A, AtA);
#if defined(DEBUG_ITER)
	std::cout << "AtA:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(AtA);
#endif
//	Float_Point_Algebraobj.IDENDTITY_MAT<float,DIAG,DIAG>(EYE);
//	Float_Point_Algebraobj.MAT_SCALAR_DOTMUL<float,DIAG,DIAG>(
//			EYE, rho, rhoEYE);
//	Float_Point_Algebraobj.MAT_ADD<float,DIAG,DIAG>(AtA,
//			rhoEYE, AtAplusrhoeye);
#if defined(DEBUG_ITER)
	std::cout << "AtAplusrhoeye:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(AtAplusrhoeye);
#endif
//	 std::cout << __FILE__ << "," << __LINE__ << std::endl;
	// LU
//	Float_Point_Algebraobj.LU_CHOLBANACHROUT<float,DIAG>(
//			AtAplusrhoeye, L, U);
#if defined(DEBUG_ITER)
	std::cout << "L:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(L);
	std::cout << "U:" << std::endl;
	DISPLAYobj.printmatrix<float,DIAG,DIAG>(U);
#endif
	// invers L and U;
//	Float_Point_Algebraobj.MAT_QRINV<float,DIAG>(L, invL);
//	Float_Point_Algebraobj.MAT_TRANS<float,DIAG,DIAG>(invL, invU);
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
	// iteration record and early termination
	float** u_hist = NULL;
	float** x_hist = NULL;
	float** z_hist = NULL;
	float* objval = NULL;
	float* r_norm = NULL;
	float* s_norm = NULL;
	float* eps_pri = NULL;
	float* eps_dual = NULL;
	u_hist = (float**) malloc(sizeof(float*)*DIAG);
	x_hist = (float**) malloc(sizeof(float*)*DIAG);
	z_hist = (float**) malloc(sizeof(float*)*DIAG);
	for(int i=0;i<DIAG;i++){
		u_hist[i] = (float*) malloc(sizeof(float)*MAX_ITER);
		x_hist[i] = (float*) malloc(sizeof(float)*MAX_ITER);
		z_hist[i] = (float*) malloc(sizeof(float)*MAX_ITER);
	}
	objval = (float*) malloc(sizeof(float)*MAX_ITER);
	r_norm = (float*) malloc(sizeof(float)*MAX_ITER);
	s_norm = (float*) malloc(sizeof(float)*MAX_ITER);
	eps_pri = (float*) malloc(sizeof(float)*MAX_ITER);
	eps_dual = (float*) malloc(sizeof(float)*MAX_ITER);

	int k = 0;
    printf("%s\t%s\t%s\t%s\t%s\t%s\n", "iter",
      "r norm", "eps pri", "s norm", "eps dual", "objective");
	for(k=0;k<MAX_ITER;k++){
		float *zminusu=NULL, *rhozminusu=NULL, *q=NULL;
		float *invLq=NULL;
		float *alphax=NULL,*oneminusalphazold=NULL,*x_hat=NULL;
		float *x_hatu=NULL,*x_hatu1=NULL,*x_hatu2=NULL;
		float *x_hatz=NULL;
		zminusu = (float*) malloc(sizeof(float)*DIAG);
		rhozminusu = (float*) malloc(sizeof(float)*DIAG);
		q = (float*) malloc(sizeof(float)*DIAG);
		invLq = (float*) malloc(sizeof(float)*DIAG);
		alphax = (float*) malloc(sizeof(float)*DIAG);
		oneminusalphazold = (float*) malloc(sizeof(float)*DIAG);
		x_hat = (float*) malloc(sizeof(float)*DIAG);
		x_hatu = (float*) malloc(sizeof(float)*DIAG);
		x_hatu1 = (float*) malloc(sizeof(float)*DIAG);
		x_hatu2 = (float*) malloc(sizeof(float)*DIAG);
		x_hatz = (float*) malloc(sizeof(float)*DIAG);
		// q = Atb + rho*(z - u);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(z, u, zminusu);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float,DIAG>(zminusu,
				rho, rhozminusu);
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(Atb,
				rhozminusu, q);
#if defined(DEBUG_ITER)
		std::cout << "q:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(q);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// x = U \ (L \ q);
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invL, q, invLq);
#if defined(DEBUG_ITER)
		std::cout << "invL:" << std::endl;
		DISPLAYobj.printmatrix<float,DIAG,DIAG>(invL);
		std::cout << "invU:" << std::endl;
		DISPLAYobj.printmatrix<float,DIAG,DIAG>(invU);
		std::cout << "invLq:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(invLq);
#endif
		Float_Point_Algebraobj.MAT_VEC_MUL<float,DIAG,DIAG>(
				invU, invLq, x);
#if defined(DEBUG_ITER)
		std::cout << "x:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
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
		DISPLAYobj.printvector<float,DIAG>(x_hat);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// z = shrinkage(x_hat + u, lambda/rho)
		// 			shrinkage(x, kappa):
		// 			z = max( 0, x - kappa ) - max( 0, -x - kappa );
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(x_hat, u, x_hatu);
#if defined(DEBUG_ITER)
		std::cout << "xhatu:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatu);
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
		DISPLAYobj.printvector<float,DIAG>(x_hatu1);
		std::cout << "xhatu2:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatu2);
#endif
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hatu1, x_hatu2, z);
#if defined(DEBUG_ITER)
		std::cout << "z:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(z);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
		// u = u + (x_hat - z);
		Float_Point_Algebraobj.VEC_SUB<float,DIAG>(x_hat, z, x_hatz);
#if defined(DEBUG_ITER)
		std::cout << "x_hatz:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(x_hatz);
#endif
		Float_Point_Algebraobj.VEC_ADD<float,DIAG>(u, x_hatz, u);
#if defined(DEBUG_ITER)
		std::cout << "u:" << std::endl;
		DISPLAYobj.printvector<float,DIAG>(u);
#endif
//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;
#if defined(RECORD_RESULT)
		float znorm;
		float Axbnorm2;
		float xznorm, rhozoldznorm;
		float xnorm;
		float rhounorm;
		float *Ax=NULL, *Axb=NULL;
		float *xz=NULL, *zoldz=NULL, *rhozoldz=NULL;
		float *rhou=NULL;
		Ax = (float*) realloc(Ax,sizeof(float)*DIAG);
		Axb = (float*) realloc(Axb,sizeof(float)*DIAG);
		xz = (float*) realloc(xz,sizeof(float)*DIAG);
		zoldz = (float*) realloc(zoldz,sizeof(float)*DIAG);
		rhozoldz = (float*) realloc(rhozoldz,sizeof(float)*DIAG);
		rhou = (float*) realloc(rhou,sizeof(float)*DIAG);

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.objval(k)  = objective(A, b, lambda, x, z);
		// p = objective(A, b, lambda, x, z)
	    //     p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) )
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(z, znorm);
		Float_Point_Algebraobj.MAT_VEC_MUL<float, DIAG, DIAG>(A, x, Ax);
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(Ax, b, Axb);
		Float_Point_Algebraobj.VEC_NORM2<float, DIAG>(Axb, Axbnorm2);
		objval[k] = 0.5 * Axbnorm2 + lambda * znorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.r_norm(k)  = norm(x - z);
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(x, z, xz);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(xz, xznorm);
		r_norm[k] = xznorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.s_norm(k)  = norm(-rho*(z - zold));
		Float_Point_Algebraobj.VEC_SUB<float, DIAG>(zold, z, zoldz);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float, DIAG>(zoldz, rho, rhozoldz);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(rhozoldz, rhozoldznorm);
		s_norm[k] = rhozoldznorm;

//		 std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(x, xnorm);
		eps_pri[k] = std::sqrt(DIAG)*ABSTOL+RELTOL*(xnorm>=znorm?xnorm:znorm);
		// history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		Float_Point_Algebraobj.VEC_SCALAR_MUL<float, DIAG>(rho, u, rhou);
		Float_Point_Algebraobj.VEC_NORM<float, DIAG>(rhou, rhounorm);
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
	DISPLAYobj.printvector<float,DIAG>(x);
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
