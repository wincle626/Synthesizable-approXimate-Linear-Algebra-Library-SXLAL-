/*
 * admm_lasso.hpp
 *
 *  Created on: 16 Sep 2019
 *      Author: yw106
 */

#ifndef SRC_ADMM_LASSO_HPP_
#define SRC_ADMM_LASSO_HPP_

#include "matplotlibcpp.hpp"
#include "eigen_algebra.hpp"
#include "xfxpt_algebra.hpp"
#include "fpt_algebra.hpp"
#include "softposit_algebra.hpp"
#include "floatx.hpp"

#define TIME_PROFILE
#define PLOT_FIGURE // plot the figure
#define SHOW_FIGURE // show the plot
//#define DEBUG_DATA // print data
//#define DEBUG_ITER // print gradient decent
#define TIME_PROFILE // print and save time profiling
//#define ALWAYS_DELETE // deleting saved file
#define RECORD_RESULT // recording result switch
#define EARLY_TERMINATE // early termination

#define QUIET      0
#define ABSTOL     1e-4
#define RELTOL     1e-2
#define MAX_ITER   1000

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

//void ADMM_LASSO_INTEGER(int A[DIAG][DIAG],
//					  int b[DIAG], int lambda,
//					  int rho, int alpha);

void ADMM_LASSO(std::string path);
#endif /* SRC_ADMM_LASSO_HPP_ */
