/*
 * admm_lasso.hpp
 *
 *  Created on: 16 Sep 2019
 *      Author: yw106
 */

#ifndef SRC_ADMM_LASSO_HPP_
#define SRC_ADMM_LASSO_HPP_

#include "eigen_algebra.hpp"
#include "xfxpt_algebra.hpp"
#include "fpt_algebra.hpp"
#include "softposit_algebra.hpp"
#include "floatx.hpp"

#define GENERAL_DOUBLE_PRECISION // general double precision switch
//#define GENERAL_FLOAT_PRECISION // general float precision switch
//#define COMSTOM_FLOAT_PRECISION // general float precision switch
//#define XILINX_FIXED_PRECISION // fixed point precision swithc
//#define SOFT_POSIT_PRECISION // SoftPosit precision

//#define TIME_PROFILE
//#define DEBUG_DATA // print data
//#define DEBUG_ITER // print gradient decent
//#define TIME_PROFILE // print and save time profiling
//#define ALWAYS_DELETE // deleting saved file
#define RECORD_RESULT // recording result switch
#define EARLY_TERMINATE // early termination

#define DIM_SIZE1 4

using fptx_admmlasso = flx::floatx<8, 32>;

void ADMM_LASSO_DOUBLE(double A[DIM_SIZE1][DIM_SIZE1],
					   double b[DIM_SIZE1], double lambda,
					   double rho, double alpha);

void ADMM_LASSO_FLOAT(float A[DIM_SIZE1][DIM_SIZE1],
					  float b[DIM_SIZE1], float lambda,
					  float rho, float alpha);

void ADMM_LASSO_FXPT(DATA_IN_T A[DIM_SIZE1][DIM_SIZE1],
					 DATA_IN_T b[DIM_SIZE1], DATA_IN_T lambda,
					 DATA_IN_T rho, float alpha);

void ADMM_LASSO_XFPT(fptx_admmlasso A[DIM_SIZE1][DIM_SIZE1],
					 fptx_admmlasso b[DIM_SIZE1], fptx_admmlasso lambda,
					 fptx_admmlasso rho, float alpha);

void ADMM_LASSO_POSIT8(posit8_t A[DIM_SIZE1][DIM_SIZE1],
					   posit8_t b[DIM_SIZE1], posit8_t lambda,
					   posit8_t rho, posit8_t alpha);

void ADMM_LASSO_POSIT16(posit16_t A[DIM_SIZE1][DIM_SIZE1],
						posit16_t b[DIM_SIZE1], posit16_t lambda,
						posit16_t rho, posit16_t alpha);

void ADMM_LASSO_POSIT32(posit32_t A[DIM_SIZE1][DIM_SIZE1],
						posit32_t b[DIM_SIZE1], posit32_t lambda,
						posit32_t rho, posit32_t alpha);

void ADMM_LASSO();
#endif /* SRC_ADMM_LASSO_HPP_ */
