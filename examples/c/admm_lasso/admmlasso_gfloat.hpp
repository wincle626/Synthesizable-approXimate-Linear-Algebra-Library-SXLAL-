/*
 * admmlasso_gfloat.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_GFLOAT_HPP_
#define SRC_ADMMLASSO_GFLOAT_HPP_

#include "admm_lasso.hpp"

void ADMM_LASSO_FLOAT(float A[DIAG][DIAG], float At[DIAG][DIAG],
		   	   	   	   float invL[DIAG][DIAG],float invU[DIAG][DIAG],
					   float b[DIAG], float lambda,
					   float rho, float alpha);

void ADMM_LASSO_FLOAT(float **A, float **At,
					   float **invL,float **invU,
					   float *b, float lambda,
					   float rho, float alpha);



#endif /* SRC_ADMMLASSO_GFLOAT_HPP_ */
