/*
 * admmlasso_gdouble.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_GDOUBLE_HPP_
#define SRC_ADMMLASSO_GDOUBLE_HPP_

#include "admm_lasso.hpp"

void ADMM_LASSO_DOUBLE(double A[DIAG][DIAG], double At[DIAG][DIAG],
		   	   	   	   double invL[DIAG][DIAG],double invU[DIAG][DIAG],
					   double b[DIAG], double lambda,
					   double rho, double alpha);

void ADMM_LASSO_DOUBLE(double **A, double **At,
					   double **invL,double **invU,
					   double *b, double lambda,
					   double rho, double alpha);

void ADMM_LASSO_DOUBLE(double **A, double *b, double lambda,
					   double rho, double alpha);

void ADMM_LASSO_DOUBLE(double **A, double *b, double lambda);



#endif /* SRC_ADMMLASSO_GDOUBLE_HPP_ */
