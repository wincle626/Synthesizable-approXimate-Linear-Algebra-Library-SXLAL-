/*
 * admmlasso_xfxpt.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_XFXPT_HPP_
#define SRC_ADMMLASSO_XFXPT_HPP_

#include "admm_lasso.hpp"

void ADMM_LASSO_FXPT(DATA_IN_T A[DIAG][DIAG], DATA_IN_T At[DIAG][DIAG],
					   DATA_IN_T invL[DIAG][DIAG],DATA_IN_T invU[DIAG][DIAG],
					   DATA_IN_T b[DIAG], DATA_IN_T lambda,
					   DATA_IN_T rho, DATA_IN_T alpha);

void ADMM_LASSO_FXPT(DATA_IN_T **A, DATA_IN_T **At,
					   DATA_IN_T **invL,DATA_IN_T **invU,
					   DATA_IN_T *b, DATA_IN_T lambda,
					   DATA_IN_T rho, DATA_IN_T alpha);



#endif /* SRC_ADMMLASSO_XFXPT_HPP_ */
