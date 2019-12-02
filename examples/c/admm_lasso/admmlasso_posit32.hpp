/*
 * admmlasso_posit32.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_POSIT32_HPP_
#define SRC_ADMMLASSO_POSIT32_HPP_

#include "admm_lasso.hpp"



void ADMM_LASSO_POSIT32(posit32_t A[DIAG][DIAG], posit32_t At[DIAG][DIAG],
		   	   	   	   posit32_t invL[DIAG][DIAG],posit32_t invU[DIAG][DIAG],
					   posit32_t b[DIAG], posit32_t lambda,
					   posit32_t rho, posit32_t alpha);

void ADMM_LASSO_POSIT32(posit32_t **A, posit32_t **At,
					   posit32_t **invL,posit32_t **invU,
					   posit32_t *b, posit32_t lambda,
					   posit32_t rho, posit32_t alpha);

#endif /* SRC_ADMMLASSO_POSIT32_HPP_ */
