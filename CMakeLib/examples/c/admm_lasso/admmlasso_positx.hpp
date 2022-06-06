/*
 * admmlasso_positx.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_POSITX_HPP_
#define SRC_ADMMLASSO_POSITX_HPP_

#include "admm_lasso.hpp"



void ADMM_LASSO_POSITX(posit_2_t A[DIAG][DIAG], posit_2_t At[DIAG][DIAG],
		   	   	   	   posit_2_t invL[DIAG][DIAG],posit_2_t invU[DIAG][DIAG],
					   posit_2_t b[DIAG], posit_2_t lambda,
					   posit_2_t rho, posit_2_t alpha);

void ADMM_LASSO_POSITX(posit_2_t **A, posit_2_t **At,
					   posit_2_t **invL,posit_2_t **invU,
					   posit_2_t *b, posit_2_t lambda,
					   posit_2_t rho, posit_2_t alpha);


#endif /* SRC_ADMMLASSO_POSITX_HPP_ */
