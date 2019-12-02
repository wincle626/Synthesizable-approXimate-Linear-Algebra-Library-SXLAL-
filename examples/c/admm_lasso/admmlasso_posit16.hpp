/*
 * admmlasso_posit16.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_POSIT16_HPP_
#define SRC_ADMMLASSO_POSIT16_HPP_

#include "admm_lasso.hpp"

void ADMM_LASSO_POSIT16(posit16_t A[DIAG][DIAG], posit16_t At[DIAG][DIAG],
		   	   	   	   posit16_t invL[DIAG][DIAG],posit16_t invU[DIAG][DIAG],
					   posit16_t b[DIAG], posit16_t lambda,
					   posit16_t rho, posit16_t alpha);

void ADMM_LASSO_POSIT16(posit16_t **A, posit16_t **At,
					   posit16_t **invL,posit16_t **invU,
					   posit16_t *b, posit16_t lambda,
					   posit16_t rho, posit16_t alpha);



#endif /* SRC_ADMMLASSO_POSIT16_HPP_ */
