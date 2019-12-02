/*
 * admmlasso_posit8.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_POSIT8_HPP_
#define SRC_ADMMLASSO_POSIT8_HPP_

#include "admm_lasso.hpp"


void ADMM_LASSO_POSIT8(posit8_t A[DIAG][DIAG], posit8_t At[DIAG][DIAG],
		   	   	   	   posit8_t invL[DIAG][DIAG],posit8_t invU[DIAG][DIAG],
					   posit8_t b[DIAG], posit8_t lambda,
					   posit8_t rho, posit8_t alpha);

void ADMM_LASSO_POSIT8(posit8_t **A, posit8_t **At,
					   posit8_t **invL,posit8_t **invU,
					   posit8_t *b, posit8_t lambda,
					   posit8_t rho, posit8_t alpha);



#endif /* SRC_ADMMLASSO_POSIT8_HPP_ */
