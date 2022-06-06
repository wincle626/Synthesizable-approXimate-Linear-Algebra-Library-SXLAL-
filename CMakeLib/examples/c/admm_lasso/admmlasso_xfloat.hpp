/*
 * admmlasso_xfloat.hpp
 *
 *  Created on: 11 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_ADMMLASSO_XFLOAT_HPP_
#define SRC_ADMMLASSO_XFLOAT_HPP_

#include "admm_lasso.hpp"

// 1 bit reserved for sign
// floatx<exponent, significand>  bitwidth = significand + exponent + 1

using fptx_admmlasso = flx::floatx<8, 19>;

void ADMM_LASSO_XFPT(fptx_admmlasso A[DIAG][DIAG], fptx_admmlasso At[DIAG][DIAG],
					   fptx_admmlasso invL[DIAG][DIAG],fptx_admmlasso invU[DIAG][DIAG],
					   fptx_admmlasso b[DIAG], fptx_admmlasso lambda,
					   fptx_admmlasso rho, fptx_admmlasso alpha);

void ADMM_LASSO_XFPT(fptx_admmlasso **A, fptx_admmlasso **At,
					   fptx_admmlasso **invL,fptx_admmlasso **invU,
					   fptx_admmlasso *b, fptx_admmlasso lambda,
					   fptx_admmlasso rho, fptx_admmlasso alpha);


#endif /* SRC_ADMMLASSO_XFLOAT_HPP_ */
