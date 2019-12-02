/*
 * pgd_xfloat.hpp
 *
 *  Created on: Aug 15, 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_XFPTHPP_
#define SRC_PGD_XFPTHPP_


#include "pgd.hpp"

//using fptx1 = flx::floatx<8, 32>;
using fptx2 = flx::floatx<5, 8>;

//void PROXIMAL_GRADIENT_DECENT_XFLOAT1(fptx1 Amatrix_c[DIAG][DIAG],
//									 fptx1 bvector_c[DIAG],
//									 fptx1 L_c);

void PROXIMAL_GRADIENT_DECENT_XFLOAT2(fptx2 Amatrix_c[DIAG][DIAG],
									 fptx2 bvector_c[DIAG],
									 fptx2 L_c);

void PROXIMAL_GRADIENT_DECENT_XFLOAT2(fptx2 **Amatrix_c,
									 fptx2 *bvector_c,
									 fptx2 L_c);

#endif /* SRC_PGD_XFPTHPP_ */
