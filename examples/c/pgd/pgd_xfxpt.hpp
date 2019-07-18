/*
 * pgd_xfxpt.hpp
 *
 *  Created on: 18 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_XFXPT_HPP_
#define SRC_PGD_XFXPT_HPP_


#include "pgd_test.hpp"

void PROXIMAL_GRADIENT_DECENT_XFXPT(DATA_IN_T Amatrix_c[DIAG][DIAG],
									DATA_IN_T bvector_c[DIAG],
									DATA_IN_T factor);


#endif /* SRC_PGD_XFXPT_HPP_ */
