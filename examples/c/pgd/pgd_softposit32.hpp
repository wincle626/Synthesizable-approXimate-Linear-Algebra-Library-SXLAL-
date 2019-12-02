/*
 * pgd_softposit32.hpp
 *
 *  Created on: Aug 22, 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_SOFTPOSIT32_HPP_
#define SRC_PGD_SOFTPOSIT32_HPP_

#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_SPOSIT32(posit32_t Amatrix_c[DIAG][DIAG],
									  posit32_t bvector_c[DIAG],
									  posit32_t L_c);

void PROXIMAL_GRADIENT_DECENT_SPOSIT32(posit32_t **Amatrix_c,
									  posit32_t *bvector_c,
									  posit32_t L_c);

#endif /* SRC_PGD_SOFTPOSIT32_HPP_ */
