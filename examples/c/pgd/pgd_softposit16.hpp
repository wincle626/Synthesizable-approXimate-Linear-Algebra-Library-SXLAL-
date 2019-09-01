/*
 * pgd_softposit16.hpp
 *
 *  Created on: 29 Aug 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_SOFTPOSIT16_HPP_
#define SRC_PGD_SOFTPOSIT16_HPP_


#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_SPOSIT16(posit16_t Amatrix_c[DIAG][DIAG],
									  posit16_t bvector_c[DIAG],
									  posit16_t L_c);

#endif /* SRC_PGD_SOFTPOSIT16_HPP_ */
