/*
 * pgd_softposit8.hpp
 *
 *  Created on: 29 Aug 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_SOFTPOSIT8_HPP_
#define SRC_PGD_SOFTPOSIT8_HPP_


#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_SPOSIT8(posit8_t Amatrix_c[DIAG][DIAG],
									  posit8_t bvector_c[DIAG],
									  posit8_t L_c);

#endif /* SRC_PGD_SOFTPOSIT8_HPP_ */
