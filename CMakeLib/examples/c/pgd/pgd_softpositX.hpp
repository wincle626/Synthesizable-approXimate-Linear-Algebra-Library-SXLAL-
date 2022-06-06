/*
 * pgd_softpositX.hpp
 *
 *  Created on: 4 Nov 2019
 *      Author: yw106
 */

#ifndef SRC_PGD_SOFTPOSITX_HPP_
#define SRC_PGD_SOFTPOSITX_HPP_


#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_SPOSITX(posit_2_t Amatrix_c[DIAG][DIAG],
									  posit_2_t bvector_c[DIAG],
									  posit_2_t L_c);

void PROXIMAL_GRADIENT_DECENT_SPOSITX(posit_2_t **Amatrix_c,
										posit_2_t *bvector_c,
										posit_2_t L_c);



#endif /* SRC_PGD_SOFTPOSITX_HPP_ */
