/*
 * pgd_gfloat.hpp
 *
 *  Created on: 18 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_GFLOAT_HPP_
#define SRC_PGD_GFLOAT_HPP_


#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_GFLOAT(float Amatrix_c[DIAG][DIAG],
									 float bvector_c[DIAG],
									 float L_c);


#endif /* SRC_PGD_GFLOAT_HPP_ */
