/*
 * pgd_gdouble.hpp
 *
 *  Created on: 18 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_GDOUBLE_HPP_
#define SRC_PGD_GDOUBLE_HPP_

#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_GDOUBLE(double Amatrix_c[DIAG][DIAG],
									  double bvector_c[DIAG],
									  double L_c);


#endif /* SRC_PGD_GDOUBLE_HPP_ */
