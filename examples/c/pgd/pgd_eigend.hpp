/*
 * pgd_eigend.hpp
 *
 *  Created on: 18 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_EIGEND_HPP_
#define SRC_PGD_EIGEND_HPP_


#include "pgd.hpp"

void PROXIMAL_GRADIENT_DECENT_EIGEND(Eigen::MatrixXd Amatrix,
									 Eigen::VectorXd bvector,
									 double L);


#endif /* SRC_PGD_EIGEND_HPP_ */
