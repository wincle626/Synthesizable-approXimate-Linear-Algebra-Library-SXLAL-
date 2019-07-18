/*
 * pgd_eigenf.hpp
 *
 *  Created on: 18 Jul 2019
 *      Author: yunwu
 */

#ifndef SRC_PGD_EIGENF_HPP_
#define SRC_PGD_EIGENF_HPP_


#include "pgd_test.hpp"

void PROXIMAL_GRADIENT_DECENT_EIGENF(Eigen::MatrixXd Amatrix,
									 Eigen::VectorXd bvector,
									 double L);


#endif /* SRC_PGD_EIGENF_HPP_ */
