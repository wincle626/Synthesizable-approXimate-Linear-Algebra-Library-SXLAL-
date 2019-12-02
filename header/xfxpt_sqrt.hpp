/*
 * xilinx_fixed_point_sqrt.hpp
 *
 *  Created on: 28 Jun 2019
 *      Author: yunwu
 */

#ifndef SRC_XFXPT_SQRT_HPP_
#define SRC_XFXPT_SQRT_HPP_

#include "xfxpt.hpp"
#include "xfxpt_algebra.hpp"

#define IN_BW   24
#define IN_IW    8
#define OUT_BW  28
#define OUT_IW   4 // ((IN_IW + 1) / 2)

// typedefs for top-level input and output fixed-point formats
typedef ap_ufixed<IN_BW,IN_IW>   in_data_t;
typedef ap_ufixed<OUT_BW,OUT_IW> out_data_t;

// Top level wrapper function - calls the core template function w/ above types
void XILINX_FXP_SQRT();


#endif /* SRC_XFXPT_SQRT_HPP_ */
