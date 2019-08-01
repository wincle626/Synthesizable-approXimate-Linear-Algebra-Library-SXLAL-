/*
 * xilinx_fixed_point.hpp
 *
 *  Created on: 28 Jun 2019
 *      Author: yunwu
 */

#ifndef SRC_XFXPT_HPP_
#define SRC_XFXPT_HPP_

#include "ap_int.h"
#include "ap_fixed.h"

#define INPUT_BIT_WIDTH 48
#define INPUT_INTE_WIDTH 12
#define INPUT_FRAC_WIDTH 12
#define OUTPUT_BIT_WIDTH 7
#define OUTPUT_INTE_WIDTH 3
#define OUTPUT_FRAC_WIDTH 4

// W  Word length in bits

// I  The number of bits used to represent the integer value
//	  (the number of bits above the decimal point)

// Q Quantization mode dictates the behavior when greater
//   precision is generated than can be defined by smallest
//   fractional bit in the variable used to store the result:
// AP_RND          Rounding to plus infinity
// AP_RND_ZERO     Rounding to zero
// AP_RND_MIN_INF  Rounding to minus infinity
// AP_RND_INF      Rounding to infinity
// AP_RND_CONV     Convergent rounding
// AP_TRN          Truncation to minus infinity
// AP_TRN_ZERO     Truncation to zero (default)

// O Overflow mode dictates the behavior when more bits
//   are generated than the variable to store the result
//   contains:
// AP_SAT       Saturation
// AP_SAT_ZERO  Saturation to zero
// AP_SAT_SYM   Symmetrical saturation
// AP_WRAP      Wrap around (default)
// AP_WRAP_SM   Sign magnitude wrap around

// N  The number of saturation bits in wrap modes.

// fixed point format: e.g. ap_[u]fixed<int W, int I, ap_q_mode Q, ap_o_mode O,ap_sat_bits N>

typedef ap_fixed<INPUT_BIT_WIDTH,INPUT_INTE_WIDTH> DATA_IN_T;
typedef ap_fixed<OUTPUT_BIT_WIDTH,OUTPUT_INTE_WIDTH> DATA_OUT_T;

#endif /* SRC_XFXPT_HPP_ */
