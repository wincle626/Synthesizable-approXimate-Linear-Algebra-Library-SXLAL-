/*
 * xilinx_fixed_point_algebra.hpp
 *
 *  Created on: 28 Jun 2019
 *      Author: yunwu
 */

#ifndef SRC_XFXPT_ALGEBRA_HPP_
#define SRC_XFXPT_ALGEBRA_HPP_

#include "fpt_algebra.hpp"
#include "xfxpt.hpp"

class Xilinx_Fixed_Point_Algebra : public Float_Point_Algebra{

public:
	template <int W2, int IW2, int W1, int IW1>
	void FXP_SQRT(ap_ufixed<W1,IW1>& in_val, ap_ufixed<W2,IW2>& result)
	{
	   enum { QW = (IW1+1)/2 + (W2-IW2) + 1 }; // derive max root width
	   enum { SCALE = (W2 - W1) - (IW2 - (IW1+1)/2) }; // scale (shift) to adj initial remainer value
	   enum { ROOT_PREC = QW - (IW1 % 2) };
	   assert((IW1+1)/2 <= IW2); // Check that output format can accommodate full result

	   ap_uint<QW> q      = 0;   // partial sqrt
	   ap_uint<QW> q_star = 0;   // diminished partial sqrt
	   ap_int<QW+2> s; // scaled remainder initialized to extracted input bits
	   if (SCALE >= 0)
	      s = in_val.range(W1-1,0) << (SCALE);
	   else
	      s = ((in_val.range(W1-1,0) >> (0 - (SCALE + 1))) + 1) >> 1;

	   // Non-restoring square-root algorithm
	   for (int i = 0; i <= ROOT_PREC; i++) {
	      if (s >= 0) {
	         s = 2 * s - (((ap_int<QW+2>(q) << 2) | 1) << (ROOT_PREC - i));
	         q_star = q << 1;
	         q = (q << 1) | 1;
	      } else {
	         s = 2 * s + (((ap_int<QW+2>(q_star) << 2) | 3) << (ROOT_PREC - i));
	         q = (q_star << 1) | 1;
	         q_star <<= 1;
	      }
	   }
	   // Round result by "extra iteration" method
	   if (s > 0)
	      q = q + 1;
	   // Truncate excess bit and assign to output format
	   result.range(W2-1,0) = ap_uint<W2>(q >> 1);
	}

	// Vector norm
	template<int W2, int IW2, int W1, int IW1, int M>
	void VEC_NORM_FXSQRT( ap_fixed<W1,IW1> V1[M], ap_fixed<W2,IW2> &S ){
		ap_fixed<W2,IW2> S1 = 0;
		for( int i=0; i<M; i++ ){
			S1 += V1[i] * V1[i];
		}
		ap_ufixed<W1,IW1> tmp1 = S1;
		ap_ufixed<W2,IW2> tmp2;
		this->FXP_SQRT<W2, IW2, W1, IW1>(tmp1, tmp2);
		S = tmp2;
	}

private:

protected:

};

#endif /* SRC_XFXPT_ALGEBRA_HPP_ */
