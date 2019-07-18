/*
 *	This is a simple example of how to call
 *	Xilinx fixed-point library to do square
 *	root operation
 *	Author: Yun Wu
 *	Created by: 2019-06-28
 *	Copyright @ Yun Wu
 *
 */

#include "xfxpt_sqrt.hpp"


void XILINX_FXP_SQRT()
{
	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebra_obj;
	in_data_t in_val;
	out_data_t result;
	in_val = 3;
	Xilinx_Fixed_Point_Algebra_obj.FXP_SQRT<OUT_BW,OUT_IW,IN_BW,IN_IW>(in_val, result);
	std::cout << "input: " << in_val << std::endl
			  << "output: " << result << std::endl
			  << std::endl << std::endl;
}
