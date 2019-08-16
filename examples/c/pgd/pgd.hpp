/*
 *	This is a simple example of how to do
 *	proximal gradient decent algorithm in
 *	C++ with Eigen library and plot with
 *	Matplotlib C++ API
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_PGD_HPP_
#define SRC_PGD_HPP_


#include "eigen_algebra.hpp"
#include "fpt_algebra.hpp"
#include "matplotlibcpp.hpp"
#include "xfxpt_algebra.hpp"
#include "floatx.hpp"

#define ITER_MAX 100000 // Iteration constraint
#define EPS_STOP 0.00001 // Early termination constraint
#define BOX_CONST 10 // Box constraint on the variable
#define ERR_STD 0 // Error standard deviation when box projection

//#define EIGEN_DOUBLE_PRECISION // eigen double precision switch
//#define EIGEN_FLOAT_PRECISION // eigen float precision switch
//#define EIGEN_INTEGER_PRECISION // eigen integer precision switch
//#define GENERAL_DOUBLE_PRECISION // general double precision switch
//#define GENERAL_FLOAT_PRECISION // general float precision switch
#define COMSTOM_FLOAT_PRECISION // general float precision switch
//#define GENERAL_INTEGER_PRECISION // general integer precision switch
//#define XILINX_FIXED_PRECISION // fixed point precision swithc

//#define DEBUG_DATA // print data
//#define DEBUG_LAFUNC // print comparison between Eigen3 and General Double
//#define DEBUG_ITER // print gradient decent

#define TIME_PROFILE // print and save time profiling
#define PLOT_FIGURE // plot the figure
//#define SHOW_FIGURE // show the plot
//#define ALWAYS_DELETE // deleting saved file
#define RECORD_RESULT // recording result switch
//#define REACTIVE_ITERATION // reactive iteration switch

void PROXIMAL_GRADIENT_DECENT();

#endif /* SRC_PGD_HPP_ */
