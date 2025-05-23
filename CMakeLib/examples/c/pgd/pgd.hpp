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
#include "softposit_algebra.hpp"
#include "floatx.hpp"

#define ITER_MAX 10000 // Iteration constraint
#define EPS_STOP 0.0001 // Early termination constraint
#define PGD_INT_SCALE 1 // Early termination constraint
#define EPS_STOP_SCALE (int)EPS_STOP*PGD_INT_SCALE*PGD_INT_SCALE // Early termination constraint
#define BOX_CONST 10*PGD_INT_SCALE*PGD_INT_SCALE // Box constraint on the variable
#define ERR_STD 0 // Error standard deviation when box projection

//#define DEBUG_DATA // print data
//#define DEBUG_LAFUNC // print comparison between Eigen3 and General Double
//#define DEBUG_ITER // print gradient decent

#define TIME_PROFILE // print and save time profiling
#define PLOT_FIGURE // plot the figure
#define SHOW_FIGURE // show the plot
//#define ALWAYS_DELETE // deleting saved file
#define RECORD_RESULT // recording result switch
//#define REACTIVE_ITERATION // reactive iteration switch

void PROXIMAL_GRADIENT_DECENT(std::string path);

#endif /* SRC_PGD_HPP_ */
