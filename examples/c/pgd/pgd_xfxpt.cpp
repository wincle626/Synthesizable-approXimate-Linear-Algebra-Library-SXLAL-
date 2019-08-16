/*
 * proximal_gradient_decent_xilinxfxpt.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: yunwu
 */

#include "pgd_xfxpt.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_XFXPT(DATA_IN_T Amatrix_c[DIAG][DIAG],
									DATA_IN_T bvector_c[DIAG],
									DATA_IN_T factor){

#ifdef XILINX_FIXED_PRECISION
	Xilinx_Fixed_Point_Algebra Xilinx_Fixed_Point_Algebra_obj;
	float error = 0;
	float error_std = ERR_STD;

	std::string clockname = "timeprofile_xfxpt.txt";
	std::string xkname = "xk_xfxpt.dat";
	std::string errorrecordname = "error_record_xfxpt.dat";
	std::string errorhistname = "error_hist_xfxpt.dat";
	std::string figurename = "ProximalGradientDecent_xfxpt.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	DATA_IN_T x_k_vec[DIAG];
   	DATA_IN_T x_k_plus1_vec[DIAG];
   	DATA_IN_T grad_g[DIAG];
   	DATA_IN_T ones_vec[DIAG];
   	std::vector<float> error_record(ITER_MAX, 0);
   	std::vector<float> error_hist(ITER_MAX, 0);
   	std::vector<std::vector<float>> x_k_record;

   	Xilinx_Fixed_Point_Algebra_obj.ONES_VEC<DATA_IN_T, DIAG>( ones_vec );
   	Xilinx_Fixed_Point_Algebra_obj.ZEROS_VEC<DATA_IN_T, DIAG>( x_k_vec );

	int k=0;
#ifdef REACTIVE_ITERATION
	float reactive_factor = 0.1;
#endif
	while( k<ITER_MAX ){

		//std::cout << __FILE__ << __LINE__ << std::endl << std::endl;
		// Save previous point
		Xilinx_Fixed_Point_Algebra_obj.VEC_EQ<DATA_IN_T, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		DATA_IN_T tmp_mul_vec1[DIAG];
		Xilinx_Fixed_Point_Algebra_obj.MAT_VEC_MUL<DATA_IN_T, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "tmp_mul_vec1["<< i << "]: " << tmp_mul_vec1[i] << std::endl << std::endl;
#endif
		Xilinx_Fixed_Point_Algebra_obj.VEC_ADD<DATA_IN_T, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "grad_g["<< i << "]: " << grad_g[i] << std::endl << std::endl;
#endif

		// new decent point
		DATA_IN_T new_point[DIAG];
		DATA_IN_T tmp_mul_vec2[DIAG];
		Xilinx_Fixed_Point_Algebra_obj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>( grad_g, factor,
									 tmp_mul_vec2 );
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "tmp_mul_vec2["<< i << "]: " << tmp_mul_vec2[i] << std::endl << std::endl;
#endif
		Xilinx_Fixed_Point_Algebra_obj.VEC_SUB<DATA_IN_T, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "new_point["<< i << "]: " << new_point[i] << std::endl << std::endl;
#endif

		// Proximal projection
		DATA_IN_T rndnoise[DIAG];
		DATA_IN_T a_ones_vec[DIAG];
		DATA_IN_T minus_a_ones_vec[DIAG];
		Xilinx_Fixed_Point_Algebra_obj.ONES_VEC<DATA_IN_T, DIAG>( a_ones_vec );
		Xilinx_Fixed_Point_Algebra_obj.ONES_VEC<DATA_IN_T, DIAG>( minus_a_ones_vec );
		DATA_IN_T tmp_min[DIAG];
		DATA_IN_T tmp_max[DIAG];
		DATA_IN_T tmp_mul_vec3[DIAG];
		Xilinx_Fixed_Point_Algebra_obj.RND_VEC<DATA_IN_T, DIAG>( rndnoise );
		Xilinx_Fixed_Point_Algebra_obj.VEC_SCALAR_MIN<DATA_IN_T, DIAG>(new_point, BOX_CONST, tmp_min);
		Xilinx_Fixed_Point_Algebra_obj.VEC_SCALAR_MAX<DATA_IN_T, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Xilinx_Fixed_Point_Algebra_obj.VEC_SCALAR_MUL<DATA_IN_T, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Xilinx_Fixed_Point_Algebra_obj.VEC_ADD<DATA_IN_T, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "tmp_mul_vec3["<< i << "]: " << tmp_mul_vec3[i] << std::endl << std::endl;
		for(int i=0;i<DIAG;i++)
			std::cout << "x_k_vec["<< i << "]: " << x_k_vec[i] << std::endl << std::endl;
		for(int i=0;i<DIAG;i++)
			std::cout << "x_k_plus1_vec["<< i << "]: " << x_k_plus1_vec[i] << std::endl << std::endl;
#endif

		// check early termination constraint
		DATA_IN_T norm1, norm2;
		DATA_IN_T tmp_sub_vec[DIAG];
		Xilinx_Fixed_Point_Algebra_obj.VEC_SUB<DATA_IN_T, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
#ifdef DEBUG_ITER
		for(int i=0;i<DIAG;i++)
			std::cout << "tmp_sub_vec["<< i << "]: " << tmp_sub_vec[i] << std::endl << std::endl;
#endif
		Xilinx_Fixed_Point_Algebra_obj.VEC_NORM_FXSQRT<DATA_IN_T,INPUT_BIT_WIDTH, INPUT_INTE_WIDTH,
				INPUT_BIT_WIDTH, INPUT_INTE_WIDTH, DIAG>(tmp_sub_vec, norm1);
		Xilinx_Fixed_Point_Algebra_obj.VEC_NORM_FXSQRT<DATA_IN_T,INPUT_BIT_WIDTH, INPUT_INTE_WIDTH,
				INPUT_BIT_WIDTH, INPUT_INTE_WIDTH, DIAG>(x_k_vec, norm2);
//		Xilinx_Fixed_Point_Algebra_obj.VEC_NORM<DATA_IN_T, DIAG>(tmp_sub_vec, norm1);
//		Xilinx_Fixed_Point_Algebra_obj.VEC_NORM<DATA_IN_T, DIAG>(x_k_vec, norm2);
#ifdef DEBUG_ITER
		std::cout << "norm1: " << norm1 << std::endl;
		std::cout << "norm1: " << norm1 << std::endl << std::endl;
#endif

		// record relative error
#ifdef REACTIVE_ITERATION
		if((float)norm1>0 && (float)norm2>0){
			error = norm1;
			error_hist.at(k) = (float)(norm1 / norm2);
			error_record.at(k) = error;
		}else{
			error = error_record.at(k-1)*1.1;
			error_hist.at(k) = error_hist.at(k-1);
			error_record.at(k) = error;
			Xilinx_Fixed_Point_Algebra_obj.VEC_SCALAR_SUB<DATA_IN_T,  DIAG>(x_k_vec, reactive_factor, x_k_vec);
		}
#else
		error = norm1;
		error_hist.at(k) = (float)(norm1 / norm2);
		error_record.at(k) = error;
#endif// endif REACTIVE_ITERATION
#ifdef RECORD_RESULT
		std::vector<float> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((float)x_k_vec[i]);
		x_k_record.push_back(xkvec);
#endif// endif RECORD_RESULT

		//std::cout << "error: " <<error << std::endl;
		if( (float)error< EPS_STOP ){
			if((float)norm1==0 && k==0){
				std::cout << "something is wrong as no sccuessful iteration. " << std::endl << std::endl;
				std::exit(0);
			}
			std::cout << "breaking error: " <<error << std::endl;
			break;
		}

#ifdef DEBUG_ITER
		std::cout << k << "-iteration" << std::endl
				  << "error: " << error << std::endl<< std::endl;
#endif

		k++;

	}
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to gradient decent"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to gradient decent"
			  << std::endl << std::endl;
#endif
	///////////////////////////////////////////////////////////



#ifdef RECORD_RESULT
	/////////////////// Print Iteration Number ////////////////
	std::cout << "Number of iterations = " << k
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-1)
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-1) << std::endl;
	///////////////////////////////////////////////////////////


	//////////////////////// Plot Figures /////////////////////
	// save data
	std::ofstream resultfile;
	std::ofstream resultfile1;
	std::ofstream resultfile2;
	resultfile.open(errorhistname);
	resultfile1.open(errorrecordname);
	resultfile2.open(xkname);
	resultfile << "# iteration relative_error \n";
	resultfile1 << "# iteration error \n";
	resultfile2 << "# optimizing x_k vector \n";
	std::vector<double> x, y;
	for(int i=0; i<ITER_MAX; i++){
		if( error_hist.at(i)>0 ){
			x.push_back(i);
			y.push_back(error_hist.at(i));
			resultfile << "  " << i << "," << error_hist.at(i) << " \n";
			resultfile1 << "  " << i << "," << error_record.at(i) << " \n";
			std::vector<float> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << xkvec.at(j) << ",";
			resultfile2 << xkvec.at(DIAG-1) << "\n";
		}
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();

#ifdef PLOT_FIGURE
	// Matplotlib plotting
	plt::named_semilogy( "simple example", x, y );
	plt::title("Proximal Gradient Decent");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Relative Error");
	plt::legend();
	plt::save(figurename);
#ifdef SHOW_FIGURE
   	plt::show();
#endif // endif PLOT_FIGURE
#endif // endif PLOT_FIGURE
#endif // endif RECORD_RESULT
#endif // endif XILINX_FIXED_PRECISION
}

