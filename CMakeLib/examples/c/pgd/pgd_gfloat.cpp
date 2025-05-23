/*
 * proximal_gradient_decent_float.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: yunwu
 */



#include "pgd_gfloat.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_GFLOAT(float Amatrix_c[DIAG][DIAG],
									 float bvector_c[DIAG],
									 float L_c){

#ifdef GENERAL_FLOAT_PRECISION
	Float_Point_Algebra Float_Point_Algebra_obj;
	float error = 0;
	float error_std = ERR_STD;

	std::string clockname = "timeprofile_gfloat.txt";
	std::string xkname = "xk_gfloat.dat";
	std::string errorrecordname = "error_record_gfloat.dat";
	std::string errorhistname = "error_hist_gfloat.dat";
	std::string figurename = "ProximalGradientDecent_gfloat.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	float x_k_vec[DIAG];
   	float x_k_plus1_vec[DIAG];
   	float grad_g[DIAG];
   	float ones_vec[DIAG];
   	std::vector<float> error_hist(ITER_MAX, 0);
   	std::vector<float> error_record(ITER_MAX, 0);
   	std::vector<std::vector<float>> x_k_record;

   	Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( ones_vec );
   	Float_Point_Algebra_obj.ZEROS_VEC<float, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<float, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		float tmp_mul_vec1[DIAG];
		Float_Point_Algebra_obj.MAT_VEC_MUL<float, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
		Float_Point_Algebra_obj.VEC_ADD<float, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		float new_point[DIAG];
		float tmp_mul_vec2[DIAG];
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<float, DIAG>( grad_g, 1/L_c,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<float, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		float rndnoise[DIAG];
		float a_ones_vec[DIAG];
		float minus_a_ones_vec[DIAG];
		Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( a_ones_vec );
		Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( minus_a_ones_vec );
		float tmp_min[DIAG];
		float tmp_max[DIAG];
		float tmp_mul_vec3[DIAG];
		Float_Point_Algebra_obj.RND_VEC<float, DIAG>( rndnoise );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<float, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<float, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<float, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<float, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		float norm1, norm2;
		float tmp_sub_vec[DIAG];
		Float_Point_Algebra_obj.VEC_SUB<float, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		Float_Point_Algebra_obj.VEC_NORM<float, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<float, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<float> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((float)x_k_vec[i]);
		x_k_record.push_back(xkvec);
#endif// endif RECORD_RESULT
#ifdef DEBUG_ITER
		std::cout << error << std::endl << std::endl;
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
	std::cout << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-2)
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-2) << std::endl;
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
#endif// endif SHOW_FIGURE
#endif// endif PLOT_FIGURE
#endif// endif RECORD_RESULT
#endif// endif FLOAT_PRECISION
}




void PROXIMAL_GRADIENT_DECENT_GFLOAT(float **Amatrix_c,
									  float *bvector_c,
									  float L_c){


#ifdef GENERAL_FLOAT_PRECISION

	Float_Point_Algebra Float_Point_Algebra_obj;
	float error = 0;
	float error_std = ERR_STD;

	std::string clockname = "timeprofile_gfloat.txt";
	std::string xkname = "xk_gfloat.dat";
	std::string errorrecordname = "error_record_gfloat.dat";
	std::string errorhistname = "error_hist_gfloat.dat";
	std::string figurename = "ProximalGradientDecent_gfloat.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	float *x_k_vec;
   	float *x_k_plus1_vec;
   	float *grad_g;
   	float *ones_vec;
   	x_k_vec = (float*) malloc(sizeof(float)*DIAG);
   	x_k_plus1_vec = (float*) malloc(sizeof(float)*DIAG);
   	grad_g = (float*) malloc(sizeof(float)*DIAG);
   	ones_vec = (float*) malloc(sizeof(float)*DIAG);
   	std::vector<float> error_hist(ITER_MAX, 0);
   	std::vector<float> error_record(ITER_MAX, 0);
   	std::vector<std::vector<float>> x_k_record;

   	Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( ones_vec );
   	Float_Point_Algebra_obj.ZEROS_VEC<float, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<float, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		float *tmp_mul_vec1;
		tmp_mul_vec1 = (float*) malloc(sizeof(float)*DIAG);
//	   	std::cout << __FILE__ << __LINE__ << std::endl;
		Float_Point_Algebra_obj.MAT_VEC_MUL<float, DIAG, DIAG>(
				Amatrix_c, x_k_vec, tmp_mul_vec1 );
//	   	std::cout << __FILE__ << __LINE__ << std::endl;
		Float_Point_Algebra_obj.VEC_ADD<float, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		float *new_point;
		float *tmp_mul_vec2;
		new_point = (float*) malloc(sizeof(float)*DIAG);
		tmp_mul_vec2 = (float*) malloc(sizeof(float)*DIAG);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<float, DIAG>( grad_g, 1/L_c,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<float, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		float *rndnoise;
		float *a_ones_vec;
		float *minus_a_ones_vec;
		rndnoise = (float*) malloc(sizeof(float)*DIAG);
		a_ones_vec = (float*) malloc(sizeof(float)*DIAG);
		minus_a_ones_vec = (float*) malloc(sizeof(float)*DIAG);
		Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( a_ones_vec );
		Float_Point_Algebra_obj.ONES_VEC<float, DIAG>( minus_a_ones_vec );
		float *tmp_min;
		float *tmp_max;
		float *tmp_mul_vec3;
		tmp_min = (float*) malloc(sizeof(float)*DIAG);
		tmp_max = (float*) malloc(sizeof(float)*DIAG);
		tmp_mul_vec3 = (float*) malloc(sizeof(float)*DIAG);
		Float_Point_Algebra_obj.RND_VEC<float, DIAG>( rndnoise );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<float, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<float, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<float, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<float, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		float norm1, norm2;
		float *tmp_sub_vec;
		tmp_sub_vec = (float*) malloc(sizeof(float)*DIAG);
		Float_Point_Algebra_obj.VEC_SUB<float, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		Float_Point_Algebra_obj.VEC_NORM<float, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<float, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<float> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((float)x_k_vec[i]);
		x_k_record.push_back(xkvec);
#endif// endif RECORD_RESULT
#ifdef DEBUG_ITER
		std::cout << error << std::endl << std::endl;
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
	std::cout << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-2)
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist.at(k-2) << std::endl;
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
	std::vector<float> x, y;
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
#endif
#endif
#endif
#endif
}



