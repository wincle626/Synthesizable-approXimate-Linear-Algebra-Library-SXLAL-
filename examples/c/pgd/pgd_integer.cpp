/*
 * pgd_integer.cpp
 *
 *  Created on: 15 Oct 2019
 *      Author: yw106
 */




#include "pgd_integer.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_INTEGER(int Amatrix_c[DIAG][DIAG],
									 int bvector_c[DIAG],
									 int L_c_inv){

#ifdef GENERAL_INTEGER_PRECISION
	Float_Point_Algebra Float_Point_Algebra_obj;
	float error = 0;
	float error_std = ERR_STD;

	std::string clockname = "timeprofile_int.txt";
	std::string xkname = "xk_int.dat";
	std::string errorrecordname = "error_record_int.dat";
	std::string errorhistname = "error_hist_int.dat";
	std::string figurename = "ProximalGradientDecent_int.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	int x_k_vec[DIAG];
   	int x_k_plus1_vec[DIAG];
   	int grad_g[DIAG];
   	std::vector<float> error_hist(ITER_MAX, 0);
   	std::vector<float> error_record(ITER_MAX, 0);
   	std::vector<std::vector<int>> x_k_record;

   	Float_Point_Algebra_obj.ZEROS_VEC<int, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<int, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		int tmp_mul_vec1[DIAG];
		Float_Point_Algebra_obj.MAT_VEC_MUL<int, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
		Float_Point_Algebra_obj.VEC_ADD<int, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g[0] << std::endl << std::endl;
#endif

		//std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// new decent point
		int new_point[DIAG];
		int tmp_mul_vec2[DIAG];
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<int, DIAG>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<int, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point[0] << ","
				  << new_point[1] << ","
				  << new_point[2] << ","
				  << new_point[3] << ","
				  << std::endl << std::endl;
#endif

		//std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// Proximal projection
		int rndnoise[DIAG];
		int tmp_min[DIAG];
		int tmp_max[DIAG];
		int tmp_mul_vec3[DIAG];
		Float_Point_Algebra_obj.RND_VEC_SCALE<int, DIAG>( rndnoise, (int)std::sqrt(PGD_INT_SCALE) );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<int, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<int, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<int, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<int, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec[0] << ","
				  << x_k_vec[1] << ","
				  << x_k_vec[2] << ","
				  << x_k_vec[3] << ","
				  << std::endl << std::endl;
#endif

		//std::cout << __FILE__ << "," << __LINE__ << std::endl;

		// check early termination constraint
		int norm1, norm2;
		int tmp_sub_vec[DIAG];
		Float_Point_Algebra_obj.VEC_SUB<int, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
#ifdef DEBUG_ITER
		std::cout << tmp_sub_vec[0] << ","
				  << tmp_sub_vec[1] << ","
				  << tmp_sub_vec[2] << ","
				  << tmp_sub_vec[3] << ","
				  << std::endl << std::endl;
#endif
		Float_Point_Algebra_obj.VEC_NORM<int, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<int, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP_SCALE && norm1>0){
			std::cout << "terminated at " << k << "th iteration: " << norm1 << std::endl;
			break;
		}

		//std::cout << __FILE__ << "," << __LINE__ << std::endl;
		//std::cout << norm1 << "," << norm2 << std::endl;
		// record relative error
#ifdef RECORD_RESULT
		error = norm1/PGD_INT_SCALE/PGD_INT_SCALE;
		error_hist.at(k) = (float)norm1 / norm2;
		error_record.at(k) = error;
		std::vector<int> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((int)x_k_vec[i]);
		x_k_record.push_back(xkvec);
#endif// endif RECORD_RESULT
#ifdef DEBUG_ITER
		std::cout << error << std::endl << std::endl;
#endif

		//std::cout << __FILE__ << "," << __LINE__ << std::endl;

		k++;
		//getchar();

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
	if(k>=2){
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
		std::vector<int> x, y;
		for(int i=0; i<ITER_MAX; i++){
			if( error_hist.at(i)>0 ){
				x.push_back(i);
				y.push_back(error_hist.at(i));
				resultfile << "  " << i << "," << error_hist.at(i) << " \n";
				resultfile1 << "  " << i << "," << error_record.at(i) << " \n";
				std::vector<int> xkvec = x_k_record.at(i);
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
	}
#endif// endif int_PRECISION
}





