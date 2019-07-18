/*
 * proximal_gradient_decent_double.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: yunwu
 */

#include "pgd_gdouble.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_GDOUBLE(double Amatrix_c[DIAG][DIAG],
									  double bvector_c[DIAG],
									  double L_c){

#ifdef GENERAL_DOUBLE_PRECISION

	Float_Point_Algebra Float_Point_Algebra_obj;
	double error = 0;
	double error_std = ERR_STD;

	std::string clockname = "timeprofile.txt";
	std::string xkname = "xk.dat";
	std::string errorrecordname = "error_record.dat";
	std::string errorhistname = "error_hist.dat";
	std::string figurename = "ProximalGradientDecent.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	double x_k_vec[DIAG];
   	double x_k_plus1_vec[DIAG];
   	double grad_g[DIAG];
   	double ones_vec[DIAG];
   	std::vector<double> error_hist(ITER_MAX, 0);
   	std::vector<double> error_record(ITER_MAX, 0);
   	std::vector<std::vector<double>> x_k_record;

   	Float_Point_Algebra_obj.ONES_VEC<double, DIAG>( ones_vec );
   	Float_Point_Algebra_obj.ZEROS_VEC<double, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<double, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		double tmp_mul_vec1[DIAG];
		Float_Point_Algebra_obj.GENERAL_MAT_VEC_MUL_BASIC<double, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
		Float_Point_Algebra_obj.VEC_ADD<double, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		double new_point[DIAG];
		double tmp_mul_vec2[DIAG];
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<double, DIAG>( grad_g, 1/L_c,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<double, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		double rndnoise[DIAG];
		double a_ones_vec[DIAG];
		double minus_a_ones_vec[DIAG];
		Float_Point_Algebra_obj.ONES_VEC<double, DIAG>( a_ones_vec );
		Float_Point_Algebra_obj.ONES_VEC<double, DIAG>( minus_a_ones_vec );
		double tmp_min[DIAG];
		double tmp_max[DIAG];
		double tmp_mul_vec3[DIAG];
		Float_Point_Algebra_obj.RND_VEC<double, DIAG>( rndnoise );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<double, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<double, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<double, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<double, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		double norm1, norm2;
		double tmp_sub_vec[DIAG];
		Float_Point_Algebra_obj.VEC_SUB<double, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		Float_Point_Algebra_obj.VEC_NORM<double, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<double, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<double> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((double)x_k_vec[i]);
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
			std::vector<double> xkvec = x_k_record.at(i);
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



