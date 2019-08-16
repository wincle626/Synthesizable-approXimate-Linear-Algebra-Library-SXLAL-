/*
 * pgd_xfloat.cpp
 *
 *  Created on: Aug 15, 2019
 *      Author: yunwu
 */

#include "pgd_xfpt.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_XFLOAT1(fptx1 Amatrix_c[DIAG][DIAG],
									  fptx1 bvector_c[DIAG],
									  fptx1 L_c){

#ifdef COMSTOM_FLOAT_PRECISION
	Float_Point_Algebra Float_Point_Algebra_obj;
	fptx1 error = 0.0;
	fptx1 error_std = ERR_STD;

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

   	fptx1 x_k_vec[DIAG];
   	fptx1 x_k_plus1_vec[DIAG];
   	fptx1 grad_g[DIAG];
   	fptx1 ones_vec[DIAG];
   	std::vector<fptx1> error_hist(ITER_MAX, 0.0);
   	std::vector<fptx1> error_record(ITER_MAX, 0.0);
   	std::vector<std::vector<fptx1>> x_k_record;

   	Float_Point_Algebra_obj.ONES_VEC<fptx1, DIAG>( ones_vec );
   	Float_Point_Algebra_obj.ZEROS_VEC<fptx1, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<fptx1, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		fptx1 tmp_mul_vec1[DIAG];
		Float_Point_Algebra_obj.MAT_VEC_MUL<fptx1, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
		Float_Point_Algebra_obj.VEC_ADD<fptx1, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		fptx1 new_point[DIAG];
		fptx1 tmp_mul_vec2[DIAG];
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<fptx1, DIAG>( grad_g, 1/L_c,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<fptx1, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		fptx1 rndnoise[DIAG];
		fptx1 a_ones_vec[DIAG];
		fptx1 minus_a_ones_vec[DIAG];
		Float_Point_Algebra_obj.ONES_VEC<fptx1, DIAG>( a_ones_vec );
		Float_Point_Algebra_obj.ONES_VEC<fptx1, DIAG>( minus_a_ones_vec );
		fptx1 tmp_min[DIAG];
		fptx1 tmp_max[DIAG];
		fptx1 tmp_mul_vec3[DIAG];
		Float_Point_Algebra_obj.RND_VEC<fptx1, DIAG>( rndnoise );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<fptx1, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<fptx1, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<fptx1, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<fptx1, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		fptx1 norm1, norm2;
		fptx1 tmp_sub_vec[DIAG];
		Float_Point_Algebra_obj.VEC_SUB<fptx1, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		Float_Point_Algebra_obj.VEC_NORM<fptx1, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<fptx1, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<fptx1> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((fptx1)x_k_vec[i]);
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
			std::vector<fptx1> xkvec = x_k_record.at(i);
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


void PROXIMAL_GRADIENT_DECENT_XFLOAT2(fptx2 Amatrix_c[DIAG][DIAG],
									  fptx2 bvector_c[DIAG],
									  fptx2 L_c){

#ifdef COMSTOM_FLOAT_PRECISION
	Float_Point_Algebra Float_Point_Algebra_obj;
	fptx2 error = 0.0;
	fptx2 error_std = ERR_STD;

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

   	fptx2 x_k_vec[DIAG];
   	fptx2 x_k_plus1_vec[DIAG];
   	fptx2 grad_g[DIAG];
   	fptx2 ones_vec[DIAG];
   	std::vector<fptx2> error_hist(ITER_MAX, 0.0);
   	std::vector<fptx2> error_record(ITER_MAX, 0.0);
   	std::vector<std::vector<fptx2>> x_k_record;

   	Float_Point_Algebra_obj.ONES_VEC<fptx2, DIAG>( ones_vec );
   	Float_Point_Algebra_obj.ZEROS_VEC<fptx2, DIAG>( x_k_vec );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		Float_Point_Algebra_obj.VEC_EQ<fptx2, DIAG>(x_k_vec, x_k_plus1_vec);

		// Compute gradient
		fptx2 tmp_mul_vec1[DIAG];
		Float_Point_Algebra_obj.MAT_VEC_MUL<fptx2, DIAG, DIAG>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
		Float_Point_Algebra_obj.VEC_ADD<fptx2, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		fptx2 new_point[DIAG];
		fptx2 tmp_mul_vec2[DIAG];
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<fptx2, DIAG>( grad_g, 1/L_c,
									 tmp_mul_vec2 );
		Float_Point_Algebra_obj.VEC_SUB<fptx2, DIAG>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		fptx2 rndnoise[DIAG];
		fptx2 a_ones_vec[DIAG];
		fptx2 minus_a_ones_vec[DIAG];
		Float_Point_Algebra_obj.ONES_VEC<fptx2, DIAG>( a_ones_vec );
		Float_Point_Algebra_obj.ONES_VEC<fptx2, DIAG>( minus_a_ones_vec );
		fptx2 tmp_min[DIAG];
		fptx2 tmp_max[DIAG];
		fptx2 tmp_mul_vec3[DIAG];
		Float_Point_Algebra_obj.RND_VEC<fptx2, DIAG>( rndnoise );
		Float_Point_Algebra_obj.VEC_SCALAR_MIN<fptx2, DIAG>(new_point, BOX_CONST, tmp_min);
		Float_Point_Algebra_obj.VEC_SCALAR_MAX<fptx2, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		Float_Point_Algebra_obj.VEC_SCALAR_MUL<fptx2, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		Float_Point_Algebra_obj.VEC_ADD<fptx2, DIAG>( tmp_max, tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		fptx2 norm1, norm2;
		fptx2 tmp_sub_vec[DIAG];
		Float_Point_Algebra_obj.VEC_SUB<fptx2, DIAG>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		Float_Point_Algebra_obj.VEC_NORM<fptx2, DIAG>(tmp_sub_vec, norm1);
		Float_Point_Algebra_obj.VEC_NORM<fptx2, DIAG>(x_k_vec, norm2);
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<fptx2> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((fptx2)x_k_vec[i]);
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
			std::vector<fptx2> xkvec = x_k_record.at(i);
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

