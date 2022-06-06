/*
 * pgd_softposit8.cpp
 *
 *  Created on: 29 Aug 2019
 *      Author: yunwu
 */


#include "pgd_softposit8.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_SPOSIT8(posit8_t Amatrix_c[DIAG][DIAG],
									  posit8_t bvector_c[DIAG],
									  posit8_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit8_t error = convertDoubleToP8(0);
	posit8_t error_std = convertDoubleToP8(ERR_STD);
	posit8_t EPS_STOP_sp = convertDoubleToP8(EPS_STOP);

	std::string clockname = "timeprofile_softposit8.txt";
	std::string xkname = "xk_softposit8.dat";
	std::string errorrecordname = "error_record_softposit8.dat";
	std::string errorhistname = "error_hist_softposit8.dat";
	std::string figurename = "ProximalGradientDecent_softposit8.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	posit8_t x_k_vec[DIAG];
   	posit8_t x_k_plus1_vec[DIAG];
   	posit8_t grad_g[DIAG];
   	posit8_t ones_vec[DIAG];
   	posit8_t spzero = convertDoubleToP8(0);
   	posit8_t spone = convertDoubleToP8(1);
   	posit8_t L_c_inv = p8_div(spone, L_c);
   	posit8_t BOX_CONST_pos = convertDoubleToP8(BOX_CONST);
   	posit8_t BOX_CONST_neg = convertDoubleToP8(-BOX_CONST);
   	std::vector<posit8_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit8_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit8_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit8_t, DIAG, convertDoubleToP8>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit8_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit8_t tmp_mul_vec1[DIAG];
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit8_t, DIAG, DIAG,
					convertDoubleToP8, p8_mul, p8_add>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit8_t, DIAG, p8_add>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit8_t new_point[DIAG];
		posit8_t tmp_mul_vec2[DIAG];
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit8_t, DIAG, p8_sub>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit8_t rndnoise[DIAG];
		posit8_t a_ones_vec[DIAG];
		posit8_t minus_a_ones_vec[DIAG];
		SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>(
				minus_a_ones_vec );
		posit8_t tmp_min[DIAG];
		posit8_t tmp_max[DIAG];
		posit8_t tmp_mul_vec3[DIAG];
		SoftPosit_Algebra_obj.RND_VEC<posit8_t, DIAG, convertDoubleToP8>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit8_t, DIAG, p8_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit8_t, DIAG, p8_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit8_t, DIAG, p8_add>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit8_t norm1, norm2;
		posit8_t tmp_sub_vec[DIAG];
		SoftPosit_Algebra_obj.VEC_SUB<posit8_t, DIAG, p8_sub>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit8_t, DIAG, convertDoubleToP8,
				p8_mul, p8_add, p8_sqrt>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit8_t, DIAG, convertDoubleToP8,
				p8_mul, p8_add, p8_sqrt>(x_k_vec, norm2);
		if( p8_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = p8_div(norm1, norm2);
		error_record.at(k) = error;
		std::vector<posit8_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit8_t)x_k_vec[i]);
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
			  << "Error: " << convertP8ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP8ToDouble(error_hist.at(k-2))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertP8ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP8ToDouble(error_hist.at(k-2)) << std::endl;
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
	for(int i=0; i<(int)x_k_record.size(); i++){
		if( !p8_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertP8ToDouble(error_hist.at(i));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertP8ToDouble(error_record.at(i)) << " \n";
			std::vector<posit8_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertP8ToDouble(xkvec.at(j)) << ",";
			resultfile2 << convertP8ToDouble(xkvec.at(DIAG-1)) << "\n";
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



void PROXIMAL_GRADIENT_DECENT_SPOSIT8(posit8_t **Amatrix_c,
									  posit8_t *bvector_c,
									  posit8_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit8_t error = convertDoubleToP8(0);
	posit8_t error_std = convertDoubleToP8(ERR_STD);
	posit8_t EPS_STOP_sp = convertDoubleToP8(EPS_STOP);

	std::string clockname = "timeprofile_softposit8.txt";
	std::string xkname = "xk_softposit8.dat";
	std::string errorrecordname = "error_record_softposit8.dat";
	std::string errorhistname = "error_hist_softposit8.dat";
	std::string figurename = "ProximalGradientDecent_softposit8.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	posit8_t *x_k_vec;
   	posit8_t *x_k_plus1_vec;
   	posit8_t *grad_g;
   	posit8_t *ones_vec;
   	x_k_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
   	x_k_plus1_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
   	grad_g = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
   	ones_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
   	posit8_t spzero = convertDoubleToP8(0);
   	posit8_t spone = convertDoubleToP8(1);
   	posit8_t L_c_inv = p8_div(spone, L_c);
   	posit8_t BOX_CONST_pos = convertDoubleToP8(BOX_CONST);
   	posit8_t BOX_CONST_neg = convertDoubleToP8(-BOX_CONST);
   	std::vector<posit8_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit8_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit8_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit8_t, DIAG, convertDoubleToP8>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit8_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit8_t *tmp_mul_vec1;
		tmp_mul_vec1 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit8_t, DIAG, DIAG,
					convertDoubleToP8, p8_mul, p8_add>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit8_t, DIAG, p8_add>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit8_t *new_point;
		posit8_t *tmp_mul_vec2;
		new_point = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		tmp_mul_vec2 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit8_t, DIAG, p8_sub>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit8_t *rndnoise;
		posit8_t *a_ones_vec;
		posit8_t *minus_a_ones_vec;
		rndnoise = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		a_ones_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		minus_a_ones_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit8_t, DIAG, convertDoubleToP8>(
				minus_a_ones_vec );
		posit8_t *tmp_min;
		posit8_t *tmp_max;
		posit8_t *tmp_mul_vec3;
		tmp_min = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		tmp_max = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		tmp_mul_vec3 = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		SoftPosit_Algebra_obj.RND_VEC<posit8_t, DIAG, convertDoubleToP8>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit8_t, DIAG, p8_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit8_t, DIAG, p8_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit8_t, DIAG, p8_add>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit8_t norm1, norm2;
		posit8_t *tmp_sub_vec;
		tmp_sub_vec = (posit8_t*) malloc(sizeof(posit8_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit8_t, DIAG, p8_sub>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit8_t, DIAG, convertDoubleToP8,
				p8_mul, p8_add, p8_sqrt>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit8_t, DIAG, convertDoubleToP8,
				p8_mul, p8_add, p8_sqrt>(x_k_vec, norm2);
		if( p8_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = p8_div(norm1, norm2);
		error_record.at(k) = error;
		std::vector<posit8_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit8_t)x_k_vec[i]);
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
			  << "Error: " << convertP8ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP8ToDouble(error_hist.at(k-2))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertP8ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP8ToDouble(error_hist.at(k-2)) << std::endl;
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
	for(int i=0; i<(int)x_k_record.size(); i++){
		if( !p8_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertP8ToDouble(error_hist.at(i));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertP8ToDouble(error_record.at(i)) << " \n";
			std::vector<posit8_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertP8ToDouble(xkvec.at(j)) << ",";
			resultfile2 << convertP8ToDouble(xkvec.at(DIAG-1)) << "\n";
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



