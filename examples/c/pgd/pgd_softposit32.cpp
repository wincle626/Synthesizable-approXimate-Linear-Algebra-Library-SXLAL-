/*
 * pgd_softposit32.cpp
 *
 *  Created on: Aug 22, 2019
 *      Author: yunwu
 */

#include "pgd_softposit32.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_SPOSIT32(posit32_t Amatrix_c[DIAG][DIAG],
									  posit32_t bvector_c[DIAG],
									  posit32_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit32_t error = convertDoubleToP32(0);
	posit32_t error_std = convertDoubleToP32(ERR_STD);
	posit32_t EPS_STOP_sp = convertDoubleToP32(EPS_STOP);

	std::string clockname = "timeprofile_softposit32.txt";
	std::string xkname = "xk_softposit32.dat";
	std::string errorrecordname = "error_record_softposit32.dat";
	std::string errorhistname = "error_hist_softposit32.dat";
	std::string figurename = "ProximalGradientDecent_softposit32.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	posit32_t x_k_vec[DIAG];
   	posit32_t x_k_plus1_vec[DIAG];
   	posit32_t grad_g[DIAG];
   	posit32_t ones_vec[DIAG];
   	posit32_t spzero = convertDoubleToP32(0);
   	posit32_t spone = convertDoubleToP32(1);
   	posit32_t L_c_inv = p32_div(spone, L_c);
   	posit32_t BOX_CONST_pos = convertDoubleToP32(BOX_CONST);
   	posit32_t BOX_CONST_neg = convertDoubleToP32(-BOX_CONST);
   	std::vector<posit32_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit32_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit32_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit32_t, DIAG, convertDoubleToP32>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit32_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit32_t tmp_mul_vec1[DIAG];
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit32_t, DIAG, DIAG,
					convertDoubleToP32, p32_mul, p32_add>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit32_t, DIAG, p32_add>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit32_t new_point[DIAG];
		posit32_t tmp_mul_vec2[DIAG];
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit32_t, DIAG, p32_sub>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit32_t rndnoise[DIAG];
		posit32_t a_ones_vec[DIAG];
		posit32_t minus_a_ones_vec[DIAG];
		SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>(
				minus_a_ones_vec );
		posit32_t tmp_min[DIAG];
		posit32_t tmp_max[DIAG];
		posit32_t tmp_mul_vec3[DIAG];
		SoftPosit_Algebra_obj.RND_VEC<posit32_t, DIAG, convertDoubleToP32>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit32_t, DIAG, p32_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit32_t, DIAG, p32_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit32_t, DIAG, p32_add>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit32_t norm1, norm2;
		posit32_t tmp_sub_vec[DIAG];
		SoftPosit_Algebra_obj.VEC_SUB<posit32_t, DIAG, p32_sub>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit32_t, DIAG, convertDoubleToP32,
				p32_mul, p32_add, p32_sqrt>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit32_t, DIAG, convertDoubleToP32,
				p32_mul, p32_add, p32_sqrt>(x_k_vec, norm2);
		if( p32_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = p32_div(norm1, norm2);
		error_record.at(k) = error;
		std::vector<posit32_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit32_t)x_k_vec[i]);
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
			  << "Error: " << convertP32ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP32ToDouble(error_hist.at(k-2))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertP32ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP32ToDouble(error_hist.at(k-2)) << std::endl;
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
		if( !p32_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertP32ToDouble(error_hist.at(i));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertP32ToDouble(error_record.at(i)) << " \n";
			std::vector<posit32_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertP32ToDouble(xkvec.at(j)) << ",";
			resultfile2 << convertP32ToDouble(xkvec.at(DIAG-1)) << "\n";
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



void PROXIMAL_GRADIENT_DECENT_SPOSIT32(posit32_t **Amatrix_c,
									  posit32_t *bvector_c,
									  posit32_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit32_t error = convertDoubleToP32(0);
	posit32_t error_std = convertDoubleToP32(ERR_STD);
	posit32_t EPS_STOP_sp = convertDoubleToP32(EPS_STOP);

	std::string clockname = "timeprofile_softposit32.txt";
	std::string xkname = "xk_softposit32.dat";
	std::string errorrecordname = "error_record_softposit32.dat";
	std::string errorhistname = "error_hist_softposit32.dat";
	std::string figurename = "ProximalGradientDecent_softposit32.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	posit32_t *x_k_vec;
   	posit32_t *x_k_plus1_vec;
   	posit32_t *grad_g;
   	posit32_t *ones_vec;
   	x_k_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
   	x_k_plus1_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
   	grad_g = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
   	ones_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
   	posit32_t spzero = convertDoubleToP32(0);
   	posit32_t spone = convertDoubleToP32(1);
   	posit32_t L_c_inv = p32_div(spone, L_c);
   	posit32_t BOX_CONST_pos = convertDoubleToP32(BOX_CONST);
   	posit32_t BOX_CONST_neg = convertDoubleToP32(-BOX_CONST);
   	std::vector<posit32_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit32_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit32_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit32_t, DIAG, convertDoubleToP32>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit32_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit32_t *tmp_mul_vec1;
		tmp_mul_vec1 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit32_t, DIAG, DIAG,
					convertDoubleToP32, p32_mul, p32_add>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit32_t, DIAG, p32_add>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit32_t *new_point;
		posit32_t *tmp_mul_vec2;
		new_point = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		tmp_mul_vec2 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit32_t, DIAG, p32_sub>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit32_t *rndnoise;
		posit32_t *a_ones_vec;
		posit32_t *minus_a_ones_vec;
		rndnoise = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		a_ones_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		minus_a_ones_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit32_t, DIAG, convertDoubleToP32>(
				minus_a_ones_vec );
		posit32_t *tmp_min;
		posit32_t *tmp_max;
		posit32_t *tmp_mul_vec3;
		tmp_min = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		tmp_max = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		tmp_mul_vec3 = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		SoftPosit_Algebra_obj.RND_VEC<posit32_t, DIAG, convertDoubleToP32>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit32_t, DIAG, p32_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit32_t, DIAG, p32_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit32_t, DIAG, p32_add>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit32_t norm1, norm2;
		posit32_t *tmp_sub_vec;
		tmp_sub_vec = (posit32_t*) malloc(sizeof(posit32_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit32_t, DIAG, p32_sub>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit32_t, DIAG, convertDoubleToP32,
				p32_mul, p32_add, p32_sqrt>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit32_t, DIAG, convertDoubleToP32,
				p32_mul, p32_add, p32_sqrt>(x_k_vec, norm2);
		if( p32_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = p32_div(norm1, norm2);
		error_record.at(k) = error;
		std::vector<posit32_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit32_t)x_k_vec[i]);
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
			  << "Error: " << convertP32ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP32ToDouble(error_hist.at(k-2))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertP32ToDouble(error) << std::endl
			  << "Relative Error: "
			  << convertP32ToDouble(error_hist.at(k-2)) << std::endl;
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
		if( !p32_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertP32ToDouble(error_hist.at(i));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertP32ToDouble(error_record.at(i)) << " \n";
			std::vector<posit32_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertP32ToDouble(xkvec.at(j)) << ",";
			resultfile2 << convertP32ToDouble(xkvec.at(DIAG-1)) << "\n";
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

