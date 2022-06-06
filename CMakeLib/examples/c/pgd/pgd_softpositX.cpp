/*
 * pgd_softpositX.cpp
 *
 *  Created on: 4 Nov 2019
 *      Author: yw106
 */



#include "pgd_softpositX.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_SPOSITX(posit_2_t Amatrix_c[DIAG][DIAG],
									  posit_2_t bvector_c[DIAG],
									  posit_2_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit_2_t error = pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS),TOTALBITS);
	posit_2_t error_std = pX1_to_pX2(convertDoubleToPX1(ERR_STD,TOTALBITS),TOTALBITS);
	posit_2_t EPS_STOP_sp = pX1_to_pX2(convertDoubleToPX1(EPS_STOP,TOTALBITS),TOTALBITS);

	std::string clockname = "timeprofile_softposit16.txt";
	std::string xkname = "xk_softposit16.dat";
	std::string errorrecordname = "error_record_softposit16.dat";
	std::string errorhistname = "error_hist_softposit16.dat";
	std::string figurename = "ProximalGradientDecent_softposit16.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

   	posit_2_t x_k_vec[DIAG];
   	posit_2_t x_k_plus1_vec[DIAG];
   	posit_2_t grad_g[DIAG];
   	posit_2_t ones_vec[DIAG];
   	posit_2_t spzero = pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS), TOTALBITS);
   	posit_2_t spone = pX1_to_pX2(convertDoubleToPX1(1,TOTALBITS), TOTALBITS);
   	posit_2_t L_c_inv = pX2_div(spone, L_c,TOTALBITS);
   	posit_2_t BOX_CONST_pos = pX1_to_pX2(convertDoubleToPX1(BOX_CONST,TOTALBITS), TOTALBITS);
   	posit_2_t BOX_CONST_neg = pX1_to_pX2(convertDoubleToPX1(-BOX_CONST,TOTALBITS), TOTALBITS);
   	std::vector<posit_2_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit_2_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit_2_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit_2_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit_2_t tmp_mul_vec1[DIAG];
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit_2_t, posit_1_t, DIAG, DIAG,
					convertDoubleToPX1, pX1_to_pX2, pX2_mul, pX2_add, TOTALBITS>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit_2_t new_point[DIAG];
		posit_2_t tmp_mul_vec2[DIAG];
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit_2_t, DIAG, pX2_sub, TOTALBITS>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit_2_t rndnoise[DIAG];
		posit_2_t a_ones_vec[DIAG];
		posit_2_t minus_a_ones_vec[DIAG];
		SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				minus_a_ones_vec );
		posit_2_t tmp_min[DIAG];
		posit_2_t tmp_max[DIAG];
		posit_2_t tmp_mul_vec3[DIAG];
		SoftPosit_Algebra_obj.RND_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit_2_t, DIAG, pX2_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit_2_t, DIAG, pX2_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit_2_t norm1, norm2;
		posit_2_t tmp_sub_vec[DIAG];
		SoftPosit_Algebra_obj.VEC_SUB<posit_2_t, DIAG, pX2_sub, TOTALBITS>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2,
				pX2_mul, pX2_add, pX2_sqrt, TOTALBITS>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2,
		pX2_mul, pX2_add, pX2_sqrt, TOTALBITS>(x_k_vec, norm2);
		if( pX2_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = pX2_div(norm1, norm2,TOTALBITS);
		error_record.at(k) = error;
		std::vector<posit_2_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit_2_t)x_k_vec[i]);
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
			  << "Error: " << convertPX1ToDouble(pX2_to_pX1(error, TOTALBITS)) << std::endl
			  << "Relative Error: "
			  << convertPX1ToDouble(pX2_to_pX1(error_hist.at(k-2), TOTALBITS))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertPX1ToDouble(pX2_to_pX1(error, TOTALBITS)) << std::endl
			  << "Relative Error: "
			  << convertPX1ToDouble(pX2_to_pX1(error_hist.at(k-2), TOTALBITS)) << std::endl;
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
		if( !pX2_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertPX1ToDouble(pX2_to_pX1(error_hist.at(i), TOTALBITS));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertPX1ToDouble(pX2_to_pX1(error_record.at(i), TOTALBITS)) << " \n";
			std::vector<posit_2_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertPX1ToDouble(pX2_to_pX1(xkvec.at(j), TOTALBITS)) << ",";
			resultfile2 << convertPX1ToDouble(pX2_to_pX1(xkvec.at(DIAG-1), TOTALBITS)) << "\n";
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



void PROXIMAL_GRADIENT_DECENT_SPOSITX(posit_2_t **Amatrix_c,
									  posit_2_t *bvector_c,
									  posit_2_t L_c){

#ifdef SOFT_POSIT_PRECISION

	SoftPosit_Algebra SoftPosit_Algebra_obj;
	posit_2_t error = pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS),TOTALBITS);
	posit_2_t error_std = pX1_to_pX2(convertDoubleToPX1(ERR_STD,TOTALBITS),TOTALBITS);
	posit_2_t EPS_STOP_sp = pX1_to_pX2(convertDoubleToPX1(EPS_STOP,TOTALBITS),TOTALBITS);

	std::string clockname = "timeprofile_softposit16.txt";
	std::string xkname = "xk_softposit16.dat";
	std::string errorrecordname = "error_record_softposit16.dat";
	std::string errorhistname = "error_hist_softposit16.dat";
	std::string figurename = "ProximalGradientDecent_softposit16.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif

	posit_2_t *x_k_vec;
	posit_2_t *x_k_plus1_vec;
	posit_2_t *grad_g;
	posit_2_t *ones_vec;
   	x_k_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
   	x_k_plus1_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
   	grad_g = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
   	ones_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
   	posit_2_t spzero = pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS), TOTALBITS);
   	posit_2_t spone = pX1_to_pX2(convertDoubleToPX1(1,TOTALBITS), TOTALBITS);
   	posit_2_t L_c_inv = pX2_div(spone, L_c,TOTALBITS);
   	posit_2_t BOX_CONST_pos = pX1_to_pX2(convertDoubleToPX1(BOX_CONST,TOTALBITS), TOTALBITS);
   	posit_2_t BOX_CONST_neg = pX1_to_pX2(convertDoubleToPX1(-BOX_CONST,TOTALBITS), TOTALBITS);
   	std::vector<posit_2_t> error_hist(ITER_MAX, spzero);
   	std::vector<posit_2_t> error_record(ITER_MAX, spzero);
   	std::vector<std::vector<posit_2_t>> x_k_record;

   	SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>( ones_vec );
   	SoftPosit_Algebra_obj.ZEROS_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>( x_k_vec );
//	std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//   	printvector(ones_vec, DIAG);
//   	printvector(x_k_vec, DIAG);

	int k=0;
	while( k<ITER_MAX ){
//	while( k<2 ){

		// Save previous point
		SoftPosit_Algebra_obj.VEC_EQ<posit_2_t, DIAG>(x_k_vec, x_k_plus1_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_plus1_vec, DIAG);

		// Compute gradient
		posit_2_t *tmp_mul_vec1;
		tmp_mul_vec1 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		SoftPosit_Algebra_obj.MAT_VEC_MUL<posit_2_t, posit_1_t, DIAG, DIAG,
					convertDoubleToPX1, pX1_to_pX2, pX2_mul, pX2_add, TOTALBITS>(
				Amatrix_c, x_k_vec,	tmp_mul_vec1 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec1, DIAG);
		SoftPosit_Algebra_obj.VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>( tmp_mul_vec1,
				bvector_c, grad_g );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(grad_g, DIAG);
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		posit_2_t *new_point;
		posit_2_t *tmp_mul_vec2;
		new_point = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		tmp_mul_vec2 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_mul_vec2, DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit_2_t, DIAG, pX2_sub, TOTALBITS>( x_k_vec, tmp_mul_vec2,
								new_point );
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(new_point, DIAG);
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		posit_2_t *rndnoise;
		posit_2_t *a_ones_vec;
		posit_2_t *minus_a_ones_vec;
		rndnoise = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		a_ones_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		minus_a_ones_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				a_ones_vec );
		SoftPosit_Algebra_obj.ONES_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				minus_a_ones_vec );
		posit_2_t *tmp_min;
		posit_2_t *tmp_max;
		posit_2_t *tmp_mul_vec3;
		tmp_min = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		tmp_max = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		tmp_mul_vec3 = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		SoftPosit_Algebra_obj.RND_VEC<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2, TOTALBITS>(
				rndnoise );
		SoftPosit_Algebra_obj.VEC_SCALAR_MIN<posit_2_t, DIAG, pX2_lt>(
				new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MAX<posit_2_t, DIAG, pX2_lt>(
				tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
		SoftPosit_Algebra_obj.VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>(
				rndnoise, error_std, tmp_mul_vec3);
		SoftPosit_Algebra_obj.VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>(
				tmp_max, tmp_mul_vec3, x_k_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(x_k_vec, DIAG);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		posit_2_t norm1, norm2;
		posit_2_t *tmp_sub_vec;
		tmp_sub_vec = (posit_2_t*) malloc(sizeof(posit_2_t)*DIAG);
		SoftPosit_Algebra_obj.VEC_SUB<posit_2_t, DIAG, pX2_sub, TOTALBITS>( x_k_vec,
				x_k_plus1_vec, tmp_sub_vec);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_sub_vec, DIAG);
		SoftPosit_Algebra_obj.VEC_NORM<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2,
				pX2_mul, pX2_add, pX2_sqrt, TOTALBITS>(tmp_sub_vec, norm1);
		SoftPosit_Algebra_obj.VEC_NORM<posit_2_t, posit_1_t, DIAG, convertDoubleToPX1, pX1_to_pX2,
		pX2_mul, pX2_add, pX2_sqrt, TOTALBITS>(x_k_vec, norm2);
		if( pX2_le(norm1, EPS_STOP_sp) ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist.at(k) = pX2_div(norm1, norm2,TOTALBITS);
		error_record.at(k) = error;
		std::vector<posit_2_t> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((posit_2_t)x_k_vec[i]);
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
			  << "Error: " << convertPX1ToDouble(pX2_to_pX1(error, TOTALBITS)) << std::endl
			  << "Relative Error: "
			  << convertPX1ToDouble(pX2_to_pX1(error_hist.at(k-2), TOTALBITS))
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << convertPX1ToDouble(pX2_to_pX1(error, TOTALBITS)) << std::endl
			  << "Relative Error: "
			  << convertPX1ToDouble(pX2_to_pX1(error_hist.at(k-2), TOTALBITS)) << std::endl;
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
		if( !pX2_lt(error_hist.at(i),spzero) ){
			x.push_back(i);
			double errorhist = convertPX1ToDouble(pX2_to_pX1(error_hist.at(i), TOTALBITS));
			y.push_back(errorhist);
			resultfile << "  " << i << "," << errorhist << " \n";
			resultfile1 << "  " << i << ","
					    << convertPX1ToDouble(pX2_to_pX1(error_record.at(i), TOTALBITS)) << " \n";
			std::vector<posit_2_t> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << convertPX1ToDouble(pX2_to_pX1(xkvec.at(j), TOTALBITS)) << ",";
			resultfile2 << convertPX1ToDouble(pX2_to_pX1(xkvec.at(DIAG-1), TOTALBITS)) << "\n";
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
