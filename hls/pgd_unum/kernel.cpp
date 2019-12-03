#include "kernel.hpp"

#if SOFT_POSIT_PRECISION==0

posit32_t Amatrix_c[DIAG][DIAG];
posit32_t bvector_c[DIAG];
posit32_t rndnoise[DIAG];
posit32_t L_c_inv;
posit32_t error_std;
posit32_t BOX_CONST_pos = convertDoubleToP32(BOX_CONST);
posit32_t BOX_CONST_neg = convertDoubleToP32(-BOX_CONST);

void kernel(posit32_t x_k_vec_in[DIAG], posit32_t x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
	posit32_t x_k_plus1_vec[DIAG];
	VEC_EQ<posit32_t, DIAG>(x_k_vec_in, x_k_plus1_vec);

	// Compute gradient
	posit32_t grad_g[DIAG];
	posit32_t tmp_mul_vec1[DIAG];
	MAT_VEC_MUL<posit32_t, DIAG, DIAG,
				convertDoubleToP32, p32_mul, p32_add>(
			Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
	VEC_ADD<posit32_t, DIAG, p32_add>( tmp_mul_vec1,
			bvector_c, grad_g );

	// new decent point
	posit32_t new_point[DIAG];
	posit32_t tmp_mul_vec2[DIAG];
	VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>( grad_g, L_c_inv,
										 tmp_mul_vec2 );
	VEC_SUB<posit32_t, DIAG, p32_sub>( x_k_plus1_vec, tmp_mul_vec2,
									new_point );
	#ifdef DEBUG_ITER
	std::cout << new_point << std::endl << std::endl;
	#endif

	// Proximal projection
	posit32_t tmp_min[DIAG];
	posit32_t tmp_max[DIAG];
	posit32_t tmp_mul_vec3[DIAG];
	VEC_SCALAR_MIN<posit32_t, DIAG, p32_lt>(
					new_point, BOX_CONST_pos, tmp_min);
	VEC_SCALAR_MAX<posit32_t, DIAG, p32_lt>(
					tmp_min, BOX_CONST_neg, tmp_max);
	VEC_SCALAR_MUL<posit32_t, DIAG, p32_mul>(
					rndnoise, error_std, tmp_mul_vec3);
	VEC_ADD<posit32_t, DIAG, p32_add>(
					tmp_max, tmp_mul_vec3, x_k_vec_out);
}

#elif SOFT_POSIT_PRECISION==1

posit16_t Amatrix_c[DIAG][DIAG];
posit16_t bvector_c[DIAG];
posit16_t rndnoise[DIAG];
posit16_t L_c_inv;
posit16_t error_std;
posit16_t BOX_CONST_pos = convertDoubleToP16(BOX_CONST);
posit16_t BOX_CONST_neg = convertDoubleToP16(-BOX_CONST);

void kernel(posit16_t x_k_vec_in[DIAG], posit16_t x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
	posit16_t x_k_plus1_vec[DIAG];
	VEC_EQ<posit16_t, DIAG>(x_k_vec_in, x_k_plus1_vec);

	// Compute gradient
	posit16_t grad_g[DIAG];
	posit16_t tmp_mul_vec1[DIAG];
	MAT_VEC_MUL<posit16_t, DIAG, DIAG,
				convertDoubleToP16, p16_mul, p16_add>(
			Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
	VEC_ADD<posit16_t, DIAG, p16_add>( tmp_mul_vec1,
			bvector_c, grad_g );

	// new decent point
	posit16_t new_point[DIAG];
	posit16_t tmp_mul_vec2[DIAG];
	VEC_SCALAR_MUL<posit16_t, DIAG, p16_mul>( grad_g, L_c_inv,
										 tmp_mul_vec2 );
	VEC_SUB<posit16_t, DIAG, p16_sub>( x_k_plus1_vec, tmp_mul_vec2,
									new_point );
#ifdef DEBUG_ITER
	std::cout << new_point << std::endl << std::endl;
#endif

	// Proximal projection
	posit16_t tmp_min[DIAG];
	posit16_t tmp_max[DIAG];
	posit16_t tmp_mul_vec3[DIAG];
	VEC_SCALAR_MIN<posit16_t, DIAG, p16_lt>(
					new_point, BOX_CONST_pos, tmp_min);
	VEC_SCALAR_MAX<posit16_t, DIAG, p16_lt>(
					tmp_min, BOX_CONST_neg, tmp_max);
	VEC_SCALAR_MUL<posit16_t, DIAG, p16_mul>(
					rndnoise, error_std, tmp_mul_vec3);
	VEC_ADD<posit16_t, DIAG, p16_add>(
					tmp_max, tmp_mul_vec3, x_k_vec_out);
}

#elif SOFT_POSIT_PRECISION==2

posit8_t Amatrix_c[DIAG][DIAG];
posit8_t bvector_c[DIAG];
posit8_t rndnoise[DIAG];
posit8_t L_c_inv;
posit8_t error_std;
posit8_t BOX_CONST_pos = convertDoubleToP8(BOX_CONST);
posit8_t BOX_CONST_neg = convertDoubleToP8(-BOX_CONST);

void kernel(posit8_t x_k_vec_in[DIAG], posit8_t x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
	posit8_t x_k_plus1_vec[DIAG];
	VEC_EQ<posit8_t, DIAG>(x_k_vec_in, x_k_plus1_vec);

	// Compute gradient
	posit8_t grad_g[DIAG];
	posit8_t tmp_mul_vec1[DIAG];
	MAT_VEC_MUL<posit8_t, DIAG, DIAG,
				convertDoubleToP8, p8_mul, p8_add>(
			Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
	VEC_ADD<posit8_t, DIAG, p8_add>( tmp_mul_vec1,
			bvector_c, grad_g );

	// new decent point
	posit8_t new_point[DIAG];
	posit8_t tmp_mul_vec2[DIAG];
	VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>( grad_g, L_c_inv,
										 tmp_mul_vec2 );
	VEC_SUB<posit8_t, DIAG, p8_sub>( x_k_plus1_vec, tmp_mul_vec2,
									new_point );
#ifdef DEBUG_ITER
	std::cout << new_point << std::endl << std::endl;
#endif

	// Proximal projection
	posit8_t tmp_min[DIAG];
	posit8_t tmp_max[DIAG];
	posit8_t tmp_mul_vec3[DIAG];
	VEC_SCALAR_MIN<posit8_t, DIAG, p8_lt>(
					new_point, BOX_CONST_pos, tmp_min);
	VEC_SCALAR_MAX<posit8_t, DIAG, p8_lt>(
					tmp_min, BOX_CONST_neg, tmp_max);
	VEC_SCALAR_MUL<posit8_t, DIAG, p8_mul>(
					rndnoise, error_std, tmp_mul_vec3);
	VEC_ADD<posit8_t, DIAG, p8_add>(
					tmp_max, tmp_mul_vec3, x_k_vec_out);
}

#elif SOFT_POSIT_PRECISION==3

posit_2_t Amatrix_c[DIAG][DIAG];
posit_2_t bvector_c[DIAG];
posit_2_t rndnoise[DIAG];
posit_2_t L_c_inv;
posit_2_t error_std;
posit_2_t BOX_CONST_pos = pX1_to_pX2(convertDoubleToPX1(BOX_CONST,TOTALBITS),TOTALBITS);
posit_2_t BOX_CONST_neg = pX1_to_pX2(convertDoubleToPX1(-BOX_CONST,TOTALBITS),TOTALBITS);


void kernel(posit_2_t x_k_vec_in[DIAG], posit_2_t x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in
	// Save previous point
	posit_2_t x_k_plus1_vec[DIAG];
	VEC_EQ<posit_2_t, DIAG>(x_k_vec_in, x_k_plus1_vec);

	// Compute gradient
	posit_2_t grad_g[DIAG];
	posit_2_t tmp_mul_vec1[DIAG];
	MAT_VEC_MUL<posit_2_t, posit_1_t, DIAG, DIAG,
				convertDoubleToPX1, pX1_to_pX2, pX2_mul, pX2_add, TOTALBITS>(
			Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
	VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>( tmp_mul_vec1,
			bvector_c, grad_g );

	// new decent point
	posit_2_t new_point[DIAG];
	posit_2_t tmp_mul_vec2[DIAG];
	VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>( grad_g, L_c_inv,
								 tmp_mul_vec2 );
	VEC_SUB<posit_2_t, DIAG, pX2_sub, TOTALBITS>( x_k_plus1_vec, tmp_mul_vec2,
							new_point );

	// Proximal projection
	posit_2_t tmp_min[DIAG];
	posit_2_t tmp_max[DIAG];
	posit_2_t tmp_mul_vec3[DIAG];
	VEC_SCALAR_MIN<posit_2_t, DIAG, pX2_lt>(
			new_point, BOX_CONST_pos, tmp_min);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_min, DIAG);
	VEC_SCALAR_MAX<posit_2_t, DIAG, pX2_lt>(
			tmp_min, BOX_CONST_neg, tmp_max);
//		std::cout << __FILE__ << ", " << __LINE__ << std::endl;
//	   	printvector(tmp_max, DIAG);
	VEC_SCALAR_MUL<posit_2_t, DIAG, pX2_mul, TOTALBITS>(
			rndnoise, error_std, tmp_mul_vec3);
	VEC_ADD<posit_2_t, DIAG, pX2_add, TOTALBITS>(
			tmp_max, tmp_mul_vec3, x_k_vec_out);

}

#endif
