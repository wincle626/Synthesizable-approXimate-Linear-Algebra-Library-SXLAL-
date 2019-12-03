#include "kernel.hpp"

DATA_IN_T Amatrix_c[DIAG][DIAG];
DATA_IN_T bvector_c[DIAG];
DATA_IN_T rndnoise[DIAG];
DATA_IN_T L_c_inv;
DATA_IN_T error_std;

void kernel(DATA_IN_T x_k_vec_in[DIAG], DATA_IN_T x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
		DATA_IN_T x_k_plus1_vec[DIAG];
		VEC_EQ<DATA_IN_T, DIAG>(x_k_vec_in, x_k_plus1_vec);

		// Compute gradient
		DATA_IN_T grad_g[DIAG];
		DATA_IN_T tmp_mul_vec1[DIAG];
		MAT_VEC_MUL<DATA_IN_T, DIAG, DIAG>(
				Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
		VEC_ADD<DATA_IN_T, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );

		// new decent point
		DATA_IN_T new_point[DIAG];
		DATA_IN_T tmp_mul_vec2[DIAG];
		VEC_SCALAR_MUL<DATA_IN_T, DIAG>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
		VEC_SUB<DATA_IN_T, DIAG>( x_k_plus1_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		DATA_IN_T tmp_min[DIAG];
		DATA_IN_T tmp_max[DIAG];
		DATA_IN_T tmp_mul_vec3[DIAG];
		VEC_SCALAR_MIN<DATA_IN_T, DIAG>(new_point, BOX_CONST, tmp_min);
		VEC_SCALAR_MAX<DATA_IN_T, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		VEC_SCALAR_MUL<DATA_IN_T, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		VEC_ADD<DATA_IN_T, DIAG>( tmp_max, tmp_mul_vec3,
							  x_k_vec_out);
}
