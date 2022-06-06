#include "kernel.hpp"

float Amatrix_c[DIAG][DIAG];
float bvector_c[DIAG];
float rndnoise[DIAG];
float L_c_inv;
float error_std;

void kernel(float x_k_vec_in[DIAG], float x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
		float x_k_plus1_vec[DIAG];
		VEC_EQ<float, DIAG>(x_k_vec_in, x_k_plus1_vec);

		// Compute gradient
		float grad_g[DIAG];
		float tmp_mul_vec1[DIAG];
		MAT_VEC_MUL<float, DIAG, DIAG>(
				Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
		VEC_ADD<float, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );

		// new decent point
		float new_point[DIAG];
		float tmp_mul_vec2[DIAG];
		VEC_SCALAR_MUL<float, DIAG>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
		VEC_SUB<float, DIAG>( x_k_plus1_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		float tmp_min[DIAG];
		float tmp_max[DIAG];
		float tmp_mul_vec3[DIAG];
		VEC_SCALAR_MIN<float, DIAG>(new_point, BOX_CONST, tmp_min);
		VEC_SCALAR_MAX<float, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		VEC_SCALAR_MUL<float, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		VEC_ADD<float, DIAG>( tmp_max, tmp_mul_vec3,
							  x_k_vec_out);
}
