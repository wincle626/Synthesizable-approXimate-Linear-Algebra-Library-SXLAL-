#include "kernel.hpp"

double Amatrix_c[DIAG][DIAG];
double bvector_c[DIAG];
double rndnoise[DIAG];
double L_c_inv;
double error_std;

void kernel(double x_k_vec_in[DIAG], double x_k_vec_out[DIAG]){
#pragma HLS INTERFACE axis register both port=x_k_vec_out
#pragma HLS INTERFACE axis register both port=x_k_vec_in

		// Save previous point
		double x_k_plus1_vec[DIAG];
		VEC_EQ<double, DIAG>(x_k_vec_in, x_k_plus1_vec);

		// Compute gradient
		double grad_g[DIAG];
		double tmp_mul_vec1[DIAG];
		MAT_VEC_MUL<double, DIAG, DIAG>(
				Amatrix_c, x_k_plus1_vec,	tmp_mul_vec1 );
		VEC_ADD<double, DIAG>( tmp_mul_vec1, bvector_c,
								grad_g );

		// new decent point
		double new_point[DIAG];
		double tmp_mul_vec2[DIAG];
		VEC_SCALAR_MUL<double, DIAG>( grad_g, L_c_inv,
									 tmp_mul_vec2 );
		VEC_SUB<double, DIAG>( x_k_plus1_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		double tmp_min[DIAG];
		double tmp_max[DIAG];
		double tmp_mul_vec3[DIAG];
		VEC_SCALAR_MIN<double, DIAG>(new_point, BOX_CONST, tmp_min);
		VEC_SCALAR_MAX<double, DIAG>(tmp_min, -BOX_CONST, tmp_max);
		VEC_SCALAR_MUL<double, DIAG>( rndnoise, error_std,
									 tmp_mul_vec3);
		VEC_ADD<double, DIAG>( tmp_max, tmp_mul_vec3,
							  x_k_vec_out);
}
