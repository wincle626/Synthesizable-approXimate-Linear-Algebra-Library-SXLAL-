#include "kernel.hpp"

void kernel(float Atb[DIAG], float x[DIAG]){

	float q[DIAG];
	float zold[DIAG];
	float zminusu[DIAG];
	float rhozminusu[DIAG];
	float invLq[DIAG];
	float alphax[DIAG];
	float oneminusalphazold[DIAG];
	float x_hat[DIAG];
	float x_hatu[DIAG];
	float x_hatu1[DIAG];
	float x_hatu2[DIAG];
	float x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<float,DIAG>(z, u, zminusu);
	VEC_SCALAR_MUL<float,DIAG>(zminusu,
			rho, rhozminusu);
	VEC_ADD<float,DIAG>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<float,DIAG,DIAG>(
			invL, q, invLq);
	MAT_VEC_MUL<float,DIAG,DIAG>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<float,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<float,DIAG>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<float,DIAG>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<float,DIAG>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<float,DIAG>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<float,DIAG>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<float,DIAG>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<float,DIAG>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<float,DIAG>(
			x_hatu1, 0, x_hatu1);
	VEC_SCALAR_MAX<float,DIAG>(
			x_hatu2, 0, x_hatu2);
	VEC_SUB<float,DIAG>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<float,DIAG>(x_hat, z, x_hatz);
	VEC_ADD<float,DIAG>(u, x_hatz, u);

}
