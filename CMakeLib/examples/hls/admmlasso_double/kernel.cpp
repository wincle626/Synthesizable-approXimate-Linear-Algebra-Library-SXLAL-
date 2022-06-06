#include "kernel.hpp"

void kernel(double Atb[DIAG], double x[DIAG]){

	double q[DIAG];
	double zold[DIAG];
	double zminusu[DIAG];
	double rhozminusu[DIAG];
	double invLq[DIAG];
	double alphax[DIAG];
	double oneminusalphazold[DIAG];
	double x_hat[DIAG];
	double x_hatu[DIAG];
	double x_hatu1[DIAG];
	double x_hatu2[DIAG];
	double x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<double,DIAG>(z, u, zminusu);
	VEC_SCALAR_MUL<double,DIAG>(zminusu,
			rho, rhozminusu);
	VEC_ADD<double,DIAG>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<double,DIAG,DIAG>(
			invL, q, invLq);
	MAT_VEC_MUL<double,DIAG,DIAG>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<double,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<double,DIAG>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<double,DIAG>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<double,DIAG>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<double,DIAG>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<double,DIAG>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<double,DIAG>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<double,DIAG>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<double,DIAG>(
			x_hatu1, 0, x_hatu1);
	VEC_SCALAR_MAX<double,DIAG>(
			x_hatu2, 0, x_hatu2);
	VEC_SUB<double,DIAG>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<double,DIAG>(x_hat, z, x_hatz);
	VEC_ADD<double,DIAG>(u, x_hatz, u);

}
