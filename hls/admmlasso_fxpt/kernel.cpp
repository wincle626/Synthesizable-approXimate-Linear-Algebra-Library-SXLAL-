#include "kernel.hpp"

void kernel(DATA_IN_T Atb[DIAG], DATA_IN_T x[DIAG]){

	DATA_IN_T q[DIAG];
	DATA_IN_T zold[DIAG];
	DATA_IN_T zminusu[DIAG];
	DATA_IN_T rhozminusu[DIAG];
	DATA_IN_T invLq[DIAG];
	DATA_IN_T alphax[DIAG];
	DATA_IN_T oneminusalphazold[DIAG];
	DATA_IN_T x_hat[DIAG];
	DATA_IN_T x_hatu[DIAG];
	DATA_IN_T x_hatu1[DIAG];
	DATA_IN_T x_hatu2[DIAG];
	DATA_IN_T x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<DATA_IN_T,DIAG>(z, u, zminusu);
	VEC_SCALAR_MUL<DATA_IN_T,DIAG>(zminusu,
			rho, rhozminusu);
	VEC_ADD<DATA_IN_T,DIAG>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
			invL, q, invLq);
	MAT_VEC_MUL<DATA_IN_T,DIAG,DIAG>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<DATA_IN_T,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<DATA_IN_T,DIAG>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<DATA_IN_T,DIAG>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<DATA_IN_T,DIAG>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<DATA_IN_T,DIAG>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<DATA_IN_T,DIAG>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<DATA_IN_T,DIAG>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<DATA_IN_T,DIAG>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<DATA_IN_T,DIAG>(
			x_hatu1, 0, x_hatu1);
	VEC_SCALAR_MAX<DATA_IN_T,DIAG>(
			x_hatu2, 0, x_hatu2);
	VEC_SUB<DATA_IN_T,DIAG>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<DATA_IN_T,DIAG>(x_hat, z, x_hatz);
	VEC_ADD<DATA_IN_T,DIAG>(u, x_hatz, u);

}
