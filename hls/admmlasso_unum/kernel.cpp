#include "kernel.hpp"

#if SOFT_POSIT_PRECISION==0

void kernel(posit32_t Atb[DIAG], posit32_t x[DIAG]){
//#pragma HLS INTERFACE ap_memory port=x
//#pragma HLS INTERFACE ap_memory port=Atb
//#pragma HLS RESOURCE variable=z core=RAM_2P
//#pragma HLS RESOURCE variable=u core=RAM_2P
//#pragma HLS RESOURCE variable=invL core=ROM_1P
//#pragma HLS RESOURCE variable=invU core=ROM_1P

	posit32_t q[DIAG];
	posit32_t zold[DIAG];
	posit32_t zminusu[DIAG];
	posit32_t rhozminusu[DIAG];
	posit32_t invLq[DIAG];
	posit32_t alphax[DIAG];
	posit32_t oneminusalphazold[DIAG];
	posit32_t x_hat[DIAG];
	posit32_t x_hatu[DIAG];
	posit32_t x_hatu1[DIAG];
	posit32_t x_hatu2[DIAG];
	posit32_t x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<posit32_t,DIAG,p32_sub>(z, u, zminusu);
	VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(zminusu,
			rho, rhozminusu);
	VEC_ADD<posit32_t,DIAG,p32_add>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(
			invL, q, invLq);
	MAT_VEC_MUL<posit32_t,DIAG,DIAG,convertDoubleToP32,p32_mul,p32_add>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<posit32_t,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<posit32_t,DIAG,p32_mul>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<posit32_t,DIAG,p32_add>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<posit32_t,DIAG,p32_add>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<posit32_t,DIAG,p32_sub>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<posit32_t,DIAG,p32_add>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<posit32_t,DIAG,convertDoubleToP32,p32_mul>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<posit32_t,DIAG,p32_lt>(
			x_hatu1, convertDoubleToP32(0), x_hatu1);
	VEC_SCALAR_MAX<posit32_t,DIAG,p32_lt>(
			x_hatu2, convertDoubleToP32(0), x_hatu2);
	VEC_SUB<posit32_t,DIAG,p32_sub>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<posit32_t,DIAG,p32_sub>(x_hat, z, x_hatz);
	VEC_ADD<posit32_t,DIAG,p32_add>(u, x_hatz, u);

}

#elif SOFT_POSIT_PRECISION==1

void kernel(posit16_t Atb[DIAG], posit16_t x[DIAG]){
//#pragma HLS INTERFACE ap_memory port=x
//#pragma HLS INTERFACE ap_memory port=Atb
//#pragma HLS RESOURCE variable=z core=RAM_2P
//#pragma HLS RESOURCE variable=u core=RAM_2P
//#pragma HLS RESOURCE variable=invL core=ROM_1P
//#pragma HLS RESOURCE variable=invU core=ROM_1P

	posit16_t q[DIAG];
	posit16_t zold[DIAG];
	posit16_t zminusu[DIAG];
	posit16_t rhozminusu[DIAG];
	posit16_t invLq[DIAG];
	posit16_t alphax[DIAG];
	posit16_t oneminusalphazold[DIAG];
	posit16_t x_hat[DIAG];
	posit16_t x_hatu[DIAG];
	posit16_t x_hatu1[DIAG];
	posit16_t x_hatu2[DIAG];
	posit16_t x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<posit16_t,DIAG,p16_sub>(z, u, zminusu);
	VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(zminusu,
			rho, rhozminusu);
	VEC_ADD<posit16_t,DIAG,p16_add>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<posit16_t,DIAG,DIAG,convertDoubleToP16,p16_mul,p16_add>(
			invL, q, invLq);
	MAT_VEC_MUL<posit16_t,DIAG,DIAG,convertDoubleToP16,p16_mul,p16_add>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<posit16_t,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<posit16_t,DIAG,p16_mul>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<posit16_t,DIAG,p16_add>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<posit16_t,DIAG,p16_add>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<posit16_t,DIAG,p16_sub>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<posit16_t,DIAG,p16_add>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<posit16_t,DIAG,convertDoubleToP16,p16_mul>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<posit16_t,DIAG,p16_lt>(
			x_hatu1, convertDoubleToP16(0), x_hatu1);
	VEC_SCALAR_MAX<posit16_t,DIAG,p16_lt>(
			x_hatu2, convertDoubleToP16(0), x_hatu2);
	VEC_SUB<posit16_t,DIAG,p16_sub>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<posit16_t,DIAG,p16_sub>(x_hat, z, x_hatz);
	VEC_ADD<posit16_t,DIAG,p16_add>(u, x_hatz, u);

}

#elif SOFT_POSIT_PRECISION==2

void kernel(posit_2_t Atb[DIAG], posit_2_t x[DIAG]){
//#pragma HLS INTERFACE ap_memory port=x
//#pragma HLS INTERFACE ap_memory port=Atb
//#pragma HLS RESOURCE variable=z core=RAM_2P
//#pragma HLS RESOURCE variable=u core=RAM_2P
//#pragma HLS RESOURCE variable=invL core=ROM_1P
//#pragma HLS RESOURCE variable=invU core=ROM_1P

	posit_2_t q[DIAG];
	posit_2_t zold[DIAG];
	posit_2_t zminusu[DIAG];
	posit_2_t rhozminusu[DIAG];
	posit_2_t invLq[DIAG];
	posit_2_t alphax[DIAG];
	posit_2_t oneminusalphazold[DIAG];
	posit_2_t x_hat[DIAG];
	posit_2_t x_hatu[DIAG];
	posit_2_t x_hatu1[DIAG];
	posit_2_t x_hatu2[DIAG];
	posit_2_t x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<posit_2_t,DIAG,p8_sub>(z, u, zminusu);
	VEC_SCALAR_MUL<posit_2_t,DIAG,p8_mul>(zminusu,
			rho, rhozminusu);
	VEC_ADD<posit_2_t,DIAG,p8_add>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<posit_2_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(
			invL, q, invLq);
	MAT_VEC_MUL<posit_2_t,DIAG,DIAG,convertDoubleToP8,p8_mul,p8_add>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<posit_2_t,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<posit_2_t,DIAG,p8_mul>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<posit_2_t,DIAG,p8_mul>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<posit_2_t,DIAG,p8_add>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<posit_2_t,DIAG,p8_add>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<posit_2_t,DIAG,p8_sub>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<posit_2_t,DIAG,p8_add>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<posit_2_t,DIAG,convertDoubleToP8,p8_mul>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<posit_2_t,DIAG,p8_lt>(
			x_hatu1, convertDoubleToP8(0), x_hatu1);
	VEC_SCALAR_MAX<posit_2_t,DIAG,p8_lt>(
			x_hatu2, convertDoubleToP8(0), x_hatu2);
	VEC_SUB<posit_2_t,DIAG,p8_sub>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<posit_2_t,DIAG,p8_sub>(x_hat, z, x_hatz);
	VEC_ADD<posit_2_t,DIAG,p8_add>(u, x_hatz, u);

}

#elif SOFT_POSIT_PRECISION==3

void kernel(posit_2_t Atb[DIAG], posit_2_t x[DIAG]){
//#pragma HLS INTERFACE ap_memory port=x
//#pragma HLS INTERFACE ap_memory port=Atb
//#pragma HLS RESOURCE variable=z core=RAM_2P
//#pragma HLS RESOURCE variable=u core=RAM_2P
//#pragma HLS RESOURCE variable=invL core=ROM_1P
//#pragma HLS RESOURCE variable=invU core=ROM_1P

	posit_2_t q[DIAG];
	posit_2_t zold[DIAG];
	posit_2_t zminusu[DIAG];
	posit_2_t rhozminusu[DIAG];
	posit_2_t invLq[DIAG];
	posit_2_t alphax[DIAG];
	posit_2_t oneminusalphazold[DIAG];
	posit_2_t x_hat[DIAG];
	posit_2_t x_hatu[DIAG];
	posit_2_t x_hatu1[DIAG];
	posit_2_t x_hatu2[DIAG];
	posit_2_t x_hatz[DIAG];

	// q = Atb + rho*(z - u);
	VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(z, u, zminusu);
	VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(zminusu,
			rho, rhozminusu);
	VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(Atb,
			rhozminusu, q);
	// x = U \ (L \ q);
	MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
			invL, q, invLq);
	MAT_VEC_MUL<posit_2_t,posit_1_t,DIAG,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,pX2_add,TOTALBITS>(
			invU, invLq, x);
	// zold = z
	VEC_EQ<posit_2_t,DIAG>(z, zold);
	//  x_hat = alpha*x + (1 - alpha)*zold;
	VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
			x, alpha, alphax);
	VEC_SCALAR_MUL<posit_2_t,DIAG,pX2_mul,TOTALBITS>(
			zold, oneminusalpha, oneminusalphazold);
	VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
			alphax, oneminusalphazold, x_hat);
	// z = max( 0, x - kappa ) - max( 0, -x - kappa );
	VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(x_hat, u, x_hatu);
	VEC_SCALAR_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(
			x_hatu, lambdadivrho, x_hatu1);
	VEC_SCALAR_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(
			x_hatu, lambdadivrho, x_hatu2);
	VEC_MINUS<posit_2_t,posit_1_t,DIAG,convertDoubleToPX1,pX1_to_pX2,pX2_mul,TOTALBITS>(
			x_hatu2, x_hatu2);
	VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
			x_hatu1, pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS),TOTALBITS), x_hatu1);
	VEC_SCALAR_MAX<posit_2_t,DIAG,pX2_lt>(
			x_hatu2, pX1_to_pX2(convertDoubleToPX1(0,TOTALBITS),TOTALBITS), x_hatu2);
	VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hatu1, x_hatu2, z);
	// u = u + (x_hat - z);
	VEC_SUB<posit_2_t,DIAG,pX2_sub,TOTALBITS>(x_hat, z, x_hatz);
	VEC_ADD<posit_2_t,DIAG,pX2_add,TOTALBITS>(u, x_hatz, u);

}
#endif
