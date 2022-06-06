/*
 *	This is a light bitset fixed point arithmetic
 *	library
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#include "fxpt_bitwise.hpp"

/*
Inputs	|	Outputs
A	B	| Cin	Cout	S
0	0	| 0		0		0
0	0	| 1		0		1
0	1	| 0		0		1
0	1	| 1		1		0
1	0	| 0		0		1
1	0	| 1		1		0
1	1	| 0		1		0
1	1	| 1		1		1
S = A ⊕ B ⊕ Cin
Cout = (A ⋅ B) + (Cin ⋅ (A ⊕ B))
*/
void BIN_DYNBITSETfulladd(bool a, bool b,
			 bool &carrier_in,
			 bool &carrier_out,
			 bool &sum){
	sum = ( a ^ b ) ^ carrier_in;
	carrier_out = ( a & b )
			   | ( carrier_in
			   & ( a ^ b ) );
}

/**
Inputs	|	Outputs
X	Y	| BIN_DYNBITSET	D	Bout
0	0	| 0		0	0
0	0	| 1		1	1
0	1	| 0		1	1
0	1	| 1		0	1
1	0	| 0		1	0
1	0	| 1		0	0
1	1	| 0		0	0
1	1	| 1		1	1
D=X⊕Y⊕BIN_DYNBITSET
Bout=X'BIN_DYNBITSET+X'Y+YBIN_DYNBITSET
 */
void BIN_DYNBITSETfullsub(bool a, bool b,
	 	 	 bool &borrow_in,
		 	 bool &borrow_out,
			 bool &sub){
	sub = ( a ^ b ) ^ borrow_in;
	borrow_out = ( !a & borrow_in )
			   | ( !a & b )
			   | ( b & borrow_in );
}

void nbitfulladd(BIN_DYNBITSET a, BIN_DYNBITSET b,
				 BIN_DYNBITSET &sum){

}

void nbitfullsub(BIN_DYNBITSET a, BIN_DYNBITSET b,
				 BIN_DYNBITSET &sub){

}

void FXPT_DYNBITSET_add(FXPT_DYNBITSET a,
		      FXPT_DYNBITSET b,
			  FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_sub(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET b,
			  FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_mul(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET b,
			  FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_div(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET b,
			  FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_sqrt(FXPT_DYNBITSET a,
			   FXPT_DYNBITSET &b){

}

void FXPT_DYNBITSET_root(FXPT_DYNBITSET a,
		  	   FXPT_DYNBITSET b,
			   FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_abs(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET &b){

}

void FXPT_DYNBITSET_exp(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET &c){

}

void FXPT_DYNBITSET_pow2(FXPT_DYNBITSET a,
			   FXPT_DYNBITSET &b){

}

void FXPT_DYNBITSET_pow(FXPT_DYNBITSET a,
			  FXPT_DYNBITSET b,
			  FXPT_DYNBITSET &c){

}

