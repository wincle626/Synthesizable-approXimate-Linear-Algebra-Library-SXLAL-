/*
 * vector.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#include <stdlib.h>
#include <time.h>
#include <cmath>

#ifndef HEADER_MATH_VECTOR_HPP_
#define HEADER_MATH_VECTOR_HPP_

// Generate all zero vector
template<class T, int M>
void ZEROS_VEC_R( T V[M] ){
	for( int i=0; i<M; i++ ){
		V[i] = (T) 0;
	}
}

// Generate all one vector
template<class T, int M>
void ONES_VEC_R( T V[M] ){
	for( int i=0; i<M; i++ ){
		V[i] = (T) 1;
	}
}

// Generate random vector
template<class T, int M>
void RND_VEC_R( T V[M] , int FLOAT_SIZE, int INTEGER_SCALE){
	srand (time(NULL));
	for( int i=0; i<M; i++ ){
		double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
		double rndnum = (double) INTEGER_SCALE * rnd ;
		V[i] = (T) rndnum;
	}
}
template<class T, int M>
void RND_VEC_SCALE_R( T V[M], T scale , int FLOAT_SIZE, int INTEGER_SCALE){
	srand (time(NULL));
	for( int i=0; i<M; i++ ){
		double rnd = 2 * (std::rand() % FLOAT_SIZE) / FLOAT_SIZE - 1;
		double rndnum = (double) scale * INTEGER_SCALE * rnd ;
		V[i] = (T) rndnum;
	}
}

// Transfer vector values
template<class T, int M>
void VEC_EQ_R( T V1[M], T V2[M] ){
	for( int i=0; i<M; i++ ){
		V2[i] = V1[i];
	}
}
template<class T1, class T2, int M>
void VEC_EQ_R( T1 V1[M], T2 V2[M] ){
	for( int i=0; i<M; i++ ){
		V2[i] = (T2) V1[i];
	}
}
template<class T1, class T2, int M>
void VEC_SCALAR_EQ_R( T1 S, T2 V2[M] ){
	for( int i=0; i<M; i++ ){
		V2[i] = (T2) S;
	}
}

// reverse the sign of vector elements
template<class T, int M>
void VEC_MINUS_R( T V1[M], T V2[M] ){
	for( int i=0; i<M; i++ ){
		V2[i] = -V1[i];
	}
}

// Vector addition
template<class T, int M>
void VEC_ADD_R( T V1[M], T V2[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] + V2[i];
	}
}
template<class T1, class T2, class T3, int M>
void VEC_ADD_R( T1 V1[M], T2 V2[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] + V2[i]);
	}
}

// Vector subtraction
template<class T, int M>
void VEC_SUB_R( T V1[M], T V2[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] - V2[i];
	}
}
// Vector subtraction
template<class T1, class T2, class T3, int M>
void VEC_SUB_R( T1 V1[M], T2 V2[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] - V2[i]);
	}
}

template<class T, int M>
void VEC_DOT_PROD(T V1[M], T V2[M], T &dotprod){
    dotprod = (T) 0;
    for(int i=0;i<M;i++){
        dotprod += V1[i] * V2[i];
    }
}

// Vector division
template<class T, int M>
void VEC_DIV_R( T V1[M], T V2[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] / V2[i];
	}
}
template<class T1, class T2, class T3, int M>
void VEC_DIV_R( T1 V1[M], T2 V2[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] / V2[i]);
	}
}

// Vector scalar addition
template<class T, int M>
void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] + S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_ADD( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] + S);
	}
}
template<class T, int M>
void VEC_SCALAR_ADD( T S, T V1[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] + S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_ADD( T1 S, T2 V1[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (S + V1[i]);
	}
}

// Vector scalar subtraction
template<class T, int M>
void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] - S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_SUB( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] - S);
	}
}
template<class T, int M>
void VEC_SCALAR_SUB( T S, T V1[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = S - V1[i];
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_SUB( T1 S, T2 V1[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (S - V1[i]);
	}
}

// Vector scalar multiplication
template<class T, int M>
void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] * S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_MUL( T1 S, T2 V1[M], T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] * S);
	}
}
template<class T, int M>
void VEC_SCALAR_MUL( T S, T V1[M], T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] * S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_MUL( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] * S);
	}
}

// Vector scalar division
template<class T, int M>
void VEC_SCALAR_DIV_R( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] / S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_DIV_R( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] / S);
	}
}

template<class T, int M>
void VEC_ABS2_R(T V1[M], int k, int m, T &S){
	S = 0;
	for( int i=k; i<m; i++ ){
		S += V1[i] * V1[i];
	}
}

// Vector l-p norm
template<class T, int M>
void VEC_L1NORM_R( T V1[M], T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += V1[i] > (T) 0 ? V1[i] : -V1[i];
	}
}
template<class T, int M>
void VEC_L2NORM_R( T V1[M], T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += V1[i] * V1[i];
	}
	double tmp = (double) S;
	S = (T) std::sqrt(tmp);
}
template<class T, int M>
void VEC_L3NORM_R( T V1[M], T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += V1[i] * V1[i] * V1[i];
	}
	double tmp = std::pow((double) S, 1/3);
	S = (T) tmp;
}
template<class T, int M>
void VEC_L4NORM_R( T V1[M], T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += V1[i] * V1[i] * V1[i] * V1[i];
	}
	double tmp = std::pow((double) S, 1/4);
	S = (T) tmp;
}
template<class T, int M>
void VEC_LPNORM_R( T V1[M], int p, T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += (T) std::pow((double) V1[i], p);
	}
	double tmp = std::pow((double) S, 1/p);
	S = (T) tmp;
}
template<class T, int M>
void VEC_LINFINORM_R( T V1[M], int p, T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		T tmp = V1[i] > (T) 0 ? V1[i] : -V1[i];
		S = S > tmp ? S : tmp;
	}
}
template<class T, int M>
void VEC_SQRTNORM_R( T V1[M], T &S ){
	S = 0;
	for( int i=0; i<M; i++ ){
		S += V1[i] * V1[i];
	}
}

// Vector scalar compare
template<class T, int M>
void VEC_SCALAR_MIN( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] < S ? V1[i] : S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_MIN( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] < S ? V1[i] : S);
	}
}
template<class T, int M>
void VEC_SCALAR_MAX( T V1[M], T S, T V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] > S ? V1[i] : S;
	}
}
template<class T1, class T2, class T3, int M>
void VEC_SCALAR_MAX( T1 V1[M], T2 S, T3 V3[M] ){
	for( int i=0; i<M; i++ ){
		V3[i] = (T3) (V1[i] > S ? V1[i] : S);
	}
}

// Vector minimum value
template<class T, int M>
void VEC_MIN( T V1[M], T &S){
	S = V1[0];
	for( int i=1; i<M; i++ ){
		S = V1[i] < S ? V1[i] : S;
	}
}

// Vector maximum value
template<class T, int M>
void VEC_MAX( T V1[M], T &S){
	S = V1[0];
	for( int i=1; i<M; i++ ){
		S = V1[i] > S ? V1[i] : S;
	}
}

#endif /* HEADER_MATH_VECTOR_HPP_ */
