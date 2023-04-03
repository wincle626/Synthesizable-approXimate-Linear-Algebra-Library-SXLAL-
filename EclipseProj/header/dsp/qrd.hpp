/*
 * qrd.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_DSP_QRD_HPP_
#define HEADER_DSP_QRD_HPP_

#include <stdlib.h>
#include <cmath>
//#include "mathfunc.hpp"
#include "matrix.hpp"
#include "vector.hpp"

// Gram-Schmidt QR Decomposition
template<class T, int M, int N>
void QRD_GS(T MatA[M][N],
		T MatQ[M][N],
		T MatR[N][N]){
	// Phase 1: initialisation
	int i = 0, j = 0;
	for(i=0;i<N;i++){
		for(j=0;j<M;j++){
			MatQ[j][i] = 0;
		}
		for(j=0;j<N;j++){
			MatR[i][j] = 0;
		}
	}

	// Phase 2: get the Q&R
	T norm[N];
	T V1[M], V2[M], V3[M], DotProd;
	for(j=0;j<N;j++){
//        v = X(:,j);
		COL_OF_MATRIX<T, M, N>(MatA, j, V1);
		for(i=0;i<j;i++){
//            R(i,j) = Q(:,i)'*X(:,j);
			COL_OF_MATRIX<T, M, N>(MatA, j, V2);
			COL_OF_MATRIX<T, M, N>(MatQ, i, V3);
			VEC_DOT_PROD<T, M>(V2, V3, DotProd);
			MatR[i][j] = DotProd;
//            v = v - R(i,j)*Q(:,i);
			VEC_SCALAR_MUL<T, M>(V3, DotProd, V3);
			VEC_SUB_R<T, M>(V1, V3, V1);
		}
//        R(j,j) = norm(v);
		VEC_L2NORM_R<T, M>(V1, norm[j]);
		MatR[j][j] = norm[j];
//        Q(:,j) = v/R(j,j);
		VEC_SCALAR_DIV_R<T, M>(V1, norm[j], V1);
		COL_INTO_MATRIX<T, M, N>(V1, MatQ, j);
	}

}

// Gram-Schmidt QR Decomposition
template<class T, int M, int N>
void QRD_MGS(T MatA[M][N],
		T MatQ[M][N],
		T MatR[N][N]){
	// Phase 1: initialisation
	int i = 0, k = 0;
	for(i=0;i<N;i++){
		for(k=0;k<M;k++){
			MatQ[k][i] = 0;
		}
		for(k=0;k<N;k++){
			MatR[i][k] = 0;
		}
	}

	// Phase 2: get the Q&R
	T norm[N];
	T V1[M], V2[M], DotProd;
	for(k=0;k<N;k++){
//        Q(:,k) = X(:,k);
		COL_OF_MATRIX<T, M, N>(MatA, k, V1);
		COL_INTO_MATRIX<T, M, N>(V1, MatQ, k);
		for(i=0;i<k;i++){
//            R(i,k) = Q(:,i)'*Q(:,k);
			COL_OF_MATRIX<T, M, N>(MatQ, i, V2);
			VEC_DOT_PROD<T, M>(V2, V1, DotProd);
			MatR[i][k] = DotProd;
//            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
			VEC_SCALAR_MUL<T, M>(V2, DotProd, V2);
			VEC_SUB_R<T, M>(V1, V2, V1);
		}
//        R(k,k) = norm(Q(:,k))';
		VEC_L2NORM_R<T, M>(V1, norm[k]);
		MatR[k][k] = norm[k];
//        Q(:,k) = Q(:,k)/R(k,k);
		VEC_SCALAR_DIV_R<T, M>(V1, norm[k], V1);
		COL_INTO_MATRIX<T, M, N>(V1, MatQ, k);
	}

}

// Gram-Schmidt QR Decomposition
template<class T, int M, int N>
void QRD_GS1(T Mat[M][N],
		T MatQ[M][N],
		T MatR[N][N]){
	// Initialisation
	int i = 0,j = 0,k = 0;
	for(i=0;i<N;i++){
		for(j=0;j<M;j++){
			MatQ[j][i] = Mat[j][i];
		}
		for(j=0;j<N;j++){
			MatR[i][j] = 0;
		}
	}

	// Phase 1: get the norm
	float norm[N];
	for(i=0;i<N;i++){
		norm[i] = 0;
		for(j=0;j<M;j++){
			float mul = Mat[j][i] * Mat[j][i];
			norm[i] = norm[i] + mul;
		}
	}

	// Phase 2: get the Q&R
	for(i=0;i<N;i++){
		// derive R
		MatR[i][i] = (T) std::sqrt(norm[i]);
		for(k=0;k<M;k++){
			MatQ[k][i]=MatQ[k][i]/MatR[i][i];
		}
		for(j=i+1;j<M;j++){
			for(k=0;k<M;k++){
				// update R
				MatR[i][j] = MatR[i][j] + MatQ[k][i] * MatQ[k][j];
			}
			for(k=0;k<M;k++){
				// update Q
				MatQ[k][j] = MatQ[k][j] - MatQ[k][i] * MatR[i][j];
			}
			// update norm
			norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
		}
	}
}

// Modified Gram-Schmidt QR Decomposition
template<class T, int M, int N>
void QRD_MGS1(T Mat[M][N],
			 T MatQ[M][N],
			 T MatR[N][N]){
	// Initialisation
	int i = 0,j = 0,k = 0;
	for(i=0;i<N;i++){
		for(j=0;j<M;j++){
			MatQ[j][i] = Mat[j][i];
		}
		for(j=0;j<N;j++){
			MatR[i][j] = 0;
		}
	}

	// Phase 1: get the norm
	float norm[N];
	for(i=0;i<N;i++){
		norm[i] = 0;
		for(j=0;j<M;j++){
			float mul = Mat[j][i] * Mat[j][i];
			norm[i] = norm[i] + mul;
		}
	}

	// Phase 2: get the Q&R
	for(i=0;i<N;i++){
		// derive R
		MatR[i][i] = (T) std::sqrt(norm[i]);
		for(k=0;k<M;k++){
			MatQ[k][i]=MatQ[k][i]/MatR[i][i];
		}
		for(j=i+1;j<M;j++){
			float tmp;
			for(k=0;k<M;k++){
				// update R
				MatR[i][j] = MatR[i][j] + MatQ[k][i] * MatQ[k][j];
			}
			for(k=0;k<M;k++){
				// update Q
				MatQ[k][j] = MatQ[k][j] - MatQ[k][i] * MatR[i][j];
			}
		// update norm: no update for QR_MGS
			// norm[j] = norm[j] - MatR[i][j] * MatR[i][j];
		}
	}

}

// Householder QR Decomposition
template<class T, int M, int N>
void QRD_HH(T Mat[M][N],
			T MatQ[M][M],
			T MatR[M][N]){
	int i,j,k;
	T MatH[M][M];
	//R=A;
	for(j=0;j<M;j++){
		for(i=0;i<N;i++){
			MatR[j][i] = Mat[j][i];
		}
		for(i=0;i<M;i++){
			//H=eye(m);
			//Q=eye(m);
			if(i==j){
				MatQ[j][i] = (T) 1;
				MatH[j][i] = (T) 1;
			}
			else{
				MatQ[j][i] = (T) 0;
				MatH[j][i] = (T) 0;
			}
		}
	}


	T g, s;
	T x[M], v[M];
	for(k=0;k<N;k++){
//        x = R(k:m,k);
		ZEROS_VEC_R<T, M>(x);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
//		printvector<T, M>(x);
//        v = x;
		VEC_EQ_R<T, M>(x, v);
//        v(1) = v(1) + sign(x(1)) * norm(x);
		VEC_L2NORM_R<T, M>(x, g);
		T signx1 = (T)(x[k]>0?1:-1);
		v[k] = x[k] + signx1*g;
//        v = v ./ v(1);
		VEC_SCALAR_DIV_R<T, M>(v, v[k], v);
//		printvector<T, M>(v);
//        h = eye(m-k+1) - 2 .* (v*v') ./ (v'*v);
//        H(k:m, k:m) = h;
		VEC_ABS2_R<T, M>(v, k, M, s);
		s=-2/s;
		for(i=k;i<M;i++){
			for(j=k;j<M;j++){
				if(i==j){
					MatH[i][j] = 1 + s * v[i] * v[j];
				}else{
					MatH[i][j] = s * v[i] * v[j];
				}
			}
		}
//		printmatrix<T, M, M>(MatH);
		T tmp1[M][N], tmp2[M][M], MatHT[M][M];
		MAT_MUL<T, M, M, N>(MatH, MatR, tmp1);
		MAT_TRANS<T, M, M>(MatH, MatHT);
		MAT_MUL<T, M, M, M>(MatQ, MatHT, tmp2);
		MAT_EQ<T, M, N>(tmp1, MatR);
//		printmatrix<T, M, N>(MatR);
		MAT_EQ<T, M, M>(tmp2, MatQ);
//		printmatrix<T, M, M>(MatQ);

		for(i=k;i<M;i++){
			for(j=k;j<M;j++){
				if(i==j)
					MatH[i][j] = 1;
				else
					MatH[i][j] = 0;
			}
		}
	}

}

// Householder QR Decomposition
template<class T, int M, int N>
void QRD_HH1(T Mat[M][N],
			T MatQ[M][N],
			T MatR[N][N]){
	int i,j,k;
	//R=A;
	for(j=0;j<M;j++){
		for(i=0;i<N;i++){
			MatR[j][i] = Mat[j][i];
			//Q=eye(m);
			if(i==j)
				MatQ[j][i] = (T) 1;
			else
				MatQ[j][i] = (T) 0;
		}
	}

	T g, s;
	T x[M], v[M], w[M], u[N];
	T tmp1[M][N], tmp2[M][N];
	for(k=0;k<M-1;k++){
		// x=zeros(m,1);
		for(j=0;j<M;j++){
			x[j] = (T) 0;
		}
		ZEROS_VEC_R<T, M>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_L2NORM_R<T, M>(x, g);
		// v=x;
		VEC_EQ_R<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = x[k] + g;
		//s=norm(v);
		VEC_L2NORM_R<T, M>(v, s);
		if(s!=0){
			// w=v/s;
			VEC_SCALAR_DIV_R<T, M>(v, s, w);
			// u=2*R'*w;
			for(i=0;i<N;i++){
				u[i] = 0;
				for(j=0;j<M;j++){
					u[i] += 2 * MatR[j][i] * w[j];
				}
			}
			// R=R-w*u'; %Product HR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatR[j][i] = MatR[j][i] - w[j] * u[i];
				}
			}
			// Q=Q-2*Q*w*w'; %Product QR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					tmp1[j][i] = w[j] * w[i];
				}
			}
			MAT_MUL<T, M, N, N>(MatQ, tmp1, tmp2);
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatQ[j][i] = MatQ[j][i] - 2 * tmp2[j][i];
				}
			}
		}
	}
}
template<class T, int M, int N>
void QRD_HH1(T **Mat,
			T **MatQ,
			T **MatR){
	int i,j,k;
	//R=A;
	for(j=0;j<M;j++){
		for(i=0;i<N;i++){
			MatR[j][i] = Mat[j][i];
			//Q=eye(m);
			if(i==j)
				MatQ[j][i] = (T) 1;
			else
				MatQ[j][i] = (T) 0;
		}
	}

	T g, s;
	T *x, *v, *w, *u;
	x = (T*) malloc(sizeof(T)*M);
	v = (T*) malloc(sizeof(T)*M);
	w = (T*) malloc(sizeof(T)*M);
	u = (T*) malloc(sizeof(T)*N);
	T **tmp1, **tmp2;
	tmp1 = (T**) malloc(sizeof(T*)*M);
	tmp2 = (T**) malloc(sizeof(T*)*M);
	for(int i=0;i<N;i++){
		tmp1[i] = (T*) malloc(sizeof(T)*N);
		tmp2[i] = (T*) malloc(sizeof(T)*N);
	}
	for(k=0;k<M-1;k++){
		// x=zeros(m,1);
		for(j=0;j<M;j++){
			x[j] = (T) 0;
		}
		ZEROS_VEC_R<T, M>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_L2NORM_R<T, M>(x, g);
		// v=x;
		VEC_EQ_R<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = x[k] + g;
		//s=norm(v);
		VEC_L2NORM_R<T, M>(v, s);
		if(s!=0){
			// w=v/s;
			VEC_DIV_R<T, M>(v, s, w);
			// u=2*R'*w;
			for(i=0;i<N;i++){
				u[i] = 0;
				for(j=0;j<M;j++){
					u[i] += 2 * MatR[j][i] * w[j];
				}
			}
			// R=R-w*u'; %Product HR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatR[j][i] = MatR[j][i] - w[j] * u[i];
				}
			}
			// Q=Q-2*Q*w*w'; %Product QR
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					tmp1[j][i] = w[j] * w[i];
				}
			}
			MAT_MUL<T, M, N, N>(MatQ, tmp1, tmp2);
			for(j=0;j<M;j++){
				for(i=0;i<N;i++){
					MatQ[j][i] = MatQ[j][i] - 2 * tmp2[j][i];
				}
			}
		}
	}
}


// Given Rotation QR Decomposition
template<class T>
void GR(T a, T b, T &c, T&s){
	if(b==0){
		c = 1;
		s = 0;
	}else{
		T tmp1 = a/b;
		T tmp2 = b/a;
		T r = (T) (b>=a||-b<-a)?tmp1:tmp2;
		T tmp3 = 1/std::sqrt(1+r*r);
		T tmp4 = tmp1*tmp3;
		T tmp5 = tmp2*tmp3;
		c = (T) (b>=a||-b<-a)? tmp4: tmp3;
		s = (T) (b>=a||-b<-a)? tmp3: tmp5;
	}
}
template<class T, int M, int N>
void QRD_GR(T MatA[M][N],
			T MatQ[M][M],
			T MatR[M][N]){
	int i, j, k;
	for(i=0;i<M;i++){
		for(j=0;j<M;j++){
			MatQ[j][i] = 0;
		}
		for(j=0;j<N;j++){
			MatR[i][j] = MatA[i][j];
		}
	}
	for(i=0;i<M;i++) {
		MatQ[i][i] = 1;
	}

	for(j=0;j<N;j++){
		for(i=M-1;i>j;i--){
//            G([i-1, i],[i-1, i]) = [c -s; s c];
			T c, s;
			GR(MatR[i-1][j], MatR[i][j], c, s);
//            Q = Q*G;
			T Qpart[M][2];
			for(k=0;k<M;k++){
				Qpart[k][0] = c*MatQ[k][i-1] + s*MatQ[k][i];
				Qpart[k][1] = -s*MatQ[k][i-1] + c*MatQ[k][i];
			}
			for(k=0;k<M;k++){
				MatQ[k][i-1] = Qpart[k][0];
				MatQ[k][i] = Qpart[k][1];
			}
//            R = G'*R;
			T Rpart[2][N];
			for(k=0;k<N;k++){
				Rpart[0][k] = c*MatR[i-1][k] + s*MatR[i][k];
				Rpart[1][k] = -s*MatR[i-1][k] + c*MatR[i][k];
			}
			for(k=0;k<N;k++){
				MatR[i-1][k] = Rpart[0][k];
				MatR[i][k] = Rpart[1][k];
			}
		}
	}
}


#endif /* HEADER_DSP_QRD_HPP_ */
