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

// Gram-Schmidt QR Decomposition
template<class T, int M, int N>
void QRD_GS(T Mat[M][N],
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
void QRD_MGS(T Mat[M][N],
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
		ZEROS_VEC<T, M>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_NORM<T, M>(x, g);
		// v=x;
		VEC_EQ<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = x[k] + g;
		//s=norm(v);
		VEC_NORM<T, M>(v, s);
		if(s!=0){
			// w=v/s;
			VEC_DIV<T, M>(v, s, w);
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
void QRD_HH(T **Mat,
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
		ZEROS_VEC<T, M>(x);
		//x(k:m,1)=R(k:m,k);
		for(j=k;j<M;j++){
			x[j] = MatR[j][k];
		}
		//g=norm(x);
		VEC_NORM<T, M>(x, g);
		// v=x;
		VEC_EQ<T, M>(x, v);
		// v(k)=x(k)+g;
		v[k] = x[k] + g;
		//s=norm(v);
		VEC_NORM<T, M>(v, s);
		if(s!=0){
			// w=v/s;
			VEC_DIV<T, M>(v, s, w);
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



#endif /* HEADER_DSP_QRD_HPP_ */
