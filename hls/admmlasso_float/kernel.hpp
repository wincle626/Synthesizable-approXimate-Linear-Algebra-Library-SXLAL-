
#include "math.h"

#define DIAG 16

float invL[DIAG][DIAG];
float invU[DIAG][DIAG];
float oneminusalpha = 0;
float lambdadivrho = 0;
float rho = 0;
float alpha = 0;
float lambda = 0;
float z[DIAG];
float u[DIAG];

template<class T, int M>
void VEC_SUB( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
#pragma HLS RESOURCE variable=V3 core=AddSub
		V3[i] = V1[i] - V2[i];
	}
}

template<class T, int M>
void VEC_SCALAR_MUL( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
#pragma HLS RESOURCE variable=V3 core=Mul
		V3[i] = V1[i] * S;
	}
}

template<class T, int M>
void VEC_ADD( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
#pragma HLS RESOURCE variable=V3 core=AddSub
		V3[i] = V1[i] + V2[i];
	}
}

template<class T, int M, int N>
void MAT_VEC_MUL(T A[M][N], T B[N], T C[M]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		C[i] = 0;
		for ( int j=0; j<N; j++ ){
			T tmp = A[i][j] * B[j];
			C[i] = C[i] + tmp;
		}
	}
}

template<class T, int M>
void VEC_EQ( T V1[M], T V2[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		V2[i] = V1[i];
	}
}

template<class T, int M>
void VEC_SCALAR_SUB( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		V3[i] = V1[i] - S;
	}
}

template<class T, int M>
void VEC_SCALAR_ADD( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] + S;
	}
}

template<class T, int M>
void VEC_MINUS( T V1[M], T V2[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		V2[i] = -V1[i];
	}
}

template<class T, int M>
void VEC_SCALAR_MAX( T V1[M], T S, T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
		V3[i] = V1[i] > S ? V1[i] : S;
	}
}

template<class T, unsigned int M, unsigned int N>
void MAT_TRANS(T Mat[M][N], T MatT[N][M]){
#pragma HLS INLINE off
	for(unsigned int i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		for(unsigned int j=0;j<N;j++){
			MatT[i][j] = Mat[j][i];
		}
	}
}

template<class T, int M, int N, int P>
void MAT_MUL(T A[M][N],
											    T B[N][P],
												T C[M][P]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		for ( int j=0; j<P; j++ ){
#pragma HLS PIPELINE II=1
			C[i][j] = 0;
			T mul, add;
			for ( int k=0; k<N; k++ ){
#pragma HLS RESOURCE variable=mul core=Mul
#pragma HLS RESOURCE variable=add core=AddSub
				mul = A[i][k] * B[k][j];
				add = C[i][j] + mul;
				C[i][j] = add;
			}
		}
	}
}

template<class T, int M, int N>
void IDENDTITY_MAT( T A[M][N] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		for( int j=0; j<N; j++ ){
			if(i==j)
				A[i][j] = 1;
			else
				A[i][j] = 0;
		}
	}
}

template<class T, int M, int N>
void MAT_SCALAR_DOTMUL(T MatA[M][N],
					   T B,
					   T MatC[M][N]){
#pragma HLS INLINE off
	for(int i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		for(int j=0;j<N;j++){
#pragma HLS RESOURCE variable=MatA core=Mul
			MatC[i][j] = MatA[i][j] * B;
		}
	}
}

template<class T, int M, int N>
void MAT_ADD(T A[M][N],
								      T B[M][N],
								      T C[M][N]){
#pragma HLS INLINE off
	for ( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		for ( int j=0; j<N; j++ ){
#pragma HLS RESOURCE variable=C core=AddSub
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

template<class T, int M>
void ZEROS_VEC( T V[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		V[i] = 0;
	}
}

template<class T, int M>
void VEC_NORM( T V1[M], T &S ){
#pragma HLS INLINE off
	S = 0;
	T mul, add;
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
#pragma HLS RESOURCE variable=mul core=Mul
#pragma HLS RESOURCE variable=add core=AddSub
		mul = V1[i] * V1[i];
		add = S + mul;
		S = add;
	}
	T tmp;
	tmp = sqrt(S);
	S = tmp;
}

template<class T, int M>
void VEC_DIV( T V1[M], T V2[M], T V3[M] ){
#pragma HLS INLINE off
	for( int i=0; i<M; i++ ){
#pragma HLS PIPELINE II=1
		V3[i] = V1[i] / V2[i];
	}
}

template<class T, unsigned int M>
void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){
#pragma HLS INLINE off

	// copy the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = 0;
		}
	}

	// decomposition in-place
	for(unsigned int j=0;j<M;j++){
#pragma HLS PIPELINE II=1
		// compute the diagonal element
		T LL = Mat[j][j];
		T mul, sub;
		for(unsigned int k=0;k<j;k++){
#pragma HLS RESOURCE variable=mul core=Mul
#pragma HLS RESOURCE variable=sub core=AddSub
			mul = MatL[j][k] * MatL[j][k];
			sub = LL - mul;
			LL = sub;
//			if(LL<0){
//				exit(-1);
//			}
		}
		MatL[j][j] = sqrt(LL);

		// compute the non-diagonal element
		T inv = 1 / MatL[j][j];
		for(unsigned int i=j+1;i<M;i++){
#pragma HLS PIPELINE II=1
			LL = Mat[i][j];
			T mul, sub;
			for(unsigned int k=0;k<j;k++){
#pragma HLS RESOURCE variable=mul core=Mul
#pragma HLS RESOURCE variable=sub core=AddSub
				mul = MatL[i][k] * MatL[j][k];
				sub = LL - mul;
				LL = sub;
			}
#pragma HLS RESOURCE variable=MatL core=Mul
			MatL[i][j] = LL * inv;
		}
	}

	// transpose L to get U
	for(unsigned int i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		for(unsigned int j=0;j<M;j++){
			MatU[i][j] = MatL[j][i];
		}
	}
}

template<class T, int M, int N>
void QRD_HH(T Mat[M][N],
			T MatQ[M][N],
			T MatR[N][N]){
#pragma HLS INLINE off
	int i,j,k;
	//R=A;
	for(j=0;j<M;j++){
#pragma HLS PIPELINE II=1
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
#pragma HLS PIPELINE II=1
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

template<class T, int M>
void UPTRIANGULARMATINV(T R[M][M],T Ri[M][M]){
#pragma HLS INLINE off
	int i=0,j=0,k=0;
	// R inversion
	for(i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		for(j=0;j<M;j++){
			Ri[i][j]=0;
		}
	}
	for(i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		Ri[i][i]=1/R[i][i];
		for(j=i+1;j<M;j++){
#pragma HLS PIPELINE II=1
			for(k=0;k<=j-1;k++){
				T div = R[k][j] / R[j][j];
				T mul = Ri[i][k] * div;
				T sub = Ri[i][j] - mul;
				Ri[i][j] = sub;
			}
		}
	}
}
template<class T, int M>
void ORTHOGONALMATINV(T Q[M][M],T Qi[M][M]){
#pragma HLS INLINE off
	int i=0,j=0;
	// Q inversion
	for(i=0;i<M;i++){
#pragma HLS PIPELINE II=1
		for(j=0;j<M;j++){
			Qi[i][j] = Q[j][i];
		}
	}
}
template<class T, int M>
void MAT_QRINV(T A[M][M], T B[M][M]){
#pragma HLS INLINE off
	T Q[M][M], R[M][M];
	T Qi[M][M], Ri[M][M];
	QRD_HH<T, M, M>(A, Q, R);
	UPTRIANGULARMATINV<T, M>(R, Ri);
	ORTHOGONALMATINV<T, M>(Q, Qi);
	MAT_MUL<T, M, M, M>(Ri, Qi, B);
}

void kernel(float Atb[DIAG], float x[DIAG]);
