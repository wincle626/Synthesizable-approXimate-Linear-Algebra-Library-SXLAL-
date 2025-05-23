/*
 *	This is C++ implementation Ali's admm
 *	Author: Yun Wu
 *	Created by: 2019-07-18
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_ADMM_HPP_
#define SRC_ADMM_HPP_

#include "eigen_algebra.hpp"
#include "fpt_algebra.hpp"
#include "fft.hpp"

#define BETA 0.001
#define MAX_ITR 4000
#define tol 0.00001
#define RHO 0.1
#define RHO2 RHO
#define DIM_SIZE 101
#define SQDIM_SIZE 256
#define BITSIZE 8

//#define DEBUG_ITER
#define DEBUG_RESULT

template<typename T, int M, int N>
void D_eigen(T U, T Dux, T Duy){
	Eigen_Algebra eigenalgebraobj;
	eigenalgebraobj.MAT_DIFF<T, M, N, 2, 1>(U,Dux,1);
	for(int i=0;i<M;i++){
		Dux(i,N-1) = U(i,0) - U(i,N-1);
	}
	eigenalgebraobj.MAT_DIFF<T, M, N, 1, 1>(U,Duy,1);
	for(int i=0;i<N;i++){
		Duy(M-1,i) = U(0,i) - U(M-1,i);
	}
}

template<typename T=double, int M, int N>
void Dd(T U[M][N],
	   T Dux[M][N],
	   T Duy[M][N]){

	Float_Point_Algebra fptalgebraobj;
	fptalgebraobj.MAT_DIFF<T, M, N, 2, 1>(U,Dux,1);
	for(int i=0;i<M;i++)
		Dux[i][N-1] = U[i][0] - U[i][N-1];
	fptalgebraobj.MAT_DIFF<T, M, N, 1, 1>(U,Duy,1);
	for(int i=0;i<N;i++)
		Duy[M-1][i] = U[0][i] - U[M-1][i];

}

template<typename T=float, int M, int N>
void Df(T U[M][N],
	   T Dux[M][N],
	   T Duy[M][N]){

	Float_Point_Algebra fptalgebraobj;
	fptalgebraobj.MAT_DIFF<T, M, N, 2, 1>(U,Dux,1);
	for(int i=0;i<M;i++)
		Dux[i][N-1] = U[i][0] - U[i][N-1];
	fptalgebraobj.MAT_DIFF<T, M, N, 1, 1>(U,Duy,1);
	for(int i=0;i<N;i++)
		Duy[M-1][i] = U[0][i] - U[M-1][i];

}

template<typename T, int M, int N>
void Dt_eigne(T X, T Y, T DtXY){

	Eigen_Algebra eigenalgebraobj;
	T dtxy_tmp1;
	T dtxy_tmp2;
	eigenalgebraobj.MAT_DIFF<T, M, N, 2, 1>(X,dtxy_tmp1,-1);
	for(int i=0;i<M;i++)
		for(int j=N-1;j>0;j--)
			dtxy_tmp1(i,j) = dtxy_tmp1(i,j-1);
	for(int i=0;i<M;i++)
		dtxy_tmp1(i,0) = X(i,N-1) - X(i,0);
	eigenalgebraobj.MAT_DIFF<T, M, N, 1, 1>(Y,dtxy_tmp2,-1);
	for(int i=0;i<N;i++)
		for(int j=M-1;j>0;j--)
			dtxy_tmp2(j,i) = dtxy_tmp2(j-1,i);
	for(int i=0;i<N;i++)
		dtxy_tmp2(0,i) = Y(N-1,i) - Y(0,i);
	eigenalgebraobj.MAT_ADD<T>(dtxy_tmp1,dtxy_tmp2,DtXY);
}

template<typename T=double, int M, int N>
void Dtd(T X[M][N],
	    T Y[M][N],
	    T DtXY[M][N]){

	Float_Point_Algebra fptalgebraobj;
	T dtxy_tmp1[M][N];
	T dtxy_tmp2[M][N];
	fptalgebraobj.MAT_DIFF<T, M, N, 2, 1>(X,dtxy_tmp1,-1);
	for(int i=0;i<M;i++)
		for(int j=N-1;j>0;j--)
			dtxy_tmp1[i][j] = dtxy_tmp1[i][j-1];
	for(int i=0;i<M;i++)
		dtxy_tmp1[i][0] = X[i][N-1] - X[i][0];
	fptalgebraobj.MAT_DIFF<T, M, N, 1, 1>(Y,dtxy_tmp2,-1);
	for(int i=0;i<N;i++)
		for(int j=M-1;j>0;j--)
			dtxy_tmp2[j][i] = dtxy_tmp2[j-1][i];
	for(int i=0;i<N;i++)
		dtxy_tmp2[0][i] = Y[N-1][i] - Y[0][i];
	fptalgebraobj.MAT_ADD<T, M, N>(dtxy_tmp1,dtxy_tmp2,DtXY);

}

template<typename T=float, int M, int N>
void Dtf(T X[M][N],
	    T Y[M][N],
	    T DtXY[M][N]){

	Float_Point_Algebra fptalgebraobj;
	T dtxy_tmp1[M][N];
	T dtxy_tmp2[M][N];
	fptalgebraobj.MAT_DIFF<T, M, N, 2, 1>(X,dtxy_tmp1,-1);
	for(int i=0;i<M;i++)
		for(int j=N-1;j>0;j--)
			dtxy_tmp1[i][j] = dtxy_tmp1[i][j-1];
	for(int i=0;i<M;i++)
		dtxy_tmp1[i][0] = X[i][N-1] - X[i][0];
	fptalgebraobj.MAT_DIFF<T, M, N, 1, 1>(Y,dtxy_tmp2,-1);
	for(int i=0;i<N;i++)
		for(int j=M-1;j>0;j--)
			dtxy_tmp2[j][i] = dtxy_tmp2[j-1][i];
	for(int i=0;i<N;i++)
		dtxy_tmp2[0][i] = Y[N-1][i] - Y[0][i];
	fptalgebraobj.MAT_ADD<T, M, N>(dtxy_tmp1,dtxy_tmp2,DtXY);

}


void ADMM_DOUBLE();
void ADMM_FLOAT();

#endif /* SRC_ADMM_HPP_ */
