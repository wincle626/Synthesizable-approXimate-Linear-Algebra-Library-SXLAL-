/*
 * tool.hpp
 *
 *  Created on: Mar 25, 2023
 *      Author: ubuntu
 */

#ifndef SRC_TOOL_HPP_
#define SRC_TOOL_HPP_

#include <cstdio>

void printmatrix(double *A, int row, int col){
	printf("\n");
	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			printf("%f\t", A[i*col+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printmatrix(double **A, int row, int col){
	printf("\n");
	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			printf("%f\t", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

template<class T, int M, int N>
void printmatrix(T A[M][N]){
	printf("\n");
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			printf("%f\t", (float)A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printvector(double *V, int row){
	printf("\n");
	for(int i=0;i<row;i++){
		printf("%f\t", V[i]);
	}
	printf("\n");
}

template<class T, int M>
void printvector(T V[M]){
	printf("\n");
	for(int i=0;i<M;i++){
		printf("%f\t", V[i]);
	}
	printf("\n");
}



#endif /* SRC_TOOL_HPP_ */
