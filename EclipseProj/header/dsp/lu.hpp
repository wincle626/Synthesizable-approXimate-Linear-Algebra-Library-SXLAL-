/*
 * lu.hpp
 *
 *  Created on: Jun 4, 2022
 *      Author: yunwu
 */

#ifndef HEADER_DSP_LU_HPP_
#define HEADER_DSP_LU_HPP_

#include <stdlib.h>
#include <cmath>

// Cholesky–Banachiewicz (row based)
// and Cholesky–Crout (column based)
// LU decomposition
template<class T, unsigned int M>
void LU_CHOLBANACHROUT(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

	// copy the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = 0;
		}
	}

	// decomposition in-place
	for(unsigned int j=0;j<M;j++){
		// compute the diagonal element
		T LL = Mat[j][j];
//			std::cout << j << "th diag update:" << std::endl;
//			std::cout << LL;
		for(unsigned int k=0;k<j;k++){
//				std::cout << "-"
//						  << MatL[j][k]
//					      << "^2";
			LL -= MatL[j][k] * MatL[j][k];
			if(LL<0){
//					std::cout << "=" LL << std::endl << std::endl;
//					std::cout << "something is wrong at ["
//							  << j << "]["
//							  << k << "]=";
				exit(-1);
			}
		}
		MatL[j][j] = (T) std::sqrt((double)LL);
//			std::cout << "=" << MatL[j][j] << std::endl << std::endl;

		// compute the non-diagonal element
//			std::cout << j << "th non-diag update:" << std::endl;
		T inv = 1 / MatL[j][j];
		for(unsigned int i=j+1;i<M;i++){
			LL = Mat[i][j];
//				std::cout << LL;
			for(unsigned int k=0;k<j;k++){
//					std::cout << "-"
//							  << MatL[i][k]
//						      << "x"
//							  << MatL[j][k];
				LL -= MatL[i][k] * MatL[j][k];
			}
			MatL[i][j] = LL * inv;
//				std::cout << "x" << inv << "=" << MatL[i][j] << std::endl;
		}
//			std::cout << std::endl;

//			std::cout << "L Matrix: " << std::endl;;
//			for(int i=0;i<4;i++){
//				for(int j=0;j<4;j++){
//					std::cout << MatL[i][j] << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
	}

	// transpose L to get U
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatU[i][j] = MatL[j][i];
		}
	}
}
template<class T, unsigned int M>
void LU_CHOLBANACHROUT(T **Mat, T **MatL, T **MatU){

	// copy the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = 0;
		}
	}

	// decomposition in-place
	for(unsigned int j=0;j<M;j++){
		// compute the diagonal element
		T LL = Mat[j][j];
//			std::cout << j << "th diag update:" << std::endl;
//			std::cout << LL;
		for(unsigned int k=0;k<j;k++){
//				std::cout << "-"
//						  << MatL[j][k]
//					      << "^2";
			LL -= MatL[j][k] * MatL[j][k];
			if(LL<0){
//					std::cout << "=" LL << std::endl << std::endl;
//					std::cout << "something is wrong at ["
//							  << j << "]["
//							  << k << "]=";
				exit(-1);
			}
		}
		MatL[j][j] = (T) std::sqrt((double)LL);
//			std::cout << "=" << MatL[j][j] << std::endl << std::endl;

		// compute the non-diagonal element
//			std::cout << j << "th non-diag update:" << std::endl;
		T inv = 1 / MatL[j][j];
		for(unsigned int i=j+1;i<M;i++){
			LL = Mat[i][j];
//				std::cout << LL;
			for(unsigned int k=0;k<j;k++){
//					std::cout << "-"
//							  << MatL[i][k]
//						      << "x"
//							  << MatL[j][k];
				LL -= MatL[i][k] * MatL[j][k];
			}
			MatL[i][j] = LL * inv;
//				std::cout << "x" << inv << "=" << MatL[i][j] << std::endl;
		}
//			std::cout << std::endl;

//			std::cout << "L Matrix: " << std::endl;;
//			for(int i=0;i<4;i++){
//				for(int j=0;j<4;j++){
//					std::cout << MatL[i][j] << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
	}

	// transpose L to get U
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatU[i][j] = MatL[j][i];
		}
	}
}

// Doolittle algorithm LU decomposition
template<class T, unsigned int M>
void LU_DOOLITTLE(T Mat[M][M], T MatL[M][M], T MatU[M][M]){

	// clean the matrix
	for(unsigned int i=0;i<M;i++){
		for(unsigned int j=0;j<M;j++){
			MatL[i][j] = 0;
			MatU[i][j] = 0;
		}
	}

	// decomposition
	for(unsigned int i=0;i<M;i++){
		// upper triangular
		for(unsigned int k=0;k<M;k++){
			T tmp = Mat[i][k];
			for(unsigned int j=0;j<i;j++){
				tmp -= MatL[i][j] * MatU[j][k];
			}
			MatU[i][k] = tmp;
		}
		// lower triangular
		MatL[i][i] = 1;
		for(unsigned int k=i+1;k<M;k++){
			T tmp = Mat[k][i];
			for(unsigned int j=0;j<i;j++){
				tmp -= MatL[k][j] * MatU[j][i];
			}
			MatL[k][i] = tmp / MatU[i][i];
		}
	}
}




#endif /* HEADER_DSP_LU_HPP_ */
