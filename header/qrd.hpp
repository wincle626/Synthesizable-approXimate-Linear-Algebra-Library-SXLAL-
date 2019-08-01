/*
 * qrd.hpp
 *
 *  Created on: Aug 1, 2019
 *      Author: yunwu
 */

#ifndef SRC_QRD_HPP_
#define SRC_QRD_HPP_

#include "common.hpp"

class QRD_EIGEN{
public:

	// QR decomposition
	template<class EigenT>
	void QRD( EigenT &Mat,
			  EigenT &Q,
			  EigenT &R,
			  int row, int col){
		Eigen::HouseholderQR<EigenT> qr(Mat);
		Q = qr.householderQ(); // get Q matrix
		if( row == col ){
			R = qr.matrixQR().template  triangularView<Eigen::Upper>();
		}else{
			Eigen::FullPivLU<EigenT>lu_decomp(Mat);
			int Rank = lu_decomp.rank(); //retrieve rank of matrix
			R = qr.matrixQR().topLeftCorner(Rank, Rank).template
				triangularView<Eigen::Upper>(); // get R matrix
		}
	}

};

class QRD_C{
public:
	// QR Decomposition
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
	template<class T, int M, int N>
	void QRD_HH(T Mat[M][N],
			    T MatQ[M][N],
			    T MatR[N][N]){

	}
	template<class T, int M, int N>
	void QRD_GR(T Mat[M][N],
			    T MatQ[M][N],
			    T MatR[N][N]){

	}
};


#endif /* SRC_QRD_HPP_ */
