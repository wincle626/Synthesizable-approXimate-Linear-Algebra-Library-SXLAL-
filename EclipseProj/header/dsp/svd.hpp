/*
 * svd.hpp
 *
 *  Created on: Mar 25, 2023
 *      Author: ubuntu
 */

#ifndef SRC_SVD_HPP_
#define SRC_SVD_HPP_

#include "qrd.hpp"

template<class T, int M, int N>
void svd_jacobi(T AM[M][N], T UM[M][M],
           T SM[M][N], T SMdiag[M],
           T VM[N][N]){

    T x[M], y[M], tmp1[M], tmp2[M], AMtmp[M][N];
    T a, b, c, e, e_sign, e_abs, t, cs, sn;

    EYE_MATRIX<double, N, N>(VM, N);
    EYE_MATRIX<double, M, N>(SM, N);
    MAT_EQ<double, M, N>(AM, AMtmp);
    for(int iter=0;iter<30;iter++){
        for(int i=0;i<N-1;i++){
            for(int j=i+1;j<N;j++){
                // Mumnwise
                COL_OF_MATRIX<double, M, N>(AMtmp, i, x);
//                pirntvector<double, M>(x);
                COL_OF_MATRIX<double, M, N>(AMtmp, j, y);
//                pirntvector<double, M>(y);
                VEC_DOT_PROD<double, M>(x, x, a);
                VEC_DOT_PROD<double, M>(y, y, b);
                VEC_DOT_PROD<double, M>(x, y, c);
                // Jacobi rotation
                e = (b - a) / (2 * c);
                e_sign = e > 0 ? 1 : -1;
                e_abs = e > 0 ? e : -e;
                t = e_sign / (e_abs + sqrt(1 + e * e));
                cs = 1 / sqrt(1 + t * t);
                sn = cs * t;
                // Update G
                VEC_SCALAR_MUL<double, M>(x, cs, tmp1);
//                pirntvector<double, M>(tmp1);
                VEC_SCALAR_MUL<double, M>(y, sn, tmp2);
//                pirntvector<double, M>(tmp2);
                VEC_SUB_R<double, M>(tmp1, tmp2, tmp1);
//                pirntvector<double, M>(tmp1);
                COL_INTO_MATRIX<double, M, N>(tmp1, AMtmp, i);
//            	pirntmatrix<double, ROW, COL>(AMtmp);
                VEC_SCALAR_MUL<double, M>(x, sn, tmp1);
//                pirntvector<double, M>(tmp1);
                VEC_SCALAR_MUL<double, M>(y, cs, tmp2);
//                pirntvector<double, M>(tmp2);
                VEC_ADD_R<double, M>(tmp1, tmp2, tmp1);
//                pirntvector<double, M>(tmp1);
                COL_INTO_MATRIX<double, M, N>(tmp1, AMtmp, j);
//            	pirntmatrix<double, ROW, COL>(AMtmp);
                // Update V
                COL_OF_MATRIX<double, M, N>(VM, i, x);
//                pirntvector<double, N>(x);
                COL_OF_MATRIX<double, M, N>(VM, j, y);
//                pirntvector<double, N>(y);
                VEC_SCALAR_MUL<double, N>(x, cs, tmp1);
//                pirntvector<double, N>(tmp1);
                VEC_SCALAR_MUL<double, N>(y, sn, tmp2);
//                pirntvector<double, N>(tmp2);
                VEC_SUB_R<double, N>(tmp1, tmp2, tmp1);
//                pirntvector<double, N>(tmp1);
                COL_INTO_MATRIX<double, N, N>(tmp1, VM, i);
//            	pirntmatrix<double, N, N>(VM);
                VEC_SCALAR_MUL<double, N>(x, sn, tmp1);
//                pirntvector<double, N>(tmp1);
                VEC_SCALAR_MUL<double, N>(y, cs, tmp2);
//                pirntvector<double, N>(tmp2);
                VEC_ADD_R<double, N>(tmp1, tmp2, tmp1);
//                pirntvector<double, N>(tmp1);
                COL_INTO_MATRIX<double, N, N>(tmp1, VM, j);
//            	pirntmatrix<double, N, N>(VM);
//            	printf("Press Any Key to Continue\n");
//            	getchar();
            }
        }
    }
//    printmatrix<double, M, N>(AMtmp);
//    printmatrix<double, N, N>(VM);
    MAT_EQ_S2L<double, M, N, M, M>(AMtmp, UM);
//    printmatrix<double, M, M>(UM);
    for(int i=0;i<N;i++){
        double norm = 0.0;
        COL_OF_MATRIX<double, M, N>(AMtmp, i, tmp1);
//		pirntvector<double, M>(tmp1);
        VEC_L2NORM_R<double, M>(tmp1, norm);
        SCALAR_INTO_MATRIX<double, M, N>(norm, SM, i, i);
//		pirntmatrix<double, M, N>(SM);
        SMdiag[i] = norm;
        VEC_SCALAR_DIV_R<double, M>(tmp1, norm, tmp2);
//		pirntvector<double, M>(tmp2);
        COL_INTO_MATRIX<double, M, M>(tmp2, UM, i);
//    	pirntmatrix<double, M, M>(UM);
    }
//	pirntmatrix<double, M, M>(UM);
//	pirntmatrix<double, M, N>(SM);
}


template<class T, int M, int N>
void svd_qrgr(T AM[M][N], T UM[M][M],
           T SM[M][N], T SMdiag[M],
           T VM[N][N]){

	T SM1[N][M], SM1T[M][N], SM2[M][N];
	T QM[M][M], RM[M][N], UM1[M][M], SM2T[N][M];
	T QM1[N][N], RM1[N][M], VMtmp[N][N];

	// u=eye(sizea(1))
	EYE_MATRIX<T, M, M>(UM, M);
//	printmatrix<T, M, M>(UM);
    // s=a';
	MAT_TRANS<T, M, N>(AM, SM1);
//	printmatrix<T, N, M>(SM1);
    // v=eye(sizea(2))
	EYE_MATRIX<T, N, N>(VM, N);
//	printmatrix<T, N, N>(VM);
    for(int i=0;i<50;i++){
    	// [q,s]=qrdhhr(s');
    	MAT_TRANS<T, N, M>(SM1, SM1T);
		QRD_GR<T, M, N>(SM1T, QM, RM);
    	MAT_EQ<T, M, N>(RM, SM2);
//    	printmatrix<T, M, N>(SM1T);
//    	printmatrix<T, M, M>(QM);
//    	printmatrix<T, M, N>(RM);
    	//  u=u*q;
		MAT_MUL<T, M, M, M>(UM, QM, UM1);
		MAT_EQ<T, M, M>(UM1, UM);
//    	printmatrix<T, M, M>(UM);

    	// [q,s]=qrdhhr(s');
		MAT_TRANS<T, M, N>(SM2, SM2T);
		QRD_GR<T, N, M>(SM2T, QM1, RM1);
		MAT_EQ<T, N, M>(RM1, SM1);
//    	printmatrix<T, N, M>(SM2T);
//    	printmatrix<T, N, N>(QM1);
//    	printmatrix<T, N, M>(RM1);
    	// v=v*q;
		MAT_MUL<T, N, N, N>(VM, QM1, VMtmp);
		MAT_EQ<T, N, N>(VMtmp, VM);
//    	printmatrix<T, N, N>(VM);
    }
    for(int i=0;i<N;i++){
    	T diag = SM1[i][i];
    	SM[i][i] = diag > 0 ? diag : -diag;
    	if(diag<0){
    		for(int j=0;j<M;j++){
    			UM[j][i] = -UM[j][i];
    		}
    	}
    }
}


template<class T, int M, int N>
void svd_qrgs(T AM[M][N], T UM[M][M],
           T SM[M][N], T SMdiag[M],
           T VM[N][N]){

	T SM1[N][M], SM1T[M][N], SM2[M][N];
	T QM[M][N], RM[M][N], UM1[M][N], SM2T[N][M];
	T QM1[N][M], RM1[N][M], VMtmp[N][M];

	// u=eye(sizea(1))
	EYE_MATRIX<T, M, M>(UM, M);
    // s=a';
	MAT_TRANS<T, M, N>(AM, SM1);
    // v=eye(sizea(2))
	EYE_MATRIX<T, N, N>(VM, N);
    for(int i=0;i<50;i++){
    	// [q,s]=qrdhhr(s');
    	MAT_TRANS<T, N, M>(SM1, SM1T);
    	QRD_GS<T, M, N>(SM1T, QM, RM);
    	MAT_EQ_S2L<T, N, N, M, N>(RM, SM2);
    	//  u=u*q;
		MAT_MUL<T, M, M, N>(UM, QM, UM1);
		MAT_EQ_S2L<T, M, N, M, M>(UM1, UM);

    	// [q,s]=qrdhhr(s');
		MAT_TRANS<T, M, N>(SM2, SM2T);
		QRD_GS<T, N, M>(SM2T, QM1, RM1);
		MAT_EQ<T, N, M>(RM1, SM1);
    	// v=v*q;
		MAT_MUL<T, N, N, M>(VM, QM1, VMtmp);
		MAT_EQ_L2S<T, N, M, N, N>(VMtmp, VM);
    }
    for(int i=0;i<N;i++){
    	T diag = SM1[i][i];
    	SM[i][i] = diag > 0 ? diag : -diag;
    	if(diag<0){
    		for(int j=0;j<M;j++){
    			UM[j][i] = -UM[j][i];
    		}
    	}
    }
}


template<class T, int M, int N>
void svd_qrmgs(T AM[M][N], T UM[M][M],
           T SM[M][N], T SMdiag[M],
           T VM[N][N]){

	T SM1[N][M], SM1T[M][N], SM2[M][N];
	T QM[M][N], RM[M][N], UM1[M][N], SM2T[N][M];
	T QM1[N][M], RM1[N][M], VMtmp[N][M];

	// u=eye(sizea(1))
	EYE_MATRIX<T, M, M>(UM, M);
    // s=a';
	MAT_TRANS<T, M, N>(AM, SM1);
    // v=eye(sizea(2))
	EYE_MATRIX<T, N, N>(VM, N);
    for(int i=0;i<50;i++){
    	// [q,s]=qrdhhr(s');
    	MAT_TRANS<T, N, M>(SM1, SM1T);
    	QRD_MGS<T, M, N>(SM1T, QM, RM);
    	MAT_EQ_S2L<T, N, N, M, N>(RM, SM2);
    	//  u=u*q;
		MAT_MUL<T, M, M, N>(UM, QM, UM1);
		MAT_EQ_S2L<T, M, N, M, M>(UM1, UM);

    	// [q,s]=qrdhhr(s');
		MAT_TRANS<T, M, N>(SM2, SM2T);
		QRD_MGS<T, N, M>(SM2T, QM1, RM1);
		MAT_EQ<T, N, M>(RM1, SM1);
    	// v=v*q;
		MAT_MUL<T, N, N, M>(VM, QM1, VMtmp);
		MAT_EQ_L2S<T, N, M, N, N>(VMtmp, VM);
    }
    for(int i=0;i<N;i++){
    	T diag = SM1[i][i];
    	SM[i][i] = diag > 0 ? diag : -diag;
    	if(diag<0){
    		for(int j=0;j<M;j++){
    			UM[j][i] = -UM[j][i];
    		}
    	}
    }
}


template<class T, int M, int N>
void svd_qrhh(T AM[M][N], T UM[M][M],
           T SM[M][N], T SMdiag[M],
           T VM[N][N]){

	T SM1[N][M], SM1T[M][N], SM2[M][N];
	T QM[M][M], RM[M][N], UM1[M][M], SM2T[N][M];
	T QM1[N][N], RM1[N][M], VMtmp[N][N];

	// u=eye(sizea(1))
	EYE_MATRIX<T, M, M>(UM, M);
//	printmatrix<T, M, M>(UM);
    // s=a';
	MAT_TRANS<T, M, N>(AM, SM1);
//	printmatrix<T, N, M>(SM1);
    // v=eye(sizea(2))
	EYE_MATRIX<T, N, N>(VM, N);
//	printmatrix<T, N, N>(VM);
    for(int i=0;i<50;i++){
    	// [q,s]=qrdhhr(s');
    	MAT_TRANS<T, N, M>(SM1, SM1T);
		QRD_HH<T, M, N>(SM1T, QM, RM);
    	MAT_EQ<T, M, N>(RM, SM2);
    	printmatrix<T, M, N>(SM1T);
    	printmatrix<T, M, M>(QM);
    	printmatrix<T, M, N>(RM);
    	//  u=u*q;
		MAT_MUL<T, M, M, M>(UM, QM, UM1);
		MAT_EQ<T, M, M>(UM1, UM);
//    	printmatrix<T, M, M>(UM);

    	// [q,s]=qrdhhr(s');
		MAT_TRANS<T, M, N>(SM2, SM2T);
		QRD_HH<T, N, M>(SM2T, QM1, RM1);
		MAT_EQ<T, N, M>(RM1, SM1);
    	printmatrix<T, N, M>(SM2T);
    	printmatrix<T, N, N>(QM1);
    	printmatrix<T, N, M>(RM1);
    	// v=v*q;
		MAT_MUL<T, N, N, N>(VM, QM1, VMtmp);
		MAT_EQ<T, N, N>(VMtmp, VM);
//    	printmatrix<T, N, N>(VM);
    }
    for(int i=0;i<N;i++){
    	T diag = SM1[i][i];
    	SM[i][i] = diag > 0 ? diag : -diag;
    	if(diag<0){
    		for(int j=0;j<M;j++){
    			UM[j][i] = -UM[j][i];
    		}
    	}
    }
}

#endif /* SRC_SVD_HPP_ */
