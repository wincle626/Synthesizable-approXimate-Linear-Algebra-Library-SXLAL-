/*
 * proximal_gradient_decent_eigenxf.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: yunwu
 */

#include "pgd_test.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_EIGENF(Eigen::MatrixXd Amatrix,
									 Eigen::VectorXd bvector,
									 double L){

#ifdef EIGEN_FLOAT_PRECISION
	Eigen_Algebra Eigen_Algebra_obj;
	float error = 0;
	float error_std = ERR_STD;

	std::string clockname = "timeprofile.txt";
	std::string xkname = "xk.dat";
	std::string errorrecordname = "error_record.dat";
	std::string errorhistname = "error_hist.dat";
	std::string figurename = "ProximalGradientDecent.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif
	//x_k = zeros(n, 1);
	Eigen::VectorXf x_k_vec( DIAG );
	Eigen::VectorXf x_k_plus1_vec( DIAG );
	Eigen::VectorXf grad_g( DIAG );
	Eigen::VectorXf ones_vec( DIAG );
	Eigen::VectorXf error_hist ( ITER_MAX );
   	std::vector<float> error_record(ITER_MAX, 0);
   	std::vector<std::vector<float>> x_k_record;

	Eigen_Algebra_obj.ONES_VEC<Eigen::VectorXf, DIAG>( ones_vec );
	Eigen_Algebra_obj.ZEROS_VEC<Eigen::VectorXf, DIAG>( x_k_vec );
	Eigen_Algebra_obj.ZEROS_VEC<Eigen::VectorXf, DIAG>( error_hist );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		x_k_plus1_vec = x_k_vec;

		// Compute gradient
		Eigen::VectorXf tmp_mul_vec1( DIAG );
		Eigen_Algebra_obj.MAT_VEC_MUL<Eigen::MatrixXf, Eigen::VectorXf>( Amatrix.cast<float>(), x_k_vec,
									tmp_mul_vec1 );
		Eigen_Algebra_obj.VEC_ADD<Eigen::VectorXf>( tmp_mul_vec1, bvector.cast<float>(),
								grad_g );
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		Eigen::VectorXf new_point( DIAG );
		Eigen::VectorXf tmp_mul_vec2( DIAG );
		Eigen_Algebra_obj.VEC_SCALAR_MUL<Eigen::VectorXf, double>( grad_g, 1/(float)L,
									 tmp_mul_vec2 );
		Eigen_Algebra_obj.VEC_SUB<Eigen::VectorXf>( x_k_vec, tmp_mul_vec2,
								new_point );
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		Eigen::VectorXf tmp_min( DIAG );
		Eigen::VectorXf tmp_max( DIAG );
		Eigen::VectorXf rndnoise( DIAG );
		Eigen::VectorXf tmp_mul_vec3( DIAG );
		Eigen_Algebra_obj.RND_VEC<Eigen::VectorXf, DIAG>( rndnoise );
		Eigen_Algebra_obj.VEC_SCALAR_MUL<Eigen::VectorXf, float>( rndnoise, error_std,
									 tmp_mul_vec3);
		Eigen_Algebra_obj.VEC_SCALAR_MIN<Eigen::VectorXf, float, Eigen::ArrayXf, DIAG>( new_point,
				BOX_CONST, tmp_min );
		Eigen_Algebra_obj.VEC_SCALAR_MAX<Eigen::VectorXf, float, Eigen::ArrayXf, DIAG>( tmp_min,
				-BOX_CONST, tmp_max );
		Eigen_Algebra_obj.VEC_ADD<Eigen::VectorXf>( tmp_max.matrix(), tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		double norm1, norm2;
		Eigen::VectorXf tmp_sub_vec( DIAG );
		Eigen_Algebra_obj.VEC_SUB<Eigen::VectorXf>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		norm1 = tmp_sub_vec.norm();
		norm2 = x_k_vec.norm();
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error = norm1;
#ifdef RECORD_RESULT
		error_hist(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<float> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((float)x_k_vec[i]);
		x_k_record.push_back(xkvec);
#endif// endif RECORD_RESULT
#ifdef DEBUG_ITER
		std::cout << error << std::endl << std::endl;
#endif

		k++;

	}
#ifdef TIME_PROFILE
	clock_t end = clock();
	double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout << "It takes "
			  << time
			  << " ms to gradient decent"
			  << std::endl << std::endl;
	TimeProfile << "It takes "
			  << time
			  << " ms to gradient decent"
			  << std::endl << std::endl;
#endif
	///////////////////////////////////////////////////////////



#ifdef RECORD_RESULT
	/////////////////// Print Iteration Number ////////////////
	std::cout << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist(k-2)
			  << std::endl << std::endl;
	TimeProfile << "Number of iterations = " << k-2
			  << " out of " << ITER_MAX << std::endl
			  << "Error: " << error << std::endl
			  << "Relative Error: " << error_hist(k-2) << std::endl;
	///////////////////////////////////////////////////////////



	//////////////////////// Plot Figures /////////////////////
	// save data
	std::ofstream resultfile;
	std::ofstream resultfile1;
	std::ofstream resultfile2;
	resultfile.open(errorhistname);
	resultfile1.open(errorrecordname);
	resultfile2.open(xkname);
	resultfile << "# iteration relative_error \n";
	resultfile1 << "# iteration error \n";
	resultfile2 << "# optimizing x_k vector \n";
	std::vector<double> x, y;
	for(int i=0; i<ITER_MAX; i++){
		if( error_hist(i)>0 ){
			x.push_back(i);
			y.push_back(error_hist(i));
			resultfile << "  " << i << "," << error_hist(i) << " \n";
			resultfile1 << "  " << i << "," << error_record.at(i) << " \n";
			std::vector<float> xkvec = x_k_record.at(i);
			for(int j=0;j<DIAG-1;j++)
				resultfile2 << xkvec.at(j) << ",";
			resultfile2 << xkvec.at(DIAG-1) << "\n";
		}
	}
	resultfile.close();
	resultfile1.close();
	resultfile2.close();

#ifdef PLOT_FIGURE
	// Matplotlib plotting
	plt::named_semilogy( "simple example", x, y );
	plt::title("Proximal Gradient Decent");
	plt::xlabel("Number of Iteration");
	plt::ylabel("Relative Error");
	plt::legend();
	plt::save(figurename);
#ifdef SHOW_FIGURE
   	plt::show();
#endif
#endif
#endif
#endif
}



