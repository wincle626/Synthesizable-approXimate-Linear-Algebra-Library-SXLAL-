/*
 * proximal_gradient_decent_eigenxd.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: yunwu
 */

#include "pgd_eigend.hpp"

namespace plt = matplotlibcpp;

void PROXIMAL_GRADIENT_DECENT_EIGEND(Eigen::MatrixXd Amatrix,
									 Eigen::VectorXd bvector,
									 double L){
#ifdef EIGEN_DOUBLE_PRECISION
	Eigen_Algebra Eigen_Algebra_obj;
	double error = 0;
	double error_std = ERR_STD;

	std::string clockname = "timeprofile_eigend.txt";
	std::string xkname = "xk_eigend.dat";
	std::string errorrecordname = "error_record_eigend.dat";
	std::string errorhistname = "error_hist_eigend.dat";
	std::string figurename = "ProximalGradientDecent_eigend.png";

#ifdef TIME_PROFILE
	clock_t start = clock();
	std::ofstream TimeProfile;
	TimeProfile.open(clockname);
#endif
	//x_k = zeros(n, 1);
	Eigen::VectorXd x_k_vec( DIAG );
	Eigen::VectorXd x_k_plus1_vec( DIAG );
	Eigen::VectorXd grad_g( DIAG );
	Eigen::VectorXd ones_vec( DIAG );
	Eigen::VectorXd error_hist ( ITER_MAX );
	std::vector<double> error_record(ITER_MAX, 0);
	std::vector<std::vector<double>> x_k_record;

	Eigen_Algebra_obj.ONES_VEC<Eigen::VectorXd, DIAG>( ones_vec );
	Eigen_Algebra_obj.ZEROS_VEC<Eigen::VectorXd, DIAG>( x_k_vec );
	Eigen_Algebra_obj.ZEROS_VEC<Eigen::VectorXd, DIAG>( error_hist );

	int k=0;
	while( k<ITER_MAX ){

		// Save previous point
		x_k_plus1_vec = x_k_vec;

		// Compute gradient
		Eigen::VectorXd tmp_mul_vec1( DIAG );
		Eigen_Algebra_obj.MAT_VEC_MUL<Eigen::MatrixXd, Eigen::VectorXd>( Amatrix, x_k_vec,
									tmp_mul_vec1 );
#if defined(DEBUG_LAFUNC)
		std::cout << "tmp_mul_vec1: " << std::endl
				  << tmp_mul_vec1
				  << std::endl << std::endl;
#endif
		Eigen_Algebra_obj.VEC_ADD<Eigen::VectorXd>( tmp_mul_vec1, bvector,
								grad_g );
#if defined(DEBUG_LAFUNC)
		std::cout << "grad_g: " << std::endl
				  << grad_g
				  << std::endl << std::endl;
#endif
#ifdef DEBUG_ITER
		std::cout << grad_g << std::endl << std::endl;
#endif

		// new decent point
		Eigen::VectorXd new_point( DIAG );
		Eigen::VectorXd tmp_mul_vec2( DIAG );
		Eigen_Algebra_obj.VEC_SCALAR_MUL<Eigen::VectorXd, double>( grad_g, 1/L,
									 tmp_mul_vec2 );
#if defined(DEBUG_LAFUNC)
		std::cout << "tmp_mul_vec2: " << std::endl
				  << tmp_mul_vec2
				  << std::endl << std::endl;
#endif
		Eigen_Algebra_obj.VEC_SUB<Eigen::VectorXd>( x_k_vec, tmp_mul_vec2,
								new_point );
#if defined(DEBUG_LAFUNC)
		std::cout << "new_point: " << std::endl
				  << new_point
				  << std::endl << std::endl;
#endif
#ifdef DEBUG_ITER
		std::cout << new_point << std::endl << std::endl;
#endif

		// Proximal projection
		Eigen::VectorXd tmp_min( DIAG );
		Eigen::VectorXd tmp_max( DIAG );
		Eigen::VectorXd rndnoise( DIAG );
		Eigen::VectorXd tmp_mul_vec3( DIAG );
		Eigen_Algebra_obj.RND_VEC<Eigen::VectorXd, DIAG>( rndnoise );
		Eigen_Algebra_obj.VEC_SCALAR_MUL<Eigen::VectorXd, double>( rndnoise, error_std,
									 tmp_mul_vec3);
		Eigen_Algebra_obj.VEC_SCALAR_MIN<Eigen::VectorXd, double, Eigen::ArrayXd, DIAG>( new_point,
				BOX_CONST, tmp_min );
		Eigen_Algebra_obj.VEC_SCALAR_MAX<Eigen::VectorXd, double, Eigen::ArrayXd, DIAG>( tmp_min,
				-BOX_CONST, tmp_max );
		Eigen_Algebra_obj.VEC_ADD<Eigen::VectorXd>( tmp_max.matrix(), tmp_mul_vec3,
								x_k_vec);
#ifdef DEBUG_ITER
		std::cout << x_k_vec << std::endl << std::endl;
#endif

		// check early termination constraint
		double norm1, norm2;
		Eigen::VectorXd tmp_sub_vec( DIAG );
		Eigen_Algebra_obj.VEC_SUB<Eigen::VectorXd>( x_k_vec, x_k_plus1_vec,
								tmp_sub_vec);
		norm1 = tmp_sub_vec.norm();
		norm2 = x_k_vec.norm();
		if( norm1<= EPS_STOP ){
			break;
		}

		// record relative error
		error_hist(k) = norm1 / norm2;
		error = norm1;
#ifdef RECORD_RESULT
		error_hist(k) = norm1 / norm2;
		error_record.at(k) = error;
		std::vector<double> xkvec;
		for(int i=0;i<DIAG;i++)
			xkvec.push_back((double)x_k_vec[i]);
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
			std::vector<double> xkvec = x_k_record.at(i);
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



