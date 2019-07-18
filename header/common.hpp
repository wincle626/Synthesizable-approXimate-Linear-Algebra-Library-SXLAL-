/*
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_COMMON_HPP_
#define SRC_COMMON_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <bitset>
#include <climits>
#include <cstdint>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/dynamic_bitset.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <boost/tuple/tuple.hpp>

#define MAXBUFSIZE  ((int) 1e6)

bool checkfileexist(std::string filename);
int writematrix_double(const Eigen::MatrixXd& inputMatrix,
			 const std::string& fileName,
			 const std::streamsize dPrec);
int writematrix_float(const Eigen::MatrixXd& inputMatrix,
			 const std::string& fileName,
			 const std::streamsize dPrec);
int writevector_double(const Eigen::VectorXd& inputVector,
			 const std::string& fileName,
			 const std::streamsize dPrec);
int writevector_float(const Eigen::VectorXf& inputVector,
			 const std::string& fileName,
			 const std::streamsize dPrec);
int writescalar_double(double num, const std::string& filename,
		const std::streamsize dPrec);
int writescalar_float(float num, const std::string& filename,
		const std::streamsize dPrec);
int writescalar_integer(int num, const std::string& filename,
		const std::streamsize dPrec);
Eigen::MatrixXd readmatrix_double(std::string filename);
Eigen::MatrixXf readmatrix_float(std::string filename);
Eigen::VectorXd readvector_double(std::string filename);
Eigen::VectorXf readvector_float(std::string filename);
double readscalar_double(std::string filename);
float readscalar_float(std::string filename);
int readscalar_integer(std::string filename);

#endif /* SRC_COMMON_HPP_ */

