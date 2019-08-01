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

class FILE_IO{
public:
	template<class T>
	// Eigen::MatrixXd, Eigen::MatrixXf,, Eigen::MatrixXi
	int writematrix_eigen(T inputMatrix,
				 const std::string fileName,
				 const std::streamsize dPrec){

		int i, j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (i = 0; i < inputMatrix.rows(); i++) {
			for (j = 0; j < inputMatrix.cols(); j++) {
				outputData << inputMatrix(i, j);
				if (j < (inputMatrix.cols() - 1))
					outputData << ",";
			}
			if (i < (inputMatrix.rows() - 1))
				outputData << std::endl;
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;

	}

	template<class T, int M, int N>
	// double, float, int, or something else
	int writematrix_c(T inputMatrix[M][N],
				 const std::string fileName,
				 const std::streamsize dPrec){

		int i, j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (i = 0; i < M; i++) {
			for (j = 0; j < N; j++) {
				outputData << inputMatrix[i][j];
				if (j < (N - 1))
					outputData << ",";
			}
			if (i < (M - 1))
				outputData << std::endl;
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;

	}

	template<class T>
	// Eigen::VectorXd, Eigen::VectorXf,, Eigen::VectorXi
	int writevector_eigen(T inputVector,
				 const std::string fileName,
				 const std::streamsize dPrec) {

		int j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (j = 0; j < inputVector.size(); j++) {
			outputData << inputVector(j);
			if (j < (inputVector.size() - 1))
				outputData << ",";
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;

	}

	template<class T, int LEN>
	// double, float, int, or something else
	int writevector_c(T inputVector[LEN],
				 const std::string fileName,
				 const std::streamsize dPrec) {

		int j;
		std::ofstream outputData;
		outputData.open(fileName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (j = 0; j < LEN; j++) {
			outputData << inputVector[j];
			if (j < (LEN - 1))
				outputData << ",";
		}
		outputData.close();
		if (!outputData)
			return -1;
		return 0;

	}

	template<class T>
	// double, float, int, or something else
	int writescalar(T num,
			const std::string filename,
			const std::streamsize dPrec){

		std::ofstream outputData;
		outputData.open(filename);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		outputData << num;
		outputData.close();
		if(!outputData)
			return -1;
		return 0;

	}

	template<class T>
	// Eigen::MatrixXd, Eigen::MatrixXf,, Eigen::MatrixXi
	void readmatrix_eigen(std::string filename,
			T &result)
	{

		int cols = 0, rows = 0;
		//double buff[MAXBUFSIZE];

		// Read numbers from file into reesult.
		std::ifstream infile;
		infile.open(filename);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);
			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				temp_cols += 1;
			}
			if (cols == 0)
				cols = temp_cols;
			rows++;
		}
		result.resize(rows, cols);
		infile.close();
		infile.open(filename);
		rows = 0;
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				//buff[cols*rows+temp_cols++] = std::stod(word);
				result(rows, temp_cols) = std::stod(word);
				temp_cols += 1;
			}

			rows++;
		}
		infile.close();

	}

	template<class T, int M, int N>
	// double, float, int, or something else
	void readmatrix_c(std::string filename,
			T result[M][N])
	{

		int cols = 0, rows = 0;
		//double buff[MAXBUFSIZE];

		// Read numbers from file into reesult.
		std::ifstream infile;
		infile.open(filename);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);
			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				temp_cols += 1;
			}
			if (cols == 0)
				cols = temp_cols;
			rows++;
		}
		infile.close();
		infile.open(filename);
		rows = 0;
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				//buff[cols*rows+temp_cols++] = std::stod(word);
				result[rows][temp_cols] = (T) std::stod(word);
				temp_cols += 1;
			}

			rows++;
		}
		infile.close();

	}

	template<class T, class T1>
	// T: Eigen::VectorXd, Eigen::VectorXf,, Eigen::VectorXi
	// T1: double, float, int, or something else
	void readvector_eigen(std::string filename,
			T &result)
	{

		int size = 0;
		T1 buff[MAXBUFSIZE];

		// Read numbers from file into buffer.
		std::ifstream infile;
		infile.open(filename);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				buff[temp_cols++] = (T1) std::stod(word);
			}

			if (size == 0)
				size = temp_cols;

		}

		// Populate matrix with numbers.
		result.resize(size);
		for (int j = 0; j < size; j++)
			result(j) = buff[ j ];

	}

	template<class T, int LEN>
	// double, float, int, or something else
	void readvector_c(std::string filename,
			T result[LEN])
	{

		int size = 0;
		T buff[MAXBUFSIZE];

		// Read numbers from file into buffer.
		std::ifstream infile;
		infile.open(filename);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_cols = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				buff[temp_cols++] = (T) std::stod(word);
			}

			if (size == 0)
				size = temp_cols;

		}

		// Populate matrix with numbers.
		for (int j = 0; j < size; j++)
			result[j] = buff[ j ];

	}

	template<class T>
	// double, float, int, or something else
	void readscalar(std::string filename,
			T &num){

		std::ifstream infile;
		infile.open(filename);
		std::string line;
		std::getline(infile, line);
		num = (T) std::stod(line);

	}
};

#endif /* SRC_COMMON_HPP_ */

