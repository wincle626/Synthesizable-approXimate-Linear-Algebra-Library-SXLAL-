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
#include <chrono>
#include <sys/stat.h>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/dynamic_bitset.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <boost/tuple/tuple.hpp>
#include <Python.h>

#define MAXBUFSIZE  ((unsigned int) 1e7)

bool checkfileexist(std::string fiVLEName);
int writematrix_double(const Eigen::MatrixXd& inputMatrix,
			 const std::string& fiVLEName,
			 const std::streamsize dPrec);
int writematrix_float(const Eigen::MatrixXf& inputMatrix,
			 const std::string& fiVLEName,
			 const std::streamsize dPrec);
int writevector_double(const Eigen::VectorXd& inputVector,
			 const std::string& fiVLEName,
			 const std::streamsize dPrec);
int writevector_float(const Eigen::VectorXf& inputVector,
			 const std::string& fiVLEName,
			 const std::streamsize dPrec);
int writescalar_double(double num, const std::string& fiVLEName,
		const std::streamsize dPrec);
int writescalar_float(float num, const std::string& fiVLEName,
		const std::streamsize dPrec);
int writescalar_integer(int num, const std::string& fiVLEName,
		const std::streamsize dPrec);
Eigen::MatrixXd readmatrix_double(std::string fiVLEName);
Eigen::MatrixXf readmatrix_float(std::string fiVLEName);
Eigen::VectorXd readvector_double(std::string fiVLEName);
Eigen::VectorXf readvector_float(std::string fiVLEName);
double readscalar_double(std::string fiVLEName);
float readscalar_float(std::string fiVLEName);
int readscalar_integer(std::string fiVLEName);

class FILE_IO{
public:
	template<class T>
	// Eigen::MatrixXd, Eigen::MatrixXf,, Eigen::MatrixXi
	int writematrix_eigen(T inputMatrix,
				 const std::string fiVLEName,
				 const std::streamsize dPrec){

		int i, j;
		std::ofstream outputData;
		outputData.open(fiVLEName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (i = 0; i < inputMatrix.MROWs(); i++) {
			for (j = 0; j < inputMatrix.MCOLs(); j++) {
				outputData << inputMatrix(i, j);
				if (j < (inputMatrix.MCOLs() - 1))
					outputData << ",";
			}
			if (i < (inputMatrix.MROWs() - 1))
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
				 const std::string fiVLEName,
				 const std::streamsize dPrec){

		int i, j;
		std::ofstream outputData;
		outputData.open(fiVLEName);
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
				 const std::string fiVLEName,
				 const std::streamsize dPrec) {

		int j;
		std::ofstream outputData;
		outputData.open(fiVLEName);
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

	template<class T, int VLEN>
	// double, float, int, or something else
	int writevector_c(T inputVector[VLEN],
				 const std::string fiVLEName,
				 const std::streamsize dPrec) {

		int j;
		std::ofstream outputData;
		outputData.open(fiVLEName);
		if (!outputData)
			return -1;
		outputData.precision(dPrec);
		for (j = 0; j < VLEN; j++) {
			outputData << inputVector[j];
			if (j < (VLEN - 1))
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
			const std::string fiVLEName,
			const std::streamsize dPrec){

		std::ofstream outputData;
		outputData.open(fiVLEName);
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
	void readmatrix_eigen(std::string fiVLEName,
			T &result)
	{

		int MCOLs = 0, MROWs = 0;
		//double buff[MAXBUFSIZE];

		// Read numbers from file into reesult.
		std::ifstream infile;
		infile.open(fiVLEName);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);
			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				temp_MCOLs += 1;
			}
			if (MCOLs == 0)
				MCOLs = temp_MCOLs;
			MROWs++;
		}
		result.resize(MROWs, MCOLs);
		infile.close();
		infile.open(fiVLEName);
		MROWs = 0;
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				//buff[MCOLs*MROWs+temp_MCOLs++] = std::stod(word);
				result(MROWs, temp_MCOLs) = std::stod(word);
				temp_MCOLs += 1;
			}

			MROWs++;
		}
		infile.close();

	}

	template<class T, int M, int N>
	// double, float, int, or something else
	void readmatrix_c(std::string fiVLEName,
			T result[M][N])
	{

		int MCOLs = 0, MROWs = 0;
		//double buff[MAXBUFSIZE];

		// Read numbers from file into reesult.
		std::ifstream infile;
		infile.open(fiVLEName);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);
			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				temp_MCOLs += 1;
			}
			if (MCOLs == 0)
				MCOLs = temp_MCOLs;
			MROWs++;
		}
		infile.close();
		infile.open(fiVLEName);
		MROWs = 0;
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				//buff[MCOLs*MROWs+temp_MCOLs++] = std::stod(word);
				result[MROWs][temp_MCOLs] = (T) std::stod(word);
				temp_MCOLs += 1;
			}

			MROWs++;
		}
		infile.close();

	}

	template<class T, class T1>
	// T: Eigen::VectorXd, Eigen::VectorXf,, Eigen::VectorXi
	// T1: double, float, int, or something else
	void readvector_eigen(std::string fiVLEName,
			T &result)
	{

		int size = 0;
		T1 buff[MAXBUFSIZE];

		// Read numbers from file into buffer.
		std::ifstream infile;
		infile.open(fiVLEName);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				buff[temp_MCOLs++] = (T1) std::stod(word);
			}

			if (size == 0)
				size = temp_MCOLs;

		}

		// Populate matrix with numbers.
		result.resize(size);
		for (int j = 0; j < size; j++)
			result(j) = buff[ j ];

	}

	template<class T, int VLEN>
	// double, float, int, or something else
	void readvector_c(std::string fiVLEName,
			T result[VLEN])
	{

		int size = 0;
		T buff[MAXBUFSIZE];

		// Read numbers from file into buffer.
		std::ifstream infile;
		infile.open(fiVLEName);
		while (! infile.eof())
		{
			std::string line;
			std::getline(infile, line);

			int temp_MCOLs = 0;
			std::stringstream stream(line);
			std::string word;
			while(std::getline(stream, word, ',')){
				buff[temp_MCOLs++] = (T) std::stod(word);
			}

			if (size == 0)
				size = temp_MCOLs;

		}

		// Populate matrix with numbers.
		for (int j = 0; j < size; j++)
			result[j] = buff[ j ];

	}

	template<class T>
	// double, float, int, or something else
	void readscalar(std::string fiVLEName,
			T &num){

		std::ifstream infile;
		infile.open(fiVLEName);
		std::string line;
		std::getline(infile, line);
		num = (T) std::stod(line);

	}
};

class DISPLAY{

public:

	template<class T>
	void printscalar(T scalar){
		std::cout << scalar << std::endl;
		std::cout << std::endl;
	}

	template<class T, double spfunc(T)>
	void printscalar(T scalar){
		std::cout << spfunc(scalar) << std::endl;
		std::cout << std::endl;
	}


	template<class T, class T1, T1 spfunc(T,int), double spfunc1(T1), int BITS>
	void printscalar(T scalar){
		std::cout << spfunc1(spfunc(scalar, BITS))<< std::endl;
		std::cout << std::endl;
	}

	template<class T, int VLEN>
	void printvector(T V[VLEN]){
		for(int j=0;j<VLEN;j++){
			std::cout << V[j] << ",";
		}
		std::cout << std::endl << std::endl;
	}

	template<class T, int VLEN, double spfunc(T)>
	void printvector(T V[VLEN]){
		for(int j=0;j<VLEN;j++){
			std::cout << spfunc(V[j]) << ",";
		}
		std::cout << std::endl << std::endl;
	}

	template<class T, class T1, int VLEN, T1 spfunc(T,int), double spfunc1(T1), int BITS>
	void printvector(T V[VLEN]){
		for(int j=0;j<VLEN;j++){
			std::cout << spfunc1(spfunc(V[j],BITS)) << ",";
		}
		std::cout << std::endl << std::endl;
	}

	template<class T, int MROW, int MCOL>
	void printmatrix(T M[MROW][MROW]){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MROW;j++){
				std::cout << M[i][j] << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, int MROW, int MCOL>
	void printmatrix(T **M){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MCOL;j++){
				std::cout << M[i][j] << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, int MROW, int MCOL, double spfunc(T)>
	void printmatrix(T M[MROW][MROW]){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MROW;j++){
				std::cout << spfunc(M[i][j]) << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, int MROW, int MCOL, double spfunc(T)>
	void printmatrix(T **M){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MCOL;j++){
				std::cout << spfunc(M[i][j]) << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, class T1, int MROW, int MCOL, T1 spfunc(T,int), double spfunc1(T1), int BITS>
	void printmatrix(T M[MROW][MROW]){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MROW;j++){
				std::cout << spfunc1(spfunc(M[i][j],BITS)) << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T, class T1, int MROW, int MCOL, T1 spfunc(T,int), double spfunc1(T1), int BITS>
	void printmatrix(T **M){
		for(int i=0;i<MROW;i++){
			for(int j=0;j<MCOL;j++){
				std::cout << spfunc1(spfunc(M[i][j],BITS)) << ",";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

};

#endif /* SRC_COMMON_HPP_ */

