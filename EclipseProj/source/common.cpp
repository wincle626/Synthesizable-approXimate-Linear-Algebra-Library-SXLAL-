/*
 *	Author: Yun Wu
 *	Created by: 2019-06-17
 *	Copyright @ Yun Wu
 *
 */

#include "../extlib/fxp/opencv2/gpu/device/common.hpp"

bool checkfileexist(std::string filename){
	std::fstream fileStream;
	fileStream.open(filename.c_str());
	if (fileStream.is_open()) {
		fileStream.close();
	    return true;
	}
	fileStream.close();
	return false;
}

int writematrix_double(const Eigen::MatrixXd& inputMatrix,
			 const std::string& fileName,
			 const std::streamsize dPrec) {

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

int writematrix_float(const Eigen::MatrixXf& inputMatrix,
			 const std::string& fileName,
			 const std::streamsize dPrec) {

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

int writevector_double(const Eigen::VectorXd& inputVector,
			 const std::string& fileName,
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

int writevector_float(const Eigen::VectorXf& inputVector,
			 const std::string& fileName,
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

int writescalar_double(double num, const std::string& filename,
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

int writescalar_float(float num, const std::string& filename,
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

int writescalar_integer(int num, const std::string& filename,
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

double readscalar_double(std::string filename){

	double num;
	std::ifstream infile;
	infile.open(filename);
	std::string line;
	std::getline(infile, line);
	num = std::stod(line);
	return num;

}

float readscalar_float(std::string filename){

	float num;
	std::ifstream infile;
	infile.open(filename);
	std::string line;
	std::getline(infile, line);
	num = std::stof(line);
	return num;

}

int readscalar_integer(std::string filename){

	int num;
	std::ifstream infile;
	infile.open(filename);
	std::string line;
	std::getline(infile, line);
	num = std::stoi(line);
	return num;

}

Eigen::MatrixXd readmatrix_double(std::string filename)
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
	Eigen::MatrixXd result(rows,cols);
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
/*
	// Populate matrix with numbers.
	Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];
*/
    return result;

}

Eigen::MatrixXf readmatrix_float(std::string filename)
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
	Eigen::MatrixXf result(rows,cols);
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
/*
	// Populate matrix with numbers.
	Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];
*/
    return result;

}

//Eigen::MatrixXf readmatrix_float(std::string filename)
//{
//
//	int cols = 0, rows = 0;
//	float buff[MAXBUFSIZE];
//
//	// Read numbers from file into buffer.
//	std::ifstream infile;
//	infile.open(filename);
//	while (! infile.eof())
//	{
//		std::string line;
//		std::getline(infile, line);
//
//		int temp_cols = 0;
//		std::stringstream stream(line);
//		std::string word;
//		while(std::getline(stream, word, ',')){
//			buff[cols*rows+temp_cols++] = std::stof(word);
//		}
//
//		if (cols == 0)
//			cols = temp_cols;
//
//		rows++;
//	}
//
//	// Populate matrix with numbers.
//	Eigen::MatrixXf result(rows,cols);
//    for (int i = 0; i < rows; i++)
//        for (int j = 0; j < cols; j++)
//            result(i,j) = buff[ cols*i+j ];
//
//    return result;
//
//}

Eigen::VectorXd readvector_double(std::string filename)
{

	int size = 0;
//	double buff[MAXBUFSIZE];

	// Read numbers from file into buffer.
	std::ifstream infile;
	infile.open(filename);
	while (! infile.eof())
	{
		std::string line;
		std::getline(infile, line);

//		int temp_cols = 0;
		std::stringstream stream(line);
		std::string word;
		while(std::getline(stream, word, ',')){
//			buff[temp_cols++] = std::stod(word);
			size += 1;
		}

//		if (size == 0)
//			size = temp_cols;

	}
	infile.close();
	infile.open(filename);
	Eigen::VectorXd result(size);
	while (! infile.eof()){
		std::string line;
		std::getline(infile, line);

		int temp_cols = 0;
		std::stringstream stream(line);
		std::string word;
		while(std::getline(stream, word, ',')){
			result(temp_cols++) = std::stod(word);
		}
	}

//	// Populate matrix with numbers.
//	Eigen::VectorXd result(size);
//	for (int j = 0; j < size; j++)
//		result(j) = buff[ j ];
//
    return result;

}

Eigen::VectorXf readvector_float(std::string filename)
{

	int size = 0;
//	double buff[MAXBUFSIZE];

	// Read numbers from file into buffer.
	std::ifstream infile;
	infile.open(filename);
	while (! infile.eof())
	{
		std::string line;
		std::getline(infile, line);

//		int temp_cols = 0;
		std::stringstream stream(line);
		std::string word;
		while(std::getline(stream, word, ',')){
//			buff[temp_cols++] = std::stod(word);
			size += 1;
		}

//		if (size == 0)
//			size = temp_cols;

	}
	infile.close();
	infile.open(filename);
	Eigen::VectorXf result(size);
	while (! infile.eof()){
		std::string line;
		std::getline(infile, line);

		int temp_cols = 0;
		std::stringstream stream(line);
		std::string word;
		while(std::getline(stream, word, ',')){
			result(temp_cols++) = std::stod(word);
		}
	}

//	// Populate matrix with numbers.
//	Eigen::VectorXd result(size);
//	for (int j = 0; j < size; j++)
//		result(j) = buff[ j ];
//
    return result;

}

//Eigen::VectorXf readvector_float(std::string filename)
//{
//
//	int size = 0;
//	float buff[MAXBUFSIZE];
//
//	// Read numbers from file into buffer.
//	std::ifstream infile;
//	infile.open(filename);
//	while (! infile.eof())
//	{
//		std::string line;
//		std::getline(infile, line);
//
//		int temp_cols = 0;
//		std::stringstream stream(line);
//		std::string word;
//		while(std::getline(stream, word, ',')){
//			buff[temp_cols++] = std::stod(word);
//		}
//
//		if (size == 0)
//			size = temp_cols;
//
//	}
//
//	// Populate matrix with numbers.
//	Eigen::VectorXf result(size);
//	for (int j = 0; j < size; j++)
//		result(j) = buff[ j ];
//
//    return result;
//
//}
