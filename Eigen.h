//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// Wrap Eigen Class with Typedefs

#pragma once

#ifndef _DEBUG
#undef EIGEN_NO_DEBUG // undef EIGEN_NO_DEBUG to remove bounds checking
#endif

#include <Eigen/Eigen/Dense>
#include <Eigen/unsupported/Eigen/KroneckerProduct>
#include <Eigen/unsupported/Eigen/MatrixFunctions>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// Define Constants to Rename Eigen Types - d for double f for float
typedef double EIGENFLOAT; // Any floating variables that go into a matrix should be declared by this typedef

typedef Eigen::MatrixXd MATRIXX;
typedef Eigen::Matrix2d MATRIX2;
typedef Eigen::Matrix3d MATRIX3;
typedef Eigen::Matrix4d MATRIX4;

typedef Eigen::VectorXd VECTORX;
typedef Eigen::Vector2d VECTOR2;
typedef Eigen::Vector3d VECTOR3;
typedef Eigen::Vector4d VECTOR4;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMATRIXX;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RVECTORX;

template<typename M>
void PrintMatrix(const M & Mat) {
	static const Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	static const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
	
	std::cout << Mat.format(CommaInitFmt) << "\n";
}

template<typename M>
void WriteMatrix(const std::string & CSVName, const M & Mat)
{
	static const Eigen::IOFormat Comma(Eigen::StreamPrecision, Eigen::DontAlignCols, ",");

	std::ofstream MatFile(CSVName, std::ofstream::trunc);
	MatFile << Mat.format(Comma);
}

template<typename M>
M ReadMatrix(const std::string & CSVName)
{
	std::ifstream indata;
	indata.open(CSVName);
	if (!indata.is_open()) throw std::runtime_error("Couldn't read file.");
	std::string line;
	std::vector<double> values;
	uint16_t rows = 0;

	while (std::getline(indata, line))
	{
		std::stringstream lineStream(line);
		std::string cell;
		while (std::getline(lineStream, cell, ','))
		{
			values.push_back(std::stod(cell));
		}
		++rows;
	}

	return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}

// Define blkdiag Matlab Function
// Concatenates Matrix1 and Matrix2 on the diagonal
MATRIXX blkdiag(const MATRIXX & Matrix1, const MATRIXX & Matrix2);

// Concatenates Matrix A and Matrix B on the diagonal and returns in A
void blkdiag_inplace(MATRIXX & A, const MATRIXX & B);

// Define repmat Matlab Function
// Create new MatrixB as DesNumRows x DesNumCols of M
// For higher dimension matrices, See Eigen::Tensor broadcast()
MATRIXX repmat(MATRIXX M, int RowReps, int ColReps);
MATRIXX repmat(VECTORX M, int RowReps, int ColReps);

// Define reshape Matlab function
// Reshape 2D M into a new 2D MatrixB of DesNumRows x DesNumCols
// For higher dimension matrices, See Eigen::Tensor reshape()
MATRIXX reshape(MATRIXX M, int DesNumRows, int DesNumCols);

// Define c2d Matlab Function for a zero-order hold
// Discretizes a Continuous State Space Model
// Works with Non-Invertible Matrices
void ZeroOrderHold(MATRIXX & MatrixA, VECTORX & MatrixB, const double & CtrlTS);
void ZeroOrderHold(MATRIX2 & MatrixA, VECTOR2 & MatrixB, const double & CtrlTS);
void ZeroOrderHold(MATRIX3 & MatrixA, VECTOR3 & MatrixB, const double & CtrlTS);
void ZeroOrderHold(MATRIX4 & MatrixA, VECTOR4 & MatrixB, const double & CtrlTS);