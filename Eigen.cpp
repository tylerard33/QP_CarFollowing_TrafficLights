//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////

#include "Eigen.h"

MATRIXX blkdiag(const MATRIXX & Matrix1, const MATRIXX & Matrix2)
{
	// Initialize Zero matrix for new size
	MATRIXX MatrixZero = MATRIXX::Zero
	(
		Matrix1.rows() + Matrix2.rows(),
		Matrix1.cols() + Matrix2.cols()
	);

	// Impose Matrix1 and Matrix2 to MatrixZero
	MatrixZero.topLeftCorner(Matrix1.rows(), Matrix1.cols()) = Matrix1;
	MatrixZero.bottomRightCorner(Matrix2.rows(), Matrix2.cols()) = Matrix2;

	return MatrixZero;
}

void blkdiag_inplace(MATRIXX & A, const MATRIXX & B)
{
	// Resize Mat1
	const int & rows1 = A.rows();
	const int & cols1 = A.cols();
	const int & rows2 = B.rows();
	const int & cols2 = B.cols();

	A.conservativeResize(rows1+rows2, cols1+cols2);

	// Lower Right is Mat2. Upper Left is already Mat1
	A.bottomRightCorner(rows2, cols2) = B;
	if (A.bottomLeftCorner(rows2, cols1).array().any()) A.bottomLeftCorner(rows2, cols1) = MATRIXX::Zero(rows2, cols1);
	if (A.topRightCorner(rows1, cols2).array().any()) A.topRightCorner(rows1, cols2) = MATRIXX::Zero(rows1, cols2);
}

MATRIXX repmat(MATRIXX M, int RowReps, int ColReps)
{
	return M.replicate(RowReps, ColReps);
}

MATRIXX repmat(VECTORX M, int RowReps, int ColReps)
{
	return M.replicate(RowReps, ColReps);
}

MATRIXX reshape(MATRIXX M, int DesNumRows, int DesNumCols)
{
	return Eigen::Map<MATRIXX>(M.data(), DesNumRows, DesNumCols);
}

void ZeroOrderHold(MATRIXX & MatrixA, VECTORX & MatrixB, const double & CtrlTS)
{
	// Find Rows and Columns
	short MAR = MatrixA.rows(); short MAC = MatrixA.cols();
	short MBR = MatrixB.rows(); short MBC = MatrixB.cols();
	
	// Form Control Hierarchy, Compute Matrix Exponential
	MATRIXX S = MATRIXX::Zero(MAR+MBC, MAC+MBC);

	S.block(0, 0, MAR, MAC) = MatrixA;
	S.block(0, MAC, MBR, MBC) = MatrixB;

	S *= CtrlTS;
	S = (S.exp()).eval(); // Avoid Aliasing With eval()

	// Separate into F and G
	MatrixA = S.block(0, 0, MAR, MAC);
	MatrixB = S.block(0, MAC, MBR, MBC);
}

void ZeroOrderHold(MATRIX2 & MatrixA, VECTOR2 & MatrixB, const double & CtrlTS)
{
	// Find Rows and Columns
	const static short MAR = 2; const static short MAC = 2;
	const static short MBR = 2; const static short MBC = 1;

	// Form Control Hierarchy, Compute Matrix Exponential
	MATRIXX S = MATRIXX::Zero(MAR+MBC, MAC+MBC);

	S.block(0, 0, MAR, MAC) = MatrixA;
	S.block(0, MAC, MBR, MBC) = MatrixB;

	S *= CtrlTS;
	S = (S.exp()).eval(); // Avoid Aliasing With eval()

						  // Separate into F and G
	MatrixA = S.block(0, 0, MAR, MAC);
	MatrixB = S.block(0, MAC, MBR, MBC);
}

void ZeroOrderHold(MATRIX3 & MatrixA, VECTOR3 & MatrixB, const double & CtrlTS)
{
	// Find Rows and Columns
	const static short MAR = 3; const static short MAC = 3;
	const static short MBR = 3; const static short MBC = 1;

	// Form Control Hierarchy, Compute Matrix Exponential
	MATRIXX S = MATRIXX::Zero(MAR+MBC, MAC+MBC);

	S.block(0, 0, MAR, MAC) = MatrixA;
	S.block(0, MAC, MBR, MBC) = MatrixB;

	S *= CtrlTS;
	S = (S.exp()).eval(); // Avoid Aliasing With eval()

						  // Separate into F and G
	MatrixA = S.block(0, 0, MAR, MAC);
	MatrixB = S.block(0, MAC, MBR, MBC);
}

void ZeroOrderHold(MATRIX4 & MatrixA, VECTOR4 & MatrixB, const double & CtrlTS)
{
	// Find Rows and Columns
	const static short MAR = 4; const static short MAC = 4;
	const static short MBR = 4; const static short MBC = 1;

	// Form Control Hierarchy, Compute Matrix Exponential
	MATRIXX S = MATRIXX::Zero(MAR+MBC, MAC+MBC);

	S.block(0, 0, MAR, MAC) = MatrixA;
	S.block(0, MAC, MBR, MBC) = MatrixB;

	S *= CtrlTS;
	S = (S.exp()).eval(); // Avoid Aliasing With eval()

	// Separate into F and G
	MatrixA = S.block(0, 0, MAR, MAC);
	MatrixB = S.block(0, MAC, MBR, MBC);
}