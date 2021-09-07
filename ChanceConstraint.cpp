//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// VehicleBase Class Functions for DLL

// To be done:
// See VehicleBase.h

/* Some Notes:
// See VehicleBase.h

*/

#include "ChanceConstraint.h"

void ChanceConstraint::GetIndex(const EIGENFLOAT & v, const VECTORX & ax, VECTORX::Index & i1, VECTORX::Index & i2) const
{
	VECTORX::Index i;
	(ax.array() - v).abs().minCoeff(&i); // Got smallest index

	// Which side is the second smallest?
	if (i == 0)
	{
		i1 = i;
		i2 = i+1;
	}
	else if (i == ax.size()-1)
	{
		i1 = i-1;
		i2 = i;
	}
	else
	{
		const EIGENFLOAT & t1 = ax(i-1);
		const EIGENFLOAT & t2 = ax(i+1);

		if (t1 < v && v < i)
		{
			i1 = i-1;
			i2 = i;
		}
		else
		{
			i1 = i;
			i2 = i+1;
		}
	}
}

EIGENFLOAT ChanceConstraint::Interp1D(const EIGENFLOAT & x, const VECTORX & ax, const VECTORX & y) const // Matlab syntax is (ax, y, x)
{
	// Check if precisely on the axis
	VECTORX::Index Index;
	if (( (ax.array() - x).abs() < 1e-12).any())
	{
		(ax.array() - x).abs().minCoeff(&Index);
		return y(Index);
	}

	// Else not on the axis
	VECTORX::Index i1, i2;
	GetIndex(x, ax, i1, i2);

	return y(i1) + (x - ax(i1)) * (y(i2)-y(i1)) / (ax(i2)-ax(i1));
}

VECTORX ChanceConstraint::Interp1D(const VECTORX & x, const VECTORX & ax, const VECTORX & y) const // Matlab syntax is (ax, y, x)
{
	VECTORX vec(x.size());
	for (int i = 0; i < x.size(); ++i)
	{
		vec(i) = Interp1D(x(i), ax, y);
	}

	return vec;
}

EIGENFLOAT ChanceConstraint::GetChanceConstr(const EIGENFLOAT & alpha, EIGENFLOAT i) const
{
	static bool NotRead = 1;
	static std::string file;

	if (NotRead)
	{
		if (CtrlTS != 1.0 && CtrlTS != 0.2 && CtrlTS != 0.5 && CtrlTS != 0.4) throw std::runtime_error("Input a valid CtrlTS for chance constraints.");
		if (CtrlTS == 1.0 && Type == Gap) file = "_dcc_10.csv";
		else if (CtrlTS == 0.2 && Type == Gap) file = "_dcc_02.csv";
		else if (CtrlTS == 0.4 && Type == Position) file = "_scc.csv";
		else if (CtrlTS == 0.2 && Type == Position) file = "_scc_02.csv";
		else if (CtrlTS == 0.5 && Type == Position) file = "_scc_05.csv";
		else if (CtrlTS == 1.0 && Type == Position) file = "_scc_10.csv";
		else throw std::runtime_error("Input a valid Type for chance constraints.");

		NotRead = 0;
	}

	static const MATRIXX d_r_MAT = ReadMatrix<MATRIXX>(file);
	static const RVECTORX alpha_ax = ReadMatrix<RVECTORX>("_alpha.csv");
	static const RVECTORX i_ax = ReadMatrix<RVECTORX>("_i.csv");

	if (!alpha_ax.size() || !i_ax.size() || !d_r_MAT.size()) throw std::runtime_error("Did not load chance constraints.");
	if (i > 10./CtrlTS) i = 10./CtrlTS; // Data only for 10s ahead
	if (i < 1) throw std::runtime_error("Accessing invalid index in chance constraint.\n");
	if (alpha > alpha_ax(97)) throw std::runtime_error("Accessing invalid index in chance constraint.\n");
	if (alpha < alpha_ax(0)) throw std::runtime_error("Accessing invalid index in chance constraint.\n");
	
	// Check if precisely on one of the axes
	// d_r_MAT is alpha x i
	VECTORX::Index Index;
	if ((i_ax.array() == i).any())
	{
		Index = i-1;
		return Interp1D(alpha, alpha_ax, d_r_MAT.col(Index));
	}
	else if (( (alpha_ax.array() - alpha).abs() < 1e-6).any())
	{
		(alpha_ax.array() - alpha).abs().minCoeff(&Index);
		return Interp1D(i, i_ax, d_r_MAT.row(Index));
	}

	// Sits within 4 quadrant points
	VECTORX::Index alpha_i1, alpha_i2, i_i1, i_i2;
	GetIndex(alpha, alpha_ax, alpha_i1, alpha_i2);
	GetIndex(i, i_ax, i_i1, i_i2);

	const EIGENFLOAT & q11 = d_r_MAT(alpha_i1, i_i1);
	const EIGENFLOAT & q21 = d_r_MAT(alpha_i2, i_i1);
	const EIGENFLOAT & q12 = d_r_MAT(alpha_i1, i_i2);
	const EIGENFLOAT & q22 = d_r_MAT(alpha_i2, i_i2);

	MATRIX2 Q;
	VECTOR2 Y;
	VECTOR2 X;

	Q(0,0) = q11; Q(0,1) = q12;
	Q(1,0) = q21; Q(1,1) = q22;
	X(0) = i_ax(i_i2) - i; 
	X(1) = i - i_ax(i_i1);
	Y(0) = alpha_ax(alpha_i2) - alpha;
	Y(1) = alpha - alpha_ax(alpha_i1);

	return 1./(i_ax(i_i2)-i_ax(i_i1))/(alpha_ax(alpha_i2)-alpha_ax(alpha_i1)) * X.transpose() * Q * Y;
}