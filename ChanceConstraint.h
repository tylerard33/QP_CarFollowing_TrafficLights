//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// Chance Constraint Class for MPCVehicle

#pragma once

#include "Eigen.h"

#include <exception>
#include <omp.h>

class ChanceConstraint
{
	private:
	double CtrlTS;
	int Type;

	void GetIndex(const EIGENFLOAT & v, const VECTORX & ax, VECTORX::Index & i1, VECTORX::Index & i2) const;

	public:
	ChanceConstraint(const double & CtrlTS = 1.0, int Type = Position) : CtrlTS(CtrlTS), Type(Type) {}

	EIGENFLOAT Interp1D(const EIGENFLOAT & x, const VECTORX & ax, const VECTORX & y) const; // Matlab syntax is (ax, y, x)
	VECTORX Interp1D(const VECTORX & x, const VECTORX & ax, const VECTORX & y) const; // Matlab syntax is (ax, y, x)
	EIGENFLOAT GetChanceConstr(const EIGENFLOAT & alpha, EIGENFLOAT i) const;

	enum Type {
		Position = 1,
		Gap = 2
	};
};