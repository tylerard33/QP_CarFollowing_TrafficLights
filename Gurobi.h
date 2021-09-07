//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// Gurobi Functions for Solving Quadratic Cost MPC

/*
	Error codes found at https://www.gurobi.com/documentation/8.1/refman/error_codes.html
	Status codes found at https://www.gurobi.com/documentation/8.1/refman/optimization_status_codes.html
*/

//////////////////////////////////////////////////////////////////////////////
/// Import Libraries and Custom Classes									   ///
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gurobi_c++.h>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

#include "Eigen.h"

#define TUNE 0 // 1 to find optimization parameters
#define REMOVEINF 1 // 1 to ignore rows with infinity on the constraint
#define FEATOL 1e-3
#define OPTTOL 1e-6

#ifdef _DLL // Autodef from /LD MSVC option
	#ifdef _DEBUG
	#pragma comment(lib, "gurobi_c++mdd2019.lib")
	#else
	#pragma comment(lib, "gurobi_c++md2019.lib")
	#endif
#else
	#ifdef _DEBUG
	#pragma comment(lib, "gurobi_c++mtd2019.lib")
	#else
	#pragma comment(lib, "gurobi_c++mt2019.lib")
	#endif
#endif

#pragma comment(lib, "gurobi91.lib")

/// Solver for Gurobi
class GurobiModel 
{
private:
	GRBEnv env;
	void WriteValue(const std::string & s);
	
public:
	GurobiModel();
	VECTORX Solve(const MATRIXX & Q, const MATRIXX & c, const MATRIXX & A, const MATRIXX & b,
		const int NumSoftCon, const int NumBinaries, std::string CType = "");
};