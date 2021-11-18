//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// MPCVehicle Inherited Class for DLL

//////////////////////////////////////////////////////////////////////////////
/// Import Libraries and Custom Classes									   ///
//////////////////////////////////////////////////////////////////////////////
# pragma once

#include "CarFollowingDefines.h"
#include "Eigen.h"
#include "Communicator.h"
#include "SensorInfo.h"
#include "ChanceConstraint.h"

#if USEGUROBI
#include "Gurobi.h"
#else
#include "CPLEXInterface.h"
#endif

#include <unordered_map>
#include <limits>
#include <cmath>
#include <string>
#include <fstream>

class MPCVehicle
{	
	private:

	enum MODE {
		FOLLOWINGCAV = 1,
		FOLLOWINGHUMAN = 2
	};

	std::ofstream LOGFILE;

	// Vehicle Control Matrix Declarations
	bool UseKinematic;
	bool RobustCase = 0;
	bool ApplyCatchUpToHuman = 0;
	double CatchUpVelocity = 0.0;
	
	int Mode;

	MATRIXX ACtrlMatrix;
	VECTORX BCtrlMatrix;
	MATRIXX CCtrlMatrix;
	VECTORX DCtrlMatrix;

	MATRIXX WOTMatrix;

	double CtrlTS;
	double TimeConstant;

	short NumSoftCon;
	short NumStates;
	short NumInputs;
	short NumOutputs;
	short NumBinaries;

	double BigM = 0.;

	VECTORX Z;
	VECTORX OtherZ;

	// Vehicle Controller Characteristics Declarations
	int N;
	double QG;
	double QA;
	double THeadway;
	
	// MPC QPMatrix Formation Variables
	MATRIXX QMatrix;
	MATRIXX PMatrix;
	MATRIXX	RMatrix;
	VECTOR4 rhoMatrix; // Default of [1e6; 1e6; 1e6; 1e7]

	// MPC QPStruct Formation Variables
	MATRIXX CurlySigma;
	MATRIXX Upsilon;
	MATRIXX Omega;
	MATRIXX Psi;
	MATRIXX bC;

	MATRIXX Phi;
	MATRIXX Gamma;

	MATRIXX GMatrix;
	MATRIXX FMatrix;
	MATRIXX TMatrix;

	MATRIXX GA;
	VECTORX CFSuff;
	Eigen::VectorXi CType;
	MATRIXX CurlyB;

	Eigen::Matrix<double, 3, 2> ConstraintMix;

	double VelocityMin = 0.;
	double VelocityMax;
	double SpeedLimitDistance;
	double SpeedLimitValue;

	double AccelerationMax;
	double AccelerationMin;

	double DMin;
	double DRef;

	ChanceConstraint Chance;

	#if USEGUROBI
    GurobiModel Gurobi;
	std::string CtoGVType(const Eigen::VectorXi & CType);
	#endif
			 
	void AppendComms(MATRIXX & OtherTrajectory);
	void AppendComms(VECTORX & OtherTrajectory, const double & VPrev);
	void FormConstraintMix(double VPrev);
	void FormVelocityCon(const SensorInfo & V);
	void BuildQPStructure();
	void SimulateDynamics(MATRIXX & Trajectory, VECTORX ZOld, VECTORX UVector);
	short GetFinite(const VECTORX & vec, short n);
	
	void AugmentS(MATRIXX & S, const MATRIXX & Upsilon, short NumSoftCon);
	void MIAugmentS(MATRIXX & S, const MATRIXX & Upsilon, const MATRIXX & CurlyB, short NumSoftCon, short NumBinaries);
	MATRIXX AugmentXi(const MATRIXX & W, const MATRIXX & cC, const VECTORX & x0, short NumSoftCon);
	MATRIXX Augmentcf(const MATRIXX & FMatrix, const MATRIXX & TMatrix, const MATRIXX & cR, const VECTORX & CFSuff, const VECTORX & x0);
	void Solve();

	public:	

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	VECTORX U;
	VECTORX eopt;
	VECTORX trajP_prob;
	VECTORX trajP_wc;
	VECTORX conS;
	VECTORX VRef;

	MPCVehicle(
		double CtrlTS = 1.0, 
		double TimeConstant = 3.6364,
		double THeadway = 1.2, int N = 16, double QG = 1, double QA = 20,
		bool UseKinematic = 0
	);

	~MPCVehicle();

	// Create Public Functions
	MATRIXX DetermineControlMove(double * Control, const int NumControls, const SensorInfo & S, const std::unordered_map <long, MATRIXX> & Communicator);
	EIGENFLOAT GetAcc(const SensorInfo & S, const MATRIXX & Trajectory);
	int CheckTime(const EIGENFLOAT & Time);
};
