//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// MPCVehicle Inherited Class Functions for DLL

#include "MPCVehicle.h"

// Build QP Matrices
void MPCVehicle::BuildQPStructure()
{
	// Control Matrices
	if (UseKinematic) {
		ACtrlMatrix = (MATRIXX(3,3) << 1, 1*CtrlTS, 0.25*CtrlTS*CtrlTS,
									   0, 1, 0.50*CtrlTS,
									   0, 0, 0.00).finished();
		BCtrlMatrix = (VECTORX(3) << 0.25*CtrlTS*CtrlTS,
									   0.50*CtrlTS,
									   1.00).finished();
	}
	else {
		ACtrlMatrix = (MATRIXX(3,3) << 0, 1, 0,
										0, 0, 1,
										0, 0, -TimeConstant).finished();
		BCtrlMatrix = (VECTORX(3) << 0,
									0,
									TimeConstant).finished();

		ZeroOrderHold(ACtrlMatrix, BCtrlMatrix, CtrlTS);
	}
	#if USESPEEDREF
	CCtrlMatrix = (MATRIXX(3,3) << 0, 0, 0,
									0, 0, 0, // Dynamic
									0, 0, 1).finished();
	DCtrlMatrix = (VECTORX(3) << 0,
								0,
								0).finished();
	#else
	CCtrlMatrix = (MATRIXX(3,3) << 1, THeadway, 0,
									0, 1, 0,
									0, 0, 1).finished();
	DCtrlMatrix = (VECTORX(3) << 0,
								0,
								0).finished();
	#endif

	NumStates = ACtrlMatrix.cols();
	NumInputs = BCtrlMatrix.cols();
	NumSoftCon = rhoMatrix.size();
	NumOutputs = CCtrlMatrix.rows();
	NumBinaries = 0;

	static const int InitialSize = 4;
	static const int StageSize = 9 - 2*TURNOFFWHEEL;
	static const int EndSize = 6 - 2*TURNOFFWHEEL - 1*TURNOFFTERMINAL;

	const MATRIX3 & A = ACtrlMatrix;
	const VECTOR3 & B = BCtrlMatrix;

	WOTMatrix = (MATRIXX(2,2) << 
		0.28500, 2.0000,
		-0.12075, 4.8300
	).finished(); // AND constraint on WOT

	AccelerationMax = WOTMatrix(0, 0)*(WOTMatrix(1, 1) - WOTMatrix(0, 1)) / (WOTMatrix(0, 0) - WOTMatrix(1, 0)) + WOTMatrix(0, 1);
	AccelerationMin = -4.50;

	DMin = 2.0;
	DRef = DMin;

	// Cost matrices for MPC
	#if USESPEEDREF
	QMatrix = (MATRIXX(3,3) << QG, 0, 0,
								0, 1, 0,
								0, 0, QA).finished();
	#else
	QMatrix = (MATRIXX(3,3) << QG, 0, 0,
								0, 10, 0,
								0, 0, QA).finished();
	#endif
	
	PMatrix = QMatrix;
	RMatrix = (MATRIXX(1,1) << QA).finished();

	rhoMatrix << 1e7, 1e10, 1e3, 1e7;
	eopt = VECTORX::Zero(4);

	// Cost matrices
	Psi = Eigen::kroneckerProduct(MATRIXX::Identity(N, N), RMatrix);
	Omega = Eigen::kroneckerProduct(MATRIXX::Identity(N-1, N-1), QMatrix);
	blkdiag_inplace(Omega, PMatrix);

	// State matrices
	MATRIXX PhiP = MATRIXX::Identity(N*NumStates, NumStates);

	for (int i = 1; i < N; ++i)
	{
		PhiP.block(i*NumStates, 0, NumStates, NumStates) =
		A * PhiP.block((i-1)*NumStates, 0, NumStates, NumStates);
	}

	Phi = PhiP * A;

	MATRIXX GammaP = MATRIXX::Identity(N*NumStates, N*NumStates);

	for (int j = N-2; j > -1; --j)
	{
		GammaP.block((j+1)*NumStates, j*NumStates, (N-j-1)*NumStates, NumStates) = 
		GammaP.block((j+1)*NumStates, (j+1)*NumStates, (N-j-1)*NumStates, NumStates) * A;
	}

	Gamma = GammaP * Eigen::kroneckerProduct(MATRIXX::Identity(N, N), B);

	#if !USESPEEDREF
	// Finalize G F T Matrices
	bC = Eigen::kroneckerProduct(MATRIXX::Identity(N, N), CCtrlMatrix);

	GMatrix = 2 * (Psi + (Gamma.transpose() * (bC.transpose()*Omega*bC) * Gamma));
	FMatrix = 2 * (Gamma.transpose() * (bC.transpose()*Omega*bC) * Phi);
	TMatrix = 2 * (Gamma.transpose() * bC.transpose()*Omega);

	GA = blkdiag(GMatrix, MATRIXX::Zero(NumSoftCon, NumSoftCon));
	#endif

	CFSuff = rhoMatrix;

	// Form Upsilon and CSigma
	static const VECTORX E0 = (VECTORX(InitialSize) << -1, 1, 1, 1).finished();
	static const VECTORX EI = (VECTORX(StageSize) << -1, 1, 0, 0, 0, 1, 1
		#if !TURNOFFWHEEL
		,0, 0
		#endif
	).finished();

	static const MATRIXX Upsilon0 = MATRIXX::Zero(InitialSize, NumSoftCon);
	static const MATRIXX UpsilonI = (MATRIXX(StageSize, NumSoftCon) << 
		0, 0, 0, 0,
		0, 0, 0, 0,
		-1, 0, 0, 0,
		0, 0, -1, 0,
		0, -1, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0
		#if !TURNOFFWHEEL
		,0, 0, 0, -1,
		0, 0, 0, -1
		#endif
	).finished();

	static const MATRIXX UpsilonN = (MATRIXX(EndSize, NumSoftCon) <<
		-1, 0, 0, 0,
		0, 0, -1, 0,
		0, -1, 0, 0
		#if !TURNOFFTERMINAL
		,-1, 0, 0, 0
		#endif
		#if !TURNOFFWHEEL
		,0, 0, 0, -1
		,0, 0, 0, -1
		#endif
	).finished();

	Upsilon = (MATRIXX(InitialSize + (N-1)*StageSize + EndSize, NumSoftCon) <<
		Upsilon0,
		Eigen::kroneckerProduct(VECTORX::Ones(N-1), UpsilonI),
		UpsilonN
	).finished();

	CurlySigma = E0;
	blkdiag_inplace(CurlySigma, Eigen::kroneckerProduct(MATRIXX::Identity(N-1, N-1), EI));
	CurlySigma.conservativeResize(CurlySigma.rows() + EndSize, Eigen::NoChange);
	CurlySigma.block(CurlySigma.rows()-EndSize, 0, EndSize, CurlySigma.cols()) = MATRIXX::Zero(EndSize, CurlySigma.cols());

	// CPLEX Variable Types
	const static int C = 2; //ILOFLOAT; // 2
	const static int I = 1; //ILOINT; // 1
	const static int Bi = 3; //ILOBOOL; // 3

	CType.conservativeResize(NumSoftCon + N);
	CType.head(NumSoftCon+N) = Eigen::VectorXi::Constant(NumSoftCon+N, C);
}

std::string MPCVehicle::CtoGVType(const Eigen::VectorXi & CType)
{
	// I = 1; // ILOINT I: GRB_INTEGER
	// C = 2; // ILOFLOAT C: GRB_CONTINOUS
	// B = 3; // ILOBOOL B: GRB_BOOL
	static char arr[] = "ICB";

	std::string vtype;
	for (int i = 0; i < CType.size(); ++i)
	{
		vtype += arr[CType(i)-1];
	}

	return vtype;
}

// Augment S for soft constraints function
void MPCVehicle::AugmentS(MATRIXX & S, const MATRIXX & Upsilon, short NumSoftCon)
{
	// S = [S, Upsilon];
	S.conservativeResize(S.rows() + NumSoftCon, S.cols() + Upsilon.cols());
	S.block(
		0, S.cols() - Upsilon.cols(),
		S.rows() - NumSoftCon, Upsilon.cols()
	) = Upsilon;

	// S = [S; zeros(NumSoftCon, S.cols())];
	S.block(
		S.rows() - NumSoftCon, 0,
		NumSoftCon, S.cols()
	) = MATRIXX::Zero(NumSoftCon, S.cols());

	// Assign Negative Diagonal to S
	S.block(
		S.rows() - NumSoftCon, S.cols() - NumSoftCon, 
		NumSoftCon, NumSoftCon
	) = -MATRIXX::Identity(NumSoftCon, NumSoftCon);
}

// Augment S for soft constraints, MIMPC function
void MPCVehicle::MIAugmentS(MATRIXX & S, const MATRIXX & Upsilon, const MATRIXX & CurlyB,
	short NumSoftCon, short NumBinaries)
{
	// S = [S, Upsilon, CurlyB];
	S.conservativeResize(S.rows() + NumSoftCon, S.cols() + Upsilon.cols() + NumBinaries);
	S.block(
		0, S.cols() - NumBinaries - Upsilon.cols(),
		S.rows() - NumSoftCon, Upsilon.cols()
	) = Upsilon;
	S.block(
		0, S.cols() - NumBinaries,
		S.rows() - NumSoftCon, NumBinaries
	) = CurlyB;

	// S = [S; zeros(NumSoftCon, S.cols())];
	S.block(
		S.rows() - NumSoftCon, 0,
		NumSoftCon, S.cols()
	) = MATRIXX::Zero(NumSoftCon, S.cols());

	// Assign Negative Diagonal to Sa
	S.block(
		S.rows() - NumSoftCon, S.cols() - NumSoftCon - NumBinaries, 
		NumSoftCon, NumSoftCon
	) = -MATRIXX::Identity(NumSoftCon, NumSoftCon);
}

// Augment Xi for soft constraints function
MATRIXX MPCVehicle::AugmentXi(const MATRIXX & W, const MATRIXX & cC, const VECTORX & x0, short NumSoftCon)
{
	MATRIXX Xi = cC + (W*x0);

	// Xi = [Xi; zeros(NumSoftCon, 1)]
	Xi.conservativeResize(Xi.rows() + NumSoftCon, Eigen::NoChange);
	Xi.block(
		Xi.rows() - NumSoftCon, 0,
		NumSoftCon, 1
	) = MATRIXX::Zero(NumSoftCon, 1);

	return Xi;
}

// Augment cf for soft constraints function
MATRIXX MPCVehicle::Augmentcf(const MATRIXX & FMatrix, const MATRIXX & TMatrix, const MATRIXX & cR, const VECTORX & CFSuff, const VECTORX & x0)
{
	MATRIXX cf = (FMatrix*x0) - (TMatrix*cR);

	// cf = [cf; CFSuff];
	cf.conservativeResize(cf.rows() + CFSuff.rows(), Eigen::NoChange);
	cf.block(
		cf.rows() - CFSuff.rows(), 0,
		CFSuff.rows(), 1
	) = CFSuff;

	return cf;
}

void MPCVehicle::Solve() {
	#if USESPEEDREF
	// Finalize G F T Matrices
	GMatrix = 2 * (Psi + (Gamma.transpose() * (bC.transpose()*Omega*bC) * Gamma));
	FMatrix = 2 * (Gamma.transpose() * (bC.transpose()*Omega*bC) * Phi);
	TMatrix = 2 * (Gamma.transpose() * bC.transpose()*Omega);

	GA = blkdiag(GMatrix, MATRIXX::Zero(NumSoftCon, NumSoftCon));
	#endif

	// Build constraints
	const VECTORX & x0 = Z;
	const VECTORX & linslp = ConstraintMix.col(0);
	const VECTORX & linmax = ConstraintMix.col(1);

	const int InitialSize = 4;
	const int StageSize = 9 - 2*TURNOFFWHEEL;
	const int EndSize = 6 - 2*TURNOFFWHEEL - 1*TURNOFFTERMINAL;

	MATRIXX cfa, Sa, Xia;

	// Build M
	const double & m1 = linslp(0);
	const double & m2 = linslp(1);
	const double & lp = linslp(2);

	const double & b1 = linmax(0);
	const double & b2 = linmax(1);
	const double & lm = linmax(2);

	const MATRIXX D0 = (MATRIXX(InitialSize, NumStates) << 
		0,  0, 0,
		0,  0, 0,
		0,-m1, 0,
		0,-m2, 0
	).finished();

	const MATRIXX b0 = (MATRIXX(InitialSize, 1) << 
		-AccelerationMin, // -ue \leq -umin
		AccelerationMax, // ue \leq umax
		b1, // m1(-v) + ue - Mb \leq b1
		b2 + BigM // m2(-v) + ue + Mb \leq b2 + M
	).finished(); 

	const MATRIXX M = (MATRIXX(StageSize, NumStates) << 
		0,  0, 0,
		0,  0, 0,
		1,  0, 0,
		0, -1, 0,
		0,  1, 0,
		0,-m1, 0,
		0,-m2, 0
		#if !TURNOFFWHEEL
		,0,-m1, 1,
		0,-m2, 1
		#endif
	).finished();

	const MATRIXX b = (MATRIXX(StageSize, 1) << 
		-AccelerationMin, // -ue \leq -umin
		AccelerationMax, // ue \leq umax
		0, // -d \leq -DMin
		-VelocityMin, // -v \leq VMin
		VelocityMax, // v \leq VMax
		b1, // m1(-v) + ue \leq b1
		b2 + BigM // m2(-v) + ue \leq b2
		#if !TURNOFFWHEEL
		,b1, // m1(-v) + aw \leq b1
		b2 + BigM // m2(-v) + aw \leq b2
		#endif
	).finished(); 

	const MATRIXX MN = (MATRIXX(EndSize, NumStates) << 
		1,  0, 0,
		0, -1, 0,
		0,  1, 0
		#if !TURNOFFTERMINAL
		,1, -lp, 0
		#endif
		#if !TURNOFFWHEEL
		,0,-m1, 1,
		0,-m2, 1
		#endif
	).finished();

	const MATRIXX bN = (MATRIXX(EndSize, 1) << 
		0, // s \leq rc - dmin
		-VelocityMin, // -v \leq VMin
		VelocityMax // v \leq VMax
		#if !TURNOFFTERMINAL
		,lm // -d + mu(-v) \leq xi
		#endif
		#if !TURNOFFWHEEL
		,b1, // m1(-v) + uw \leq b1
		b2 + BigM // m2(-v) + uw \leq b2
		#endif
	).finished(); 

	// Prediction stage reference and constraints
	#if USESPEEDREF
	const MATRIXX cR = 
		Eigen::kroneckerProduct(trajP_prob.block(1, 0, N, 1) - MATRIXX::Constant(N, 1, DRef), (MATRIXX(NumOutputs, 1) << 1, 0, 0).finished()) + 
		Eigen::kroneckerProduct(VRef, (MATRIXX(NumOutputs, 1) << 0, 1, 0).finished());
	#else
	const MATRIXX cR = Eigen::kroneckerProduct(trajP_prob.block(1, 0, N, 1) - MATRIXX::Constant(N, 1, DRef), (MATRIXX(NumOutputs, 1) << 1, 0, 0).finished()); // c(0) = r_prob - DRef;
	#endif

	MATRIXX cM = (MATRIXX(InitialSize + (N-1)*StageSize, (N-1)*NumStates) <<
		MATRIXX::Zero(InitialSize, (N-1)*NumStates), 
		Eigen::kroneckerProduct(MATRIXX::Identity(N-1, N-1), M)
	).finished();
	blkdiag_inplace(cM, MN);

	static MATRIXX cD = MATRIXX::Zero(StageSize, 1);
	cD(2, 0) = 1; // b(2) = r_wc - DMin - Dcc(x);

	MATRIXX ccDi = MATRIXX::Zero(N-1, 1);
	MATRIXX ccDN = MATRIXX::Zero(EndSize, 1);

	#if USECC
	const VECTORX Alpha = VECTORX::LinSpaced(N, 0.999999, 0.5);

	int i;
	if (Mode == FOLLOWINGHUMAN) {
		for (i = 1; i < N; ++i) {
			ccDi(i-1) = (trajP_wc(i) - DMin - Chance.GetChanceConstr(Alpha(i-1), (double)i) < trajP_wc(0)) ? trajP_wc(i) - DMin : trajP_wc(i) - DMin - Chance.GetChanceConstr(Alpha(i-1), (double)i);
			ccDi(i-1) = fmin(ccDi(i-1), conS(i-1));
		}
		ccDN(0, 0) = (trajP_wc(i) - DMin - Chance.GetChanceConstr(Alpha(i-1), (double)i) < trajP_wc(0)) ? trajP_wc(i) - DMin : trajP_wc(i) - DMin - Chance.GetChanceConstr(Alpha(i-1), (double)i); // r_wc - DMin - Dcc(x)
		ccDN(0, 0) = fmin(ccDN(0, 0), conS(N-1));
	}
	else {
		for (i = 1; i < N; ++i) {
			ccDi(i-1) = trajP_wc(i) - DMin;
			ccDi(i-1) = fmin(ccDi(i-1), conS(i-1));
		}
		ccDN(0, 0) = trajP_wc(i) - DMin; // r_wc - DMin
		ccDN(0, 0) = fmin(ccDN(0, 0), conS(N-1));
	}
	#endif

	static MATRIXX cV = MATRIXX::Zero(StageSize, 1);
	cV(4, 0) = 1; // b(4) = conV(x)

	MATRIXX cVN = MATRIXX::Zero(EndSize, 1);
	cVN(2, 0) = VelocityMax;

	const MATRIXX cC = (MATRIXX(InitialSize + (N-1)*StageSize + EndSize, 1) <<
		b0,
		Eigen::kroneckerProduct(MATRIXX::Ones(N-1, 1), b) + Eigen::kroneckerProduct(ccDi, cD), //Eigen::kroneckerProduct(MATRIXX::Ones(N-1, 1), b) + Eigen::kroneckerProduct(ccDi, cD) + Eigen::kroneckerProduct(conS.block(0, 0, N-1, 1), cV),
		bN + ccDN //+ cVN
	).finished();

	// Build D Matrix
	const MATRIXX D = (MATRIXX(InitialSize + (N-1)*StageSize + EndSize, NumStates) << 
		D0,
		MATRIXX::Zero((N-1)*StageSize + EndSize, NumStates)
	).finished();

	// Build QP Matrices
	Sa = (cM * Gamma) + CurlySigma;
	if (!NumBinaries)
	{
		AugmentS(Sa, Upsilon, NumSoftCon);
	}
	else
	{
		MIAugmentS(Sa, Upsilon, CurlyB, NumSoftCon, NumBinaries);
	}

	Xia = AugmentXi(-((cM * Phi) + D), cC, x0, NumSoftCon);
	cfa = Augmentcf(FMatrix, TMatrix, cR, CFSuff, x0);

	// Call Solver
	#if USEGUROBI
	VECTORX U0 = Gurobi.Solve(QMOD*GA, cfa, Sa, Xia, NumSoftCon, NumBinaries, CtoGVType(CType));
	#else
	VECTORX U0 = CPLEXInterface(QMOD*GA, cfa, Sa, Xia, CType, NumSoftCon, NumBinaries);
	#endif

	// Postprocess
	if (U0.size() == 0)
	{
		throw std::runtime_error("Gurobi did not return a solution.");
	}
	else 
	{
		U = U0.head(N).array();
		eopt = U0.segment(N, NumSoftCon);
	}
}

void MPCVehicle::AppendComms(MATRIXX & OtherTrajectory)
{
	// Initialize Local Variables
	int OtherCols = OtherTrajectory.cols();
	int NumAppendedCols = N+1 - OtherCols;
	VECTORX ZOld = OtherTrajectory.col(OtherCols-1);
	VECTORX UZero = VECTORX::Zero(NumAppendedCols);

	// Find Appended Dynamics if Other Vehicle Applies No Input
	MATRIXX AppendedMat(NumStates, NumAppendedCols+1);
	SimulateDynamics(AppendedMat, ZOld, UZero);
	
	// Append Matrices
	OtherTrajectory.conservativeResize(Eigen::NoChange, N+1);
	OtherTrajectory.block(0, N+1-NumAppendedCols, NumStates, NumAppendedCols) = AppendedMat.block(0, 1, NumStates, NumAppendedCols);
}

void MPCVehicle::AppendComms(VECTORX & OtherTrajectory, const double & VPrev)
{
	// Initialize Local Variables
	int OtherCols = OtherTrajectory.size();
	int NumAppendedCols = N+1 - OtherCols;
	VECTOR3 ZOld;
	ZOld(0) = OtherTrajectory(OtherCols-1);
	ZOld(1) = VPrev;
	ZOld(2) = 0.;

	VECTORX UZero = VECTORX::Constant(NumAppendedCols, 0.);

	// Find Appended Dynamics if Other Vehicle Applies No Input
	MATRIXX AppendedMat(NumStates, NumAppendedCols+1);
	SimulateDynamics(AppendedMat, ZOld, UZero);

	// Append Matrices
	OtherTrajectory.conservativeResize(N+1);
	VECTORX AppendedVec = AppendedMat.row(0);
	OtherTrajectory.tail(NumAppendedCols) = AppendedVec.segment(1, NumAppendedCols);
}

short MPCVehicle::GetFinite(const VECTORX & vec, short n)
{
	if (n) {
		if (!std::isfinite(vec(n))) GetFinite(vec, n-1);
		return n;
	}
	else {
		throw std::runtime_error("Error with trajP_wc and terminal constraint.");
	}
}

void MPCVehicle::FormVelocityCon(const SensorInfo & S) {
	// Initialize constraint vector and predict ego's trajectory
	conS = VECTORX::Constant(N, 2000);

	#if SPEEDLIMITDISTANCECAL && !USESPEEDREF
	MATRIXX Forward = MATRIXX::Zero(NumStates, N+1);
	SimulateDynamics( Forward, Z, VECTORX::Constant(N, Z(2)) );

	// Speed Limits
	// if (S.SpeedLimitDistance > -5.0 && S.SpeedLimitDistance < SPEEDLIMITDISTANCECAL)
	// {
	// 	const double prefA = -2.0;
	// 	const double deltaX = abs((S.SpeedLimitValue*S.SpeedLimitValue - VelocityMax*VelocityMax) / (2.0*prefA));

	// 	// Build constraint vector
	// 	if (VelocityMax < S.SpeedLimitValue) // Will be speeding up
	// 	{
	// 		const double prefA = AccelerationMax;
	// 		const double deltaX = abs((S.SpeedLimitValue*S.SpeedLimitValue - VelocityMax * VelocityMax) / (2.0*prefA));
	// 		for (int i = 0; i < N; ++i)
	// 		{
	// 			if (Forward(0, i) >= S.SpeedLimitDistance + deltaX) conV(i) = S.SpeedLimitValue;
	// 			else if (Forward(0, i) < S.SpeedLimitDistance) conV(i) = VelocityMax; 
	// 			else if (Forward(0, i) < S.SpeedLimitDistance + deltaX) conV(i) = pow(VelocityMax*VelocityMax + 2.0*prefA*(Forward(0, i) - S.SpeedLimitDistance), 0.5);
	// 			else throw std::runtime_error("Should have a condition met?");
	// 		}
	// 	}
	// 	else // Will be slowing down
	// 	{
	// 		const double prefA = -2.0;
	// 		const double deltaX = abs((S.SpeedLimitValue*S.SpeedLimitValue - VelocityMax * VelocityMax) / (2.0*prefA));
	// 		for (int i = 0; i < N; ++i) 
	// 		{
	// 			if (Forward(0, i) <= S.SpeedLimitDistance - deltaX) conV(i) = VelocityMax;
	// 			else if (Forward(0, i) < S.SpeedLimitDistance) conV(i) = pow(VelocityMax*VelocityMax + 2.0*prefA*(deltaX - (S.SpeedLimitDistance-Forward(0, i))), 0.5);
	// 			else if (Forward(0, i) >= S.SpeedLimitDistance || S.DriverDataVehVelocity*CtrlTS >= S.SpeedLimitDistance) conV(i) = S.SpeedLimitValue;
	// 			else throw std::runtime_error("Should have a condition met?");
	// 		}
	// 	}
	// }

	// Traffic Lights - Check first light is not green and within range
	if (S.Lights.size() > 0) {
		if (S.Lights.front().SignalDistance > -5.0 && S.Lights.front().SignalDistance < SPEEDLIMITDISTANCECAL && S.Lights.front().SignalState != GREEN) {
			const double prefA = -2.0;
			const double deltaX = abs((0 - VelocityMax*VelocityMax) / (2.0*prefA));

			long SS = S.Lights.front().SignalState;
			if (SS == AMBER && S.DriverDataVehVelocity*(5.0 - S.Lights.front().SignalStateStart) < S.Lights.front().SignalDistance) SS = RED; // Cannot cruise through light

			if (SS == RED) {
				const double & d0 = DRef;
				for (int i = 0; i < N; ++i) 
				{
					conS(i) = S.Lights.front().SignalDistance - 5.0;
					// if (Forward(0, i) <= S.Lights.front().SignalDistance - deltaX) conS(i) = VelocityMax;
					// else if (S.Lights.front().SignalDistance < d0 || S.DriverDataVehVelocity*CtrlTS + d0 >= S.Lights.front().SignalDistance) conS(i) = 0;
					// else if (Forward(0, i) < S.Lights.front().SignalDistance) conS(i) = pow(VelocityMax*VelocityMax + 2.0*prefA*(deltaX - (S.Lights.front().SignalDistance-Forward(0, i))), 0.5);
					// else if (Forward(0, i) >= S.Lights.front().SignalDistance) conS(i) = 0;
					// else conS(i) = 0;
				}
			}
		}
		#endif

		#if USESPEEDREF
		const double prefA = -2.0;
		const double deltaX = abs((0 - VelocityMax*VelocityMax) / (2.0*prefA)) + DRef;
		
		// Traffic Lights - Check first light is within range
		if (S.Lights.front().SignalDistance > -5.0 && S.Lights.front().SignalDistance < deltaX && S.Lights.front().SignalState != GREEN) {
			// Deal with amber light
			long SS = S.Lights.front().SignalState;
			if (SS == AMBER && S.DriverDataVehVelocity*(5.0 - S.Lights.front().SignalStateStart) < S.Lights.front().SignalDistance) {
				SS = RED; // Cannot cruise through light
			}
			else {
				SS = GREEN;
			}

			if (SS == GREEN) { // Treat as green light
				CCtrlMatrix(1,0) = 0;
				CCtrlMatrix(1,1) = 1;

				bC = Eigen::kroneckerProduct(MATRIXX::Identity(N, N), CCtrlMatrix);

				VRef = conS;
			}
			else { // Treat as red light
				CCtrlMatrix(1,0) = 1;
				CCtrlMatrix(1,1) = (2.0*prefA);

				bC = Eigen::kroneckerProduct(MATRIXX::Identity(N, N), CCtrlMatrix);

				VRef = VECTORX::Constant(N, S.Lights.front().SignalDistance);
			}
		} // S.Lights.front().SignalDistance > -5.0 && S.Lights.front().SignalDistance < deltaX
		else { // Treat as green light
			CCtrlMatrix(1,0) = 0;
			CCtrlMatrix(1,1) = 1;

			bC = Eigen::kroneckerProduct(MATRIXX::Identity(N, N), CCtrlMatrix);

			VRef = conS;
		}
		#endif
	}
}

void MPCVehicle::FormConstraintMix(double VPrev)
{
	// First Part is the WOT Matrix
	ConstraintMix.block(0, 0, WOTMatrix.rows(), WOTMatrix.cols()) = WOTMatrix;

	// Decide s
	const double & OAM = AccelerationMin;

	const double Tc = (abs(AccelerationMin - OAM) > FLT_MIN) ? (VPrev - VelocityMax) - (AccelerationMin - OAM) : -1.0;
	const double Vc = OAM * Tc + VPrev;
	
	short i = GetFinite(trajP_wc, N);
	const double s2 = trajP_wc(i) - DMin;

	const double s1 = (VelocityMax > VPrev && Vc > 0 && Tc > 0) ?
		s2 + 0.5 * (pow(VPrev - VelocityMax, 2) / (AccelerationMin - OAM)) :
		s2 - 0.5 * (pow(VPrev, 2) / OAM - pow(VelocityMax, 2) / AccelerationMin);

	// Decide v
	if (abs(VPrev - VelocityMax) < 1e-12) VPrev = VelocityMax-0.5;
	const double v1 = VelocityMax;
	const double v2 = (OAM > AccelerationMin) ? VPrev : VPrev * pow(AccelerationMin / OAM, 0.5);

	// Finish Constraint Mix
	const double Mu = (s2 - s1) / (v2 - v1);
	const double Xi = s2 - (Mu * VPrev);

	ConstraintMix(2, 0) = Mu;
	ConstraintMix(2, 1) = Xi;
}

// Returns Trajectory as nxN+1 Matrix, Given UVector is of size pxN
void MPCVehicle::SimulateDynamics(MATRIXX & Trajectory, VECTORX Z, VECTORX UVector)
{
	Trajectory.block(0, 0, NumStates, 1) = Z;

	// Calculate Dynamics
	for (int x = 1; x <= UVector.size(); ++x)
	{
		const double & u = (abs(UVector(x-1)) > 0.01) ? UVector(x-1) : 0.;
		Trajectory.block(0, x, NumStates, 1) = ACtrlMatrix*Z + BCtrlMatrix*UVector(x-1);
		Z = Trajectory.col(x);

		// Saturate and break if Trajectory is moving backwards or moving at max velocity
		const double VM = (Z(0) > SpeedLimitDistance && SpeedLimitValue > 1) ? SpeedLimitValue : VelocityMax;
		if (Z(1) < 0)
		{
			Z(1) = 0;
			Z(2) = 0;
			UVector.setConstant(0.);

			Trajectory.block(0, x, NumStates, 1) = Z;
		}
		else if (Z(1) > VM)
		{
			Z(1) = VM;
			Z(2) = 0;
			UVector.setConstant(0.);

			Trajectory.block(0, x, NumStates, 1) = Z;
		}
	}
}

/* -------------------------------------------------------------------- */

MPCVehicle::MPCVehicle(
		double CtrlTS, double TimeConstant, double THeadway, int N, 
		double QG, double QA, bool UseKinematic
	) : Chance(CtrlTS, ChanceConstraint::Position),
	CtrlTS(CtrlTS), TimeConstant(TimeConstant), THeadway(THeadway), 
	N(N), QG(QG), QA(QA), UseKinematic(UseKinematic) 
{
	// Logging
	#if MPCLOGGING
	if(!LOGFILE.is_open()) LOGFILE.open(MPCLOGNAME, std::ios::trunc);
	#endif

	// OPENMP
	#if NUMTHREADS != 1 && NUMTHREADS != 0
	static bool setOMP = 1;
	if (setOMP) {
		omp_set_num_threads(NUMTHREADS);
		Eigen::initParallel();
		Eigen::setNbThreads(NUMTHREADS);

		setOMP = 0;
	}
	#endif

	// Build Constant QP Matrices
	BuildQPStructure();
}

MPCVehicle::~MPCVehicle(){}

MATRIXX MPCVehicle::DetermineControlMove(double * Control, const int NumControls, const SensorInfo & S, const std::unordered_map <long, MATRIXX> & Communicator) {
	// Handle sensing information
	VelocityMax = S.DriverDataDesVelocity;
	if (S.NVehs.size() > 0 && (ApplyCatchUpToHuman || S.NVehs.front().NDataVehType != HUMAN)) VelocityMax += CatchUpVelocity;

	// Road Information
	SpeedLimitDistance = S.SpeedLimitDistance;
	SpeedLimitValue = S.SpeedLimitValue;

	// Check if Other Vehicle Exists
	double OtherDistance, RelativeVelocity, OtherVelocity, OtherAcceleration;
	long OtherID;

	if (S.NVehs.size() > 0 && S.NVehs.front().NDataVehID > 0)
	{
		// Relative Vehicle Values
		OtherDistance = S.NVehs.front().NDataRelDistance - S.NVehs.front().NDataVehLength;
		RelativeVelocity = S.NVehs.front().NDataRelVelocity;
		OtherAcceleration = S.NVehs.front().NDataAcceleration;

		// Other Vehicle Values
		OtherID = S.NVehs.front().NDataVehID;
	}
	else // Free Flow Variables Otherwise
	{
		// Relative Vehicle Values
		OtherDistance = 2*S.DriverDataDesVelocity;
		RelativeVelocity = S.DriverDataVehVelocity - S.DriverDataDesVelocity; // OtherVelocity = DesiredVelocity
		VelocityMax -= CatchUpVelocity;
		OtherAcceleration = 2.;

		// Other Vehicle Values
		OtherID = NOVEHICLE;											
	}

	if (OtherDistance > 225) OtherDistance = 225.;
	OtherVelocity = -RelativeVelocity + S.DriverDataVehVelocity; // Relative is Ego - Other

	Z = (VECTORX(3) << 0.0, // For relative coordinates to the ego
					   S.DriverDataVehVelocity,
					   S.DriverDataVehAcceleration).finished();

	OtherZ = (VECTORX(3) << OtherDistance,
							OtherVelocity,
							OtherAcceleration).finished();
	
	// Handle preceding vehicle trajectory
	MATRIXX OtherTrajectory; // Started with the Current State of the Vehicle, Should be nxN+1

	if (CheckComms<MATRIXX>(Communicator, OtherID)) // Following an MPC Vehicle in which Comms have been Sent
	{
		OtherTrajectory = ReceiveComms<MATRIXX>(Communicator, OtherID);
		if (N+1 > OtherTrajectory.cols()) AppendComms(OtherTrajectory); // Other Trajectory should be at least 1 old stage plus N Stages

		trajP_prob = OtherTrajectory.row(0).array() + OtherDistance;
		trajP_wc = trajP_prob;

		Mode = FOLLOWINGCAV;

		if (LOGFILE.is_open()) LOGFILE << "FOLLOWING CAV: ";
	}
	else // Not Following an MPC Vehicle or Opening Part of Simulation
	{
		OtherTrajectory = MATRIXX::Zero(NumStates, N+1);

		SimulateDynamics( OtherTrajectory, OtherZ, VECTORX::Constant(N, OtherZ(2)) );
		trajP_prob = OtherTrajectory.row(0);
		trajP_wc = trajP_prob; // s_PV_panic_vec

		Mode = FOLLOWINGHUMAN;

		if (LOGFILE.is_open()) LOGFILE << "FOLLOWING HUMAN: ";
	}

	// Build time-varying velocity constraint
	FormVelocityCon(S);

	// Robust constraint
	const double & VPrev = OtherTrajectory(1, N); 
	FormConstraintMix(VPrev);

	// Build QP and solve
	Solve();

	// Calculate Future States Given Admissible Control
	MATRIXX Trajectory(NumStates, N+1); 
	if ((S.DriverDataVehVelocity*CtrlTS >= SpeedLimitDistance && SpeedLimitValue < 1.0) || 
		(S.Lights.size() > 0 && S.Lights.front().SignalState == RED && S.DriverDataVehVelocity*CtrlTS >= S.Lights.front().SignalDistance)) // Force slow down at traffic light for safety
	{
		const double & brake = -2.0;
		SimulateDynamics(Trajectory, Z, VECTORX::Constant(N, brake));

		// Return Immediate Move
		Control[0] = Z(1) + brake*0.1;
		Control[1] = brake;
	}
	else
	{
		SimulateDynamics(Trajectory, Z, U);

		// Return Immediate Move
		Control[0] = Z(1) + U(0)*0.1; // Velocity
		Control[1] = U(0); // Acceleration

		if (LOGFILE.is_open()) {
			for (int i = 0; i < 16; ++i) LOGFILE << trajP_prob(i+1) - Trajectory(0,i) << ",\t";
			LOGFILE << "\n";
		}
	}

	return Trajectory;
}

EIGENFLOAT MPCVehicle::GetAcc(const SensorInfo & S, const MATRIXX & Trajectory) {
	return Chance.Interp1D(S.DriverDataTimeStep, (VECTORX(2) << 0., CtrlTS).finished(), (VECTORX(2) << Trajectory(2,0), Trajectory(2,1)).finished()); 
}

int MPCVehicle::CheckTime(const EIGENFLOAT & Time)
{
	static int res = 10;

	int CT = round(Time*res);
	int TS = round(CtrlTS*res);

	return (CT % TS);
}