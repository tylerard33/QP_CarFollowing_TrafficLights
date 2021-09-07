//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// Gurobi Functions for Solving Quadratic Cost MPC

#include "Gurobi.h"

GurobiModel::GurobiModel() : env() {} // If error here, likely licensing issue

void GurobiModel::WriteValue(const std::string & s) {
	std::cout << s << "\n";
}

// Ga, cfa, Sa, Xia == Q, c, A, b
// Note that it is expected Matrices are in Colwise Order
VECTORX GurobiModel::Solve(const MATRIXX & Q, const MATRIXX & c, const MATRIXX & A, const MATRIXX & b,
	const int NumSoftCon, const int NumBinaries, std::string CType)
{
	// Initialize Counting Variables
	const int NumVars = Q.rows();
	const int NumCons = b.rows();

	// Create Model
	GRBModel model = GRBModel(env);

	// Initialize Variables
	static std::vector<double> lb;
	static std::vector<GRBLinExpr> lhs;
	static GRBQuadExpr obj;
	obj = 0.;

	if (CType.empty()) // Assume continuous model if no CType argument supplied
	{
		for (int i = 0; i < NumVars; ++i) CType += "C";
	}

	int Size = lb.size();
	if (Size < NumVars)
	{
		for (int i = 0; i < NumVars - Size; ++i)
		{
			lb.push_back(-GRB_INFINITY); // No Lower Bound
		}
	}

	GRBVar * vars = model.addVars(lb.data(), NULL, NULL, CType.c_str(), NULL, NumVars);

	// Q Matrix:
	GRBVar * rowvars = new GRBVar[Q.size()];
	GRBVar * colvars = new GRBVar[Q.size()];

	for (int j = 0; j < Q.cols(); ++j)
	{
		for (int i = 0; i < Q.rows(); ++i)
		{
			rowvars[i + j*Q.cols()] = vars[i];
			colvars[i + j*Q.rows()] = vars[j];
		}
	}

	obj.addTerms(Q.data(), rowvars, colvars, Q.size());

	// c Vector:
	obj.addTerms(c.data(), vars, NumVars);
	model.setObjective(obj, GRB_MINIMIZE);

	// Add Constraints to Model: A <= b:
	Size = lhs.size();
	if (Size < NumCons)
	{
		for (int i = 0; i < NumCons - Size; ++i)
		{
			lhs.push_back(0.); // Initialize to 0, reset later
		}
	}

	for (int j = 0; j < NumVars; ++j)
	{
		for (int i = 0; i < NumCons; ++i)
		{
			#if REMOVEINF
			if (b(i, 0) >= GRB_INFINITY) continue;
			#endif

			if (A(i, j))
			{
				lhs[i] += A(i, j) * vars[j];
			}
		}
	}

	for (int i = 0; i < NumCons; ++i)
	{
		#if REMOVEINF
		if (b(i, 0) >= GRB_INFINITY) continue;
		#endif
		model.addConstr(lhs[i], GRB_LESS_EQUAL, b(i, 0));
		lhs[i] = 0.; // Reset lhs
	}

	// Set Parameters
	#if TUNE
	static int Count = 0;
	Count++;
	if (Count == 100){
		model.getEnv().set(GRB_DoubleParam_TuneTimeLimit, 30.);
		model.tune();
	}
	else if (Count > 100) std::cin.ignore();
	#endif

	model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	#ifdef NDEBUG
	model.set(GRB_IntParam_DualReductions, 1); // 0 to turn off and read if model is infeasible or unbounded
	#else
	model.set(GRB_IntParam_DualReductions, 0);
	#endif
	model.set(GRB_IntParam_Presolve, 1); // 1 to turn on presolving model
	model.set(GRB_DoubleParam_FeasibilityTol, FEATOL); // 1e-6 default
	model.set(GRB_DoubleParam_OptimalityTol, OPTTOL); // 1e-6 default
	if (NumBinaries)
	{
		model.set(GRB_DoubleParam_Heuristics, 0.);
		model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
	}
	else
	{
		model.set(GRB_IntParam_PreDual, 1);
		model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
	}

	// Extract Solution
	model.optimize();

	int status = model.get(GRB_IntAttr_Status);
	if (status == GRB_OPTIMAL) {}
	else if (status == GRB_INFEASIBLE) WriteValue("Model is infeasible.\n");
	else if (status == GRB_UNBOUNDED) WriteValue("Model is unbounded.\n");
	else if (status == GRB_INF_OR_UNBD) WriteValue("Model is inf or unb.\n");
	else if (status == GRB_CUTOFF) WriteValue("Objective exceeded Cutoff.\n");
	else if (status == GRB_ITERATION_LIMIT) WriteValue("Barrier/Simplex iterations exceeded BarIterLimit.\n");
	else if (status == GRB_NODE_LIMIT) WriteValue("Nodes explored exceeded NodeLimit.\n");
	else if (status == GRB_TIME_LIMIT) WriteValue("Time exceeded TimeLimit.\n");
	else if (status == GRB_SOLUTION_LIMIT) WriteValue("Number of solutions reached SolutionLimit.\n");
	else if (status == GRB_INTERRUPTED) WriteValue("Optimization was terminated by user.\n");
	else if (status == GRB_SUBOPTIMAL) WriteValue("Unable to satisfy optimality tolerances.\n");
	else if (status == GRB_INPROGRESS) WriteValue("Solution in progress.\n");
	else if (status == GRB_USER_OBJ_LIMIT) WriteValue("Objective reached limit specified by user.\n");

	VECTORX Solution(NumVars);
	for (int i = 0; i < NumVars; ++i)
	{
		Solution(i) = vars[i].get(GRB_DoubleAttr_X);
	}

	delete[] vars;
	delete[] rowvars;
	delete[] colvars;

	return Solution;
}