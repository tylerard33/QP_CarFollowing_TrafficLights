//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////

#pragma once

// Constants
constexpr auto MPH2MS = 0.44704;

// Multithreading
#define NUMTHREADS 2

// Road Geometry
#define SPEEDLIMITDISTANCECAL 1000 // Nonzero to create time varying speed constraint based on SpeedLimitDistance and SpeedLimitValue
#define USESPEEDREF 0 // 1 to use speed limit references instead of constraints

// Solver and MPC Building
constexpr auto QMOD = 0.5; // 1.0 or 0.5 for coefficient on Q in QP
#define USEGUROBI 1 // 1 to use Gurobi, 0 to use CPLEX
#define WRITEMATRIX 0 // 1 to write matrices from build
#define TESTCONSTRAINT 0 // 1 to simulate trajectory from U and compare to constraints on AU <= b and output from SimDyn

#define USECC 1 // 1 to use chance constraint in S
#define TURNOFFTERMINAL 1 // 1 to turn off terminal constraint on velocity
#define TURNOFFWHEEL 1 // 1 to turn off constraint on acceleration at the wheel

#define MPCLOGNAME "data/MPCLog.txt"
#define MPCLOGGING false