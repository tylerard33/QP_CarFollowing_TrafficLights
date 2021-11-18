# Important Features
## Constructor Inputs
CtrlTS [s] should be the fastest multiple of 0.10 that runs online

VehCat can be HGV or SUV (choosing CAR will default to SUV) - this choice affects the default braking capability of the ego and the asymmetric model used in dynamics calculations for simulation

MPCType should be MIPCCONSERVATIVE or MPCCONSERVATIVE - for HGV or SUV cases respectively

1/Tau should be chosen as determined by the acceleration data for the physical vehicles, or left as the default in simulation

MyAccMin limits the braking capabilities of the Ego - input -0.0 to use max braking

OtherAccMin limits the braking capabilities of the PV - it has a safe range of [-0.78, -8.5]

CommDelayConstant [s] is the delay MPC will compensate to the PV communicated position by linearly extrapolating the communicated message by the sensed velocity of the PV - (e.g. S_PV += Delay*VOther)

THeadway [s] is a control switch between gap tracking and acceleration smoothening - 0.0 to use *old code*, or input a time headway for MPC to track with *new code*

N, QG, QA are overrides to the control horizon, gap weighting, and acceleration weighting in the cost function

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}"/>

# Car-Following MPC Access Methods

Access structures for the MPC C++ code are:

* MPCVehicle
* VissimInfo

The public methods are:

    MPCVehicle::MPCVehicle(CtrlTS, VehCat, MPCType, 1/Tau, MyAccMin, OtherAccMin, CommDelayConstant, THeadway, N, QG, QA)

    MPCVehicle::Update(VissimInfo)

    MyTrajectory = MPCVehicle::DetermineControlMove(Control, NumControls, Communicator)
    MPCVehicle::DetermineControlMove(Control, NumControls, MyTrajectory, MyVN)
    MPCVehicle::DetermineControlMove(Control, NumControls, MyTrajectory, MyVN, OtherTrajectory, OtherVN)

The public variables are:

    MPCVehicle::eopt; // Vector of length 4
    MPCVehicle::constr_mix; // Vector of length 2
    MPCVehicle::wc; // Vector of length N+1
    MPCVehicle::prob; // Vector of length N+1
    
# Compilation
## Dependencies:
* Eigen // Eigen is a header-only library that needs to be added to the include path
* CPLEX // On Windows, this requires linking libs "cplex1280.lib" "concert.lib" "ilocplex.lib" and containing cplex1280.dll in the same directory as the binary file


## Multithreading:
* Enable multithreading for the matrix options with additional command line argument:
    
    -fopenmp

## Code Optimization:
* Enable code optimization for speed during compilation with the compiler flag 

    -O3

* If this is slower than without the flag, run -O2 instead

## Code to run:

    /* -- Preprocessor -- */
    #include "MPCVehicle.h"
    #include "VissimInfo.h"
    #include "Eigen.h"
    #include "Timer.h"

    #define USECONNECTED 1 // 1 for ConnectedMPC, 0 for Unconnected
    #if USECONNECTED
        #define MPCTYPE MPCAGGRESSIVE
        #define STAGES 17 // N
        #define QA 1530
        #define QG 1
    #else 
        #define MPCTYPE MPCCONSERVATIVE
        #define STAGES 16 // N
        #define QA 850
        #define QG 1
    #endif

    // These are explained above
    #define MPCTIMESTEP 1.0 // 1s discretization step
    #define MPCTAUINV 3.6364 // Substitute unique value for Leaf and Mazda
    #define MYACCMIN -0.0 // 0 for default
    #define COMMSDELAY 0.0 
    #define TIMEHEADWAY 0.0 // 0 for original cost function

    /* -- Main code -- */
    int main()
    {
        double Controls[2];

        MPCVehicle Ego(MPCTIMESTEP, SUV, MPCTYPE, MPCTIMESTEP, MPCTAUINV, MYACCMIN, COMMSDELAY, TIMEHEADWAY, STAGES, QA, QG); // This will be Aggressive or Conservative depending on connected or unconnected cases, as defined by USECONNECTED

        VECTORX MyTrajectory(STAGES+1);
        VECTORX OtherTrajectory = // Vector received from communication
        double MyVPrev;
        double OtherVPrev = // Value received from communication

        VissimInfo EgoInfo; // Populate RelDistance, RelVelocity, DriverVelocity, DriverAcceleration, NVehID>0, NVehCategory=SUV - Populate NVehLength if RelDistance is front-to-front, make NVehLength 0 if RelDistance is back-to-front

        Ego.Update(EgoInfo);
        #if USECONNECTED
        Ego.DetermineControlMove(Controls, 2, MyTrajectory, MyVPrev, OtherTrajectory, OtherVPrev);
        #else
        Ego.DetermineControlMove(Controls, 2, MyTrajectory, MyVPrev);
        #endif

        // Communicate MyTrajectory, MyVPrev

        return 0;
    }