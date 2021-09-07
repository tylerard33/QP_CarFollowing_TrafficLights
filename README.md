# QP_CarFollowing_TrafficLights

Model predictive control formulation for a car following controller that follows behind another vehicle with desired time headway *TH*. Constant acceleration prediction for the preceding vehicle is made over the control horizon, so chance constraints are then used to mitigate risk against prediction mismatch.

## Dependencies
Include header and library path for Gurobi optimizer. 9.12 assumed in "Gurobi.h". For example: /IC:/gurobi912/win64/include/ /link /LIBPATH:C:/gurobi912/win64/lib/
Include header path for Eigen3 base and unsupported libraries.

# Usage

#include "CarFollowing/MPCVehicle.h"
#include "CarFollowing/Eigen.h" // typedef of MATRIXX contained here
#include "CarFollowing/SensorInfo.h"

#include <vector>
#include <unordered_map>

int main() {
    MPCVehicle Anti(1.0, 3.6364, 1.2, 16); // MPC mode following behind an unconnected vehicle
    MPCVehicle Conn(1.0, 3.6364, 0.6, 17); // MPC mode following behind a connected MPC vehicle
    
    SensorInfo Sensor; // Data structure containing sensor inputs MPC will consume
    std::unordered_map<long, MATRIXX> Communications; // Map containing long Vehicle ID of surrounding vehicles, MATRIXX structure of their planned intentions

    double Control[2];
  
    /* Populate SensorInfo */
    /* Populate Communications */
  
    MPCVehicle & Ego = (S.NVehs.size() > 0 && CheckComms<MATRIXX>(Communications, S.NVehs.back().NDataVehID)) ? Conn : Anti; // Assumed S.NVehs.back() is the preceding vehicle
    Ego.DetermineControlMove(Control, 2, Sensor, Communications);
  
    const double & VelCommand = Control[0];
    const double & AccCommand = Control[1];
  
    return 0;
}
