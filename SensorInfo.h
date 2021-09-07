//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University										   ///
//////////////////////////////////////////////////////////////////////////////
// VissimInfo Public Class for DLL Information Pass Between Vissim and VehicleBase

#pragma once

#include <vector>

// Define Constants to Rename Vissim Types
enum VissimTypes
{
	// VehicleType Definitions
	HUMAN = 101,
	ANTICIPATIVE = 102,
	CONNECTED = 103,
	ANLHORIZON = 104,
	CYCLE = 105,

	// Vehicle Definitions
	NOVEHICLE = 0
};

// Interaction Targets
enum InteractionTargetType {
	 /* no target = 0, real vehicle = 1, signal head = 2,       */
    /*         priority rule = 3, conflict area = 4, reduced speed     */
    /*         area = 5, stop sign = 6, parking lot = 7, PT stop = 8 */
	NOTARGET = 0,
	REALVEHICLE = 1, 
	SIGNALHEAD = 2,
	PRIORITYRULE = 3,
	CONFLICTAREA = 4,
	REDUCEDSPEEDAREA = 5,
	STOPSIGN = 6,
	PARKINGLOT = 7,
	PTSTOP = 8
};

// Traffic Lights
enum VissimLightStatus {
	/* long:   red = 1, amber = 2, green = 3, red/amber = 4, */
    /*         amber flashing = 5, off = 6, green arrow = 7  */
    GREEN = 3,
    AMBER = 2,
    RED = 1
};

enum ANLLightStatus {
	ANLGREEN = 2,
	ANLAMBER = 1,
	ANLRED = 0
};

enum ObjClass {
    CAR = 0,
    SUV = 1,
    HGV = 2,
    PED = 3,
    STOP = 4,
    LIGHT = 5
};

enum TurningIndicator {
	BRAKELIGHT = 2,
	LEFTTURN = 1,
	RIGHTTURN = -1,
	NOTURN = 0,
    BRAKING = 3,
    NORMAL = 4
};

struct LightInfo {
	// Signal Info
	double SignalDistance;					// Distance to Upcoming Signal
	double SignalStateStart;				// Time that the Signal Turned to its Current State
	long SignalState;						// Current State of the Signal - red = 1, amber = 2, green = 3

	long IntacTargetID;
};

struct NVehInfo {
	// Nearby Vehicle Characteristics
	long NDataVehID = 0;					// Neighboring Unique ID
	long NDataVehType = 0;					// Neighboring Vehicle Type		(MPCCon..)
	long NDataVehCategory = 0;				// Neighboring Vehicle Category (CAR, HGV)

	double NDataVehLength = 0.0;			// Length of Neighboring Vehicle
	double NDataVehWidth = 0.0;

	// Nearby Vehicle Kinematic Values
	double NDataRelDistance = 0.0;			// Neighboring Relative Distance (Front end to Front end)
	double NDataRelVelocity = 0.0;			// Neighboring Relative Velocity (Veh Speed - NVeh Speed) => VehSpeed - RelVelocity = NVeh Speed
	double NDataAcceleration = 0.0;

	long NDataTurningIndicator = 0;			// Neighboring Turning Indicator
	double NDataLanePosition = 0.0;			// Neighboring Longitudinal Position from Their Lane
};

struct SensorInfo // DriverInfo
{
	// Init Parameters
	double DriverDataTime = 0.00;			// Total time
	double DriverDataTimeStep = 0.10;		// Simulation Time Step

	long IntacTargetType = 0;				// InteractionTargetType enum
	long IntacTargetID = 0;
	double IntacTargetHeadway = 250.0;		// [m] defaults to 250 if nothing ahead

	// Current Vehicle Characteristics
	long DriverDataVehID = 0;				// Unique ID
	long DriverDataVehType = 0;				// VehicleType
	long DriverDataVehCategory = 0;			// VehicleCategory

	double DriverDataVehLength = 0.0;		// Length
	double DriverDataVehWidth = 0.0;		// Width
	double DriverDataVehWeight = 0.0;		// Weight

	// Current Vehicle Kinematic Values
	double DriverDataVehPosition = 0.0;
	double DriverDataVehVelocity = 0.0;		// Instantaneous Velocity
	double DriverDataVehAcceleration = 0.0;	// Instantaneous Acceleration
	double DriverDataMaxAcceleration = 0.0;
	double DriverDataDesVelocity = 0.0;		// Road Speed Limit

	long DriverDataVehLane = 0;				// Current Lane - +1 Rightmost
	long DriverDataPrefLane = 0;			// Preferred Lane by Routing - Negative to the Right
	double DriverDataLanePosition = 0.0;	// Instantaneous Lateral Position - Front of the Vehicle From the Middle of the Lane
	double DriverDataLaneAngle = 0.0;		// Instantaneous Lane Angle - CCW From the Tangent of the Lane
	double DriverDataDesAngle = 0.0;		// Desired Lane Angle (May need this to follow the lanes)

	// Vehicle Control Values
	double Command = 0.;                    // Engine command u
	double Traction = 0.;                   // Wheel command aw
	double Command_v = 0.;

	// Surrounding Vehicle Information
	std::vector<NVehInfo> NVehs;

	// Road Information
	double SpeedLimitDistance = 2000.0;		// [m] Distance to speed limit/reduced speed area
	double SpeedLimitValue = 100.0;			// [m/s] Value of speed limit/reduced speed area

	double LaneWidth;						// Lane Width
	long NumberOfLanes;						// Number of Lanes of Current Link

	std::vector<LightInfo> Lights;
};