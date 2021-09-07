//////////////////////////////////////////////////////////////////////////////
/// EMC2 Lab, Clemson University, Tyler Ard								   ///
//////////////////////////////////////////////////////////////////////////////
// Mapping Function for MPC Vehicle Communication

/* Some Notes:
	// Handles mapping of vehicle communication by VehicleID
	// Map is intended to be updated after all vehicles have calculated their move
*/

# pragma once

#include <map>
#include <unordered_map>

// Define Entry Point to Check if Comms are Available from Other Vehicle
template<typename T>
bool CheckComms(const std::unordered_map <long, T> & Map, long OtherVehicleID)
{
	return !(Map.find(OtherVehicleID) == Map.end());
}

template<typename T>
bool CheckComms(const std::map <long, T> & Map, long OtherVehicleID)
{
	return !(Map.find(OtherVehicleID) == Map.end());
}

// Define Entry Point to Read Map
template<typename T>
T ReceiveComms(const std::unordered_map <long, T> & Map, long OtherVehicleID)
{
	return Map.at(OtherVehicleID);
}

template<typename T>
T ReceiveComms(const std::map <long, T> & Map, long OtherVehicleID)
{
	return Map.at(OtherVehicleID);
}

// Define Entry Point to Update Map
template<typename T>
void UpdateComms(std::unordered_map <long, T> & Map, long VehicleID, const T & ControlHorizon)
{
	Map[VehicleID] = ControlHorizon;
}

template<typename T>
void UpdateComms(std::map <long, T> & Map, long VehicleID, const T & ControlHorizon)
{
	Map[VehicleID] = ControlHorizon;
}

// Define Entry Point to Remove Mapping
template<typename T>
void RemoveComms(std::unordered_map <long, T> & Map, long VehicleID)
{
	if (CheckComms<T>(Map, VehicleID))
	{
		Map.erase(VehicleID);
	}
}

template<typename T>
void RemoveComms(std::map <long, T> & Map, long VehicleID)
{
	if (CheckComms<T>(Map, VehicleID))
	{
		Map.erase(VehicleID);
	}
}