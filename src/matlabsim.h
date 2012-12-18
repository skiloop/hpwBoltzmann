
#ifndef MATLAB_SIM_H
#define MATLAB_SIM_H
#include "microdef.h"

#ifdef MATLAB_SIMULATION
#include "engine.h"
#endif//MATLAB_SIMULATION

void InitMatlabEngine();
void MatlabSimulation();
void EndSimulation();


#endif//MATLAB_SIM_H
