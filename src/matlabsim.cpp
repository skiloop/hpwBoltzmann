
#include <iostream>
#include "matlabsim.h"
#include "commondata.h"

using namespace std;

void MatlabSimulation() {

    //Ex.PlotArrays();
    //Ey.PlotArrays();
    Hz.PlotArrays();
    Ex.PlotArrays();
    Ey.PlotArrays();
    if(ifWithDensity) {
	Em.PlotArrays();
        Ue.PlotArrays();
        Ne.PlotArrays();
    }

}

void EndSimulation() {

    //Ex.ClearSim();
    //Ey.ClearSim();
    Ex.ClearSim();
    Ey.ClearSim();
    Hz.ClearSim();
    if(ifWithDensity) {
        Ne.ClearSim();
        Em.ClearSim();
        Ue.ClearSim();
    }
}

void InitMatlabEngine() {
    MyStruct::InitMatlabEngine();

    //Ex.InitPlot("Ex");
    //Ey.InitPlot("Ey");
    Ex.InitPlot();
    Ey.InitPlot();
    Hz.InitPlot();
    if(ifWithDensity) {
        Ne.InitPlot();
        Em.InitPlot();
        Ue.InitPlot();
    }
}

