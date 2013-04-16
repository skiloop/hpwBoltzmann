
#include <iostream>
#include "matlabsim.h"
#include "commondata.h"

using namespace std;

void MatlabSimulation() {

    //Ex.PlotArrays();
    //Ey.PlotArrays();
	if(IsTMz){
		Hz.PlotArrays();
		Ex.PlotArrays();
		Ey.PlotArrays();
	}
	if(IsTEz){
		Ez.PlotArrays();
		Hx.PlotArrays();
		Hy.PlotArrays();
	}
    if(ifWithDensity) {
	Em.PlotArrays();
        Ue.PlotArrays();
        Ne.PlotArrays();
    }

}

void EndSimulation() {

    //Ex.ClearSim();
    //Ey.ClearSim();
	if(IsTMz){
		Hz.ClearSim();
		Ex.ClearSim();
		Ey.ClearSim();
	}
	if(IsTEz){
		Ez.ClearSim();
		Hx.ClearSim();
		Hy.ClearSim();
	}

    if(ifWithDensity) {
        Ne.ClearSim();
        Em.ClearSim();
        Ue.ClearSim();
    }
}

void InitMatlabEngine() {
    MyStruct::InitMatlabEngine();

 	if(IsTMz){
		Hz.InitPlot();
		Ex.InitPlot();
		Ey.InitPlot();
	}
	if(IsTEz){
		Ez.InitPlot();
		Hx.InitPlot();
		Hy.InitPlot();
	}

    if(ifWithDensity) {
        Ne.InitPlot();
        Em.InitPlot();
        Ue.InitPlot();
    }
}

