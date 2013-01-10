
#include <iostream>
#define SD(l,i,j) l.data[(i)*l.ny+j]
#define DD(s,p,i,j) (s.data[(i)*s.ny+j]+p.data[(i)*p.ny+j])
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <iomanip>

#include "fdtdloop.h"
#include "datastruct.h"
#include "updatefields.h"
#include "commondata.h"
#include "bounddata.h"
#include "pmlbound.h"
#include "connectingdef.h"
#include "matlabsim.h"

#ifdef MATLAB_SIMULATION
extern Engine *ep;
#endif

extern unsigned tpis, tpjs, tpie, tpje;
using namespace std;
unsigned int Step = 0;
unsigned int MTimeStep = 1, MultiSize;

void printsize(int xpos, int ypos, int sxpos, int sypos);

void fdtdloop() {

    unsigned int CurTimeStep; //current time step

    //FILE *fne;
    //unsigned int i,j;
    //MyDataF *CapEF;
    MyDataF cdtF;
    //int MTimeStep=1,MultiSize;
    //int cnt=0;
    int step_per_half_ns8 = (int) (0.5 + 0.125e-9 / dt_F);
    clock_t t_start, t_end;

    //MyDataF RealEz,rhtb;
    /*
        unsigned int xpos, ypos, sxpos, sypos;
        sxpos = (int) (0.5 + ny / 2); //x position of sources
        sypos = (int) (0.5 + ny / 2);
        xpos = (int) (0.5 + nx / 2 + lamda / dx); //x position of field to be captured
        ypos = (int) (0.5 + ny / 2 + lamda / dy); //y position of field to be captured

        sxpos = (int) (0.5 + nx / 2 + 0.125 * lamda / dx); //tpis+3;//tpis-SCATTER_FIELD_DOMAIN_BND_SIZE/2;//x position of sources
        sypos = (int) (0.5 + (tpjs + tpje) / 2); //(ny/2); //y position of sources
    */
    MultiSize = (int) (0.5 + dt_F / dt);

    cout<<endl;
    cout<<setw(20)<<"m:"<<m<<endl;
    cout<<setw(20)<<"MultiSize:"<< MultiSize<<endl;
    cout<<setw(20)<<"Total time steps:"<<TotalTimeStep<<endl;
    cout<<setw(20)<<"Density Time Steps:"<<Density_Time_Step * m<<endl;
    cout<<setw(20)<<"Step per half ns:"<< step_per_half_ns8<<endl;
    cout<<endl;
    cout<<"============================================================="<<endl;
    //system("pause");

    CurTime = -half_dt;

    InitMatlabEngine();
    InitConnectingInterface(phi);
    t_start = clock();
    for(CurTimeStep=1; CurTimeStep<=Density_Time_Step&&TotalTimeStep>=Step; CurTimeStep++)
    {
        for(cdtF = 0,MTimeStep = 1; MTimeStep<=MultiSize&&cdtF < dt_F&&TotalTimeStep>=Step; cdtF=cdtF+dt,MTimeStep++,Step++)
        {
            //MyDataF ezl,ezlr;
            CurTime += half_dt;
            //E Field
            UpdateMField();
            //Hz.PlotArrays(ep);
            ApplyConnectingM(CurTime);
            UpdMagFldForPML_TMz(Hz, Ex, Ey);
            //U Field
            if(ifWithDensity)UpdateUField();
            //Erms
            CurTime += half_dt;
            UpdateEField();
            ApplyConnectingE(CurTime);
            UpdEltFldForPML_TMz(Ex, Ey, Hz);

            if(ifWithDensity) {
                if(denFormula==4)CalEmax();
                else UpdateUeField();
            }

            if (Step % PlotStep == 0)
                MatlabSimulation();
            if (Step % CStep == 0)
                CapFields(Step / CStep);
        }
        if(ifWithDensity){
			InterpEmax();
			Pem=Em;
			Em.CaptData(CStep);
		}
        if(ifWithDensity) UpdateDensity();
        MatlabSimulation();

        cout << Step << '\t' << CurTime / 1e-9 << endl;
    }//END FOR
    //SaveD(cnt);

    EndSimulation();
    t_end = clock();
    cout<<"========================================"<<endl;
    cout<<"Total time used : "<< t_end - t_start<<endl;
    //system("pause");

}//END OF FDTDLOOP

#undef SD
