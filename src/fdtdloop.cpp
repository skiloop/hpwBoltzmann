
#include <iostream>
#define SD(l,i,j) l.data[(i)*l.ny+j]
#define DD(s,p,i,j) (s.data[(i)*s.ny+j]+p.data[(i)*p.ny+j])
#include <cstdio>
#include <cstdlib>
#include <time.h>

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

    //unsigned int CurTimeStep; //current time step
    unsigned int xpos, ypos, sxpos, sypos;
    //FILE *fne;
    //unsigned int i,j;
    //MyDataF *CapEF;
    //MyDataF cdtF;
    //int MTimeStep=1,MultiSize;
    //int cnt=0;
    int step_per_half_ns8 = (int) (0.5 + 0.125e-9 / dt_F);
    clock_t t_start, t_end;
    //MyDataF RealEz,rhtb;
    sxpos = (int) (0.5 + ny / 2); //x position of sources
    sypos = (int) (0.5 + ny / 2);
    //fne=fopen("cnep.dat","w");
    //y position of sources

    xpos = (int) (0.5 + nx / 2 + lamda / dx); //x position of field to be captured
    ypos = (int) (0.5 + ny / 2 + lamda / dy); //y position of field to be captured

    MultiSize = (int) (0.5 + dt_F / dt);

    printf("\nm\t=\t%u\nMultiSize\t=\t%u", m, MultiSize);
    printf("\nNumber of time steps: %u\n", TotalTimeStep);
    printf("Density Time Step: %u\n", Density_Time_Step * m);
    printf("Step per half ns: %u\n", step_per_half_ns8);
    printf("\n*************************************************************\n");
    //system("pause");

    sxpos = (int) (0.5 + nx / 2 + 0.125 * lamda / dx); //tpis+3;//tpis-SCATTER_FIELD_DOMAIN_BND_SIZE/2;//x position of sources
    sypos = (int) (0.5 + (tpjs + tpje) / 2); //(ny/2); //y position of sources

    CurTime = -half_dt;

    InitMatlabEngine();
    InitConnectingInterface(phi);
    t_start = clock();
    for (Step = 1; TotalTimeStep >= Step; Step++) {

        //MyDataF ezl,ezlr;
        CurTime += half_dt;

        //E Field
        UpdateMField();
        //Hz.PlotArrays(ep);
        ApplyConnectingM(CurTime);
        UpdMagFldForPML_TMz(Hz, Ex, Ey);

        //U Field
        UpdateUField();
        //Erms
        CurTime += half_dt;
        UpdateEField();
        ApplyConnectingE(CurTime);
        UpdEltFldForPML_TMz(Ex, Ey, Hz);
        UpdateUeField();

        if (Step % PlotStep == 0)
            MatlabSimulation();
        if (Step % CStep == 0)
            CapFields(Step / CStep);
        UpdateDensity();

        cout << "Step = " << Step << '\t' << CurTime / 1e-9 << " ns" << endl;
    }//END FOR
    //SaveD(cnt);
    EndSimulation();
    t_end = clock();
    printf("Total time used : %ld\n", t_end - t_start);
    //system("pause");

}//END OF FDTDLOOP

#undef SD
