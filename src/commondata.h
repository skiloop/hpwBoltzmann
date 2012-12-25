

#ifndef COMMOM_DATA_INCLUDE
#define COMMOM_DATA_INCLUDE

#include<fstream>

#include "datastruct.h"
#include "microdef.h"
#include "FileInterp.h"
using namespace std;
//Feild components
extern MyStruct Ex;
extern MyStruct Ey;
extern MyStruct Hz;
extern MyStruct Ux;
extern MyStruct Uy;
extern MyStruct Ue; //
extern MyStruct Ne; //density field component
extern MyStruct Pex; //Ey at previous step
extern MyStruct Pey; //Ey at previous step
extern MyStruct Pne; //Density Ne at previous step
extern MyStruct Pux;
extern MyStruct Puy;
extern MyStruct ppne,Deff,Niu_i,Niu_a;
extern int denFormula;
//update coefficients
extern MyDataF chzex, chzey; //coefficients updating Hz
extern MyDataF cexhz, cexux; //coefficients updating Ex
extern MyDataF ceyhz, ceyuy; //coefficients updating Ey

//Fields need to interpolate
extern FileInterp vi;
extern FileInterp va;
extern FileInterp vc;
extern FileInterp lossE;


const MyDataF N_air = 2.44e25; // in m^(-3)
extern FileInterp VcDivN;
extern FileInterp ViDivN;
extern FileInterp EnergyDivN;
extern MyStruct Em;


extern MyDataF De;
extern MyDataF Da;
extern MyDataF D_kasi_max;
extern MyDataF mu_e;
extern MyDataF mu_i;
extern MyDataF MueDivMui;
extern MyDataF p;

//extern const unsigned Density_Time_Step;
const MyDataF e = 1.602e-19;

//speed of light
const MyDataF c = 2.998E8;

extern MyDataF GasDen;
//electron mass
const MyDataF me = 9.110e-31;
const MyDataF mu_0 = 1.257e-6;
const MyDataF eps_0 = 8.854e-12;
const MyDataF vm = 760 * 5.3e9;
const MyDataF length = 1.0;
const MyDataF CFL_factor = CFL_FACTOR;
const MyDataF CourantFactor = COURANT_FACTOR;
const unsigned Density_Time_Step = DEN_TIME_STEPS;


extern MyDataF rei;


//FDTD DATA
extern MyDataF dt, dt_F, dt_M;
extern MyDataF half_dt;
extern MyDataF dt2; //dt*2
extern MyDataF ds_F, ds_M;
extern MyDataF dx;
extern MyDataF dy;
extern MyDataF ds_Pow_2;//ds_F*ds_F
extern unsigned m, m2;
extern MyDataF dt_me_e;
extern MyDataF dt_me_e_2;
extern MyDataF miu2DivE; //2*miu_e/e
extern MyDataF dt_ds2_2; //2*dt/ds_F/ds_F
extern MyDataF eps_m_e_miu;//eps_0 / (e * (mu_e + mu_i))
extern MyDataF DtfDivDsfs;
//DOMAIN DATA
extern unsigned nx, nxp1, nxm1;
extern unsigned ny, nyp1, nym1;
extern unsigned nbound;
extern unsigned NumOfWaveLength;
//PARAMETERS OF INCIDENT WAVE
extern MyDataF f; //frequency
extern MyDataF k; //
extern MyDataF T; //
extern MyDataF E0, H0;
extern MyDataF Hx0, Hz0, Hy0, Ez0, Ex0, Ey0;
extern MyDataF Ratio_x, Ratio_y;
extern MyDataF lamda;
extern MyDataF omega;
extern MyDataF phi; //incidence wave inject angle on x-axis
extern MyDataF TotalTime;
extern int IsTMz;
extern int IsTEz;

extern unsigned TotalTimeStep;
extern unsigned int CStep; //Capture steps
extern unsigned scatwidth;
extern unsigned PlotStep;
extern unsigned NumOfCellPerWaveLen;
extern MyDataF CurTime; //current time

extern std::ofstream denfile;
extern std::ofstream deff_file;

extern unsigned Deff_Store_Index_x[10],Deff_Store_Index_y[10];
extern unsigned minSI,minSJ,maxSI,maxSJ;
extern unsigned midi,midj,pci;
extern int ifWithDensity; // if with density or not
#endif

