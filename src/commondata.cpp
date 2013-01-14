

#ifndef COMMOM_DATA_INCLUDE
#define COMMOM_DATA_INCLUDE

#include <string>
#include <fstream>

#include "microdef.h"
#include "datastruct.h"
#include "FileInterp.h"


using std::string;

//Feild components
MyStruct Ex;
MyStruct Ey;
MyStruct Hz;
MyStruct Ux;
MyStruct Uy;
MyStruct Ue; //
MyStruct Ne; //density field component
MyStruct Pex; //Ex at previous step
MyStruct Pey; //Ey at previous step
MyStruct Pne; //Density Ne at previous step
MyStruct ppne;
MyStruct beta;

MyStruct Pez; //Ez at previous step
MyStruct Ceex,Ceey,Ceez;
MyStruct Cevx,Cevy,Cevz;
MyStruct Cehx,Cehy,Cehz;
//MyDataF Chxez,Chyez,Chzey,Chzex;
MyStruct Hx,Hy,Ez,Vex,Vez,Vey,Uz;

MyStruct Deff, Niu_i, Niu_a;
int denFormula;
//Gas Density at sea level
MyDataF GasDen = 2.7e25;

//file names need to read before initials
const string svix("interp/vi3700x.txt"); //file name for vi_x
const string sviy("interp/vi3700y.txt"); //file name for vi_y
const string svax("interp/va3700x.txt"); //file name for va_x
const string svay("interp/va3700y.txt"); //file name for va_y
const string svcx("interp/vc3700x.txt"); //file name for vc_x
const string svcy("interp/vc3700y.txt"); //file name for vc_y
const string slossx("interp/lossE3700x.txt"); //file name for loss_x
const string slossy("interp/lossE3700y.txt"); //file name for loss_y
//array count for interpolated fields
const int arrCount = 3700;

//Fields need to interpolate
FileInterp vi(arrCount, svix, sviy);
FileInterp va(arrCount, svax, svay);
FileInterp vc(arrCount, svcx, svcy);
FileInterp lossE(arrCount, slossx, slossy);

// for projects 2012 10

const int ArrayCount = 400;
const string EmDivNTxt("boltzmann/EmDivN.txt");
const string EnergyTxt("boltzmann/m_energy.txt");
const string vcTxt("boltzmann/vcDivN.txt");
const string viTxt("boltzmann/viDivN.txt");
FileInterp VcDivN(ArrayCount,EmDivNTxt,vcTxt);
FileInterp ViDivN(ArrayCount,EmDivNTxt,viTxt);
FileInterp EnergyDivN(ArrayCount,EmDivNTxt,EnergyTxt);
MyStruct Em;
MyStruct Pem;

MyDataF chxez,chyez,Cve,alpha;
//update coefficients
MyDataF chzex, chzey; //coefficients updating Hz
MyDataF cexhz, cexux; //coefficients updating Ex
MyDataF ceyhz, ceyuy; //coefficients updating Ey
MyDataF cezhx,cezhy, cezuz; //coefficients updating Ez

MyDataF De;
MyDataF Da;
MyDataF D_kasi_max;
MyDataF mu_e;
MyDataF mu_i;
MyDataF p = DEFAULT_AIR_PRESSURE;

MyDataF rei = DEFAULT_REI;

unsigned int CStep;
//FDTD DATA
MyDataF dt, dt_F, dt_M;
MyDataF half_dt;
MyDataF ds_F, ds_M;
MyDataF dx;
MyDataF dy;
MyDataF dt_me_e; //dt*e/me
MyDataF dt_me_e_2; //dt*e/me*2
MyDataF dt_ds2_2; //2*dt/ds_F/ds_F
MyDataF miu2DivE; //2*mu_e/e
MyDataF DtfDivDsfs; //dt_F/ds_F/ds_F
MyDataF eps_m_e_miu; //eps_0 / (e * (mu_e + mu_i))
MyDataF dt2; //dt*2
MyDataF ds_Pow_2; //ds_F*ds_F
MyDataF MueDivMui=DEFAULT_MIU_DIV;

unsigned m, m2;

//DOMAIN DATA
unsigned nx, nxp1, nxm1;
unsigned ny, nyp1, nym1;
unsigned nbound;
unsigned NumOfWaveLength;
//PARAMETERS OF INCIDENT WAVE
MyDataF f; //frequency
MyDataF k; //
MyDataF T; //
MyDataF E0, H0;
MyDataF Hx0, Hz0, Hy0, Ez0, Ex0, Ey0;
MyDataF Ratio_x, Ratio_y;
MyDataF lamda;
MyDataF omega;
MyDataF phi; //incidence wave inject angle on x-axis


int IsTMz;
int IsTEz;
MyDataF TotalTime; //in nane seconds
unsigned TotalTimeStep;
unsigned NumOfCellPerWaveLen;

MyDataF CurTime; //current time

/// scatter  bound width
unsigned scatwidth;

//plot step
unsigned PlotStep;
//file to store density information
std::ofstream denfile;

//file to store density information
std::ofstream deff_file;

unsigned Deff_Store_Index_x[10], Deff_Store_Index_y[10];
unsigned minSI, minSJ, maxSI, maxSJ;
unsigned midi, midj, pci;

int ifWithDensity = IF_WITH_DENSITY;// if with density or not

#endif
