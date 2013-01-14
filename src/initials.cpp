

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "commondata.h"
#include "microdef.h"
#include "initials.h"

using namespace std;

void InitStoreIndex() {
    ifstream cfile("cdeff.txt");
    if(!cfile.is_open()) {
        cerr<<"Cannot open file cdeff.txt!"<<endl;
        exit(-1);
    }
    try {
        MyDataF cx,cy;
        unsigned c = 0;
        cout<<"Read Capture files..."<<endl;
        minSI = Ne.nx;
        minSJ = Ne.ny;
        maxSI = 0;
        maxSJ = 0;

        while(cfile.good()&&c<10)
        {
            cfile>>cx>>cy;
            Deff_Store_Index_x[c] = (cx-1.5)*lamda/ds_F+midi;
            Deff_Store_Index_y[c] = (cy-1.5)*lamda/ds_F+midj;
            cout<<c<<'\t'<<Deff_Store_Index_x[c]<<'\t'<<Deff_Store_Index_y[c]<<endl;

            //index range for computation area
            if(minSI>Deff_Store_Index_x[c])minSI=Deff_Store_Index_x[c];
            if(minSJ>Deff_Store_Index_y[c])minSJ=Deff_Store_Index_y[c];
            if(maxSI<Deff_Store_Index_x[c])maxSI=Deff_Store_Index_x[c];
            if(maxSJ<Deff_Store_Index_y[c])maxSJ=Deff_Store_Index_y[c];

            c++;
        }
        cout<<"MinSI\tMinSJ\tMaxSI\t\tMaxSJ\t"<<endl;
        cout<<minSI<<'\t'<<minSJ<<'\t'<<maxSI<<'\t'<<maxSJ<<endl;
        cfile.close();
    } catch(exception &e) {
        cerr<<"Caught exception "<<e.what()<<endl;
        cfile.close();
        exit(-1);
    }
}
void OpenFiles()
{
    denfile.open("resdata/den.txt");
    if (!denfile.is_open()) {
        cerr << "open denfile failed!" << endl;
        exit(-1);
    } else cout<<"denfile is opened"<<endl;
    deff_file.open("resdata/deff.txt");
    if (!deff_file.is_open()) {
        cerr << "open denfile failed!" << endl;
        exit(-1);
    } else cout<<"deff file is opened"<<endl;
}
/*
 * @brief Initial Problem Domain
 *
 */
void InitDomain() {

    int nxi, nyi;
    //int nxip1, nyip1;

    nxp1 = nx = nxm1 = ny = nyp1 = nym1 = 0;

    nxi = (unsigned int) (0.5 + (unsigned int) (0.5 + lamda / dx) * NumOfWaveLength);
    nyi = (unsigned int) (0.5 + (unsigned int) (0.5 + lamda / dy) * NumOfWaveLength);

    //    nxip1 = nxi + 1;
    //    nyip1 = nyi + 1;

    nx = nxi + 2 * (nbound + scatwidth);
    ny = nyi + 2 * (nbound + scatwidth);

    nxp1 = nx + 1;
    nyp1 = ny + 1;
    nxm1 = nx - 1;
    nym1 = ny - 1;
}

/*
 * @brief Initial FDTD problem
 *
 */
void InitFDTDProblem(const string &fname) {

    if (fname != "") {
        string line;

        ifstream infile(fname.c_str());
        int width = 40;
        //int vwidth = 15;
        if (!infile.is_open()) {
            cerr << "file open failed!\t" << fname << endl;
            exit(-1);
        }
        try {
            cout << "========================<Read Paremeters From File>===========================" << endl;
            // read frequency
            infile >> f;
            cout << setw(width) << "incident wave frequency:" << f << endl;
            // read amplitude of electric component
            infile >> E0;
            cout << setw(width) << "electric component amplitude:" << E0 << endl;
            // read phrase angle
            infile >> phi;
            cout << setw(width) << "incident phrase angle:" << phi << endl;
            infile >> TotalTime;
            cout << setw(width) << "Total Propagation Time:" << TotalTime << endl;
            infile >> NumOfWaveLength;
            cout << setw(width) << "Number of wavelength in domain:" << NumOfWaveLength << endl;
            infile >> NumOfCellPerWaveLen;
            cout << setw(width) << "Number of coarse cells per wavelength:" << NumOfCellPerWaveLen << endl;
            infile >> m; //number of fine cells per coarse cell
            cout << setw(width) << "Number of fine cells per coarse cell:" << m << endl;
            infile >> rei; //number of fine cells per coarse cell
            cout << setw(width) << "rei:" << rei << endl;
            infile >> nbound; //number of coarse cells in PML boundary
            cout << setw(width) << "PML width:" << nbound << endl;
            infile >> scatwidth; // number of coarse cells between connecting interface and PML boundary
            cout << setw(width) << "scatter field width:" << scatwidth << endl;
            infile >> PlotStep;
            cout << setw(width) << "PlotStep:" << PlotStep << endl;
            infile >> CStep;
            cout << setw(width) << "Capture Step:" << CStep << endl;
            infile >> denFormula;
            cout << setw(width) << "Density Formula:" << denFormula << endl;
            infile >> ifWithDensity;
            cout << setw(width) << "if with Density:"<<ifWithDensity <<endl;
            infile.close();
        } catch (exception &e) {
            cerr << e.what() << endl;
            exit(-1);
        }
        if (f < 10)
            //initial wave parameters
            f = FREQUENCY;
        lamda = c / f;
        T = 1 / f;
        omega = 2 * M_PI*f;
        k = omega / c;

        if (E0 == 0)
            E0 = E_0;
        H0 = E0 * sqrt(eps_0 / mu_0);
        //wave incident angle
        phi = phi * M_PI;
        // wave amplitudes
        if (IsTMz) {
            Ratio_x = sin(phi); //sin(0.5*M_PI);
            Ratio_y = -cos(phi); //cos(0.5*M_PI);

            Ex0 = Ratio_x*E0;
            Hz0 = H0;
            Ey0 = Ratio_y*E0;
        }
        if (IsTEz) {
            Ratio_x = sin(phi); //sin(0.5*M_PI);
            Ratio_y = cos(phi); //cos(0.5*M_PI);

            Hx0 = Ratio_x*H0;
            Hy0 = Ratio_y*H0;
            Ez0 = E0;
        }

        //wave type
        IsTMz = _SOURCE_TMZ_;
        IsTEz = _SOURCE_TEZ_;

        //number of wavelength in domain
        if (NumOfWaveLength > 1000)
            NumOfWaveLength = NUMBER_OF_WAVELENGTHS_IN_DOMAIN;

        //initial FDTD grids
        if (NumOfCellPerWaveLen <= 0)
            NumOfCellPerWaveLen = NUMCELLPERWAVELEN;
        ds_M = dx = dy = lamda / NumOfCellPerWaveLen; //

        //initial  FDTD time step
        dt = 0.5 * dx / c;

        //initial total iteration step
        if (TotalTime < 0)
            TotalTime = TOTAL_TIME;
        TotalTimeStep = (int) (0.5 + TotalTime / dt);

        //initial fine grid
        if (m > 10000)
            m = FINE_GRID_SIZE;

        //set PML boundary width
        if (nbound < 4 || nbound > 30)
            nbound = NUMBER_OF_CELLS_IN_PML_BOUND;

        //set free space width
        if (scatwidth < 2 || scatwidth > 20)
            scatwidth = 4; //in Yee cells
        // plot step
        if (PlotStep < 1 || PlotStep > TotalTimeStep)
            PlotStep = 5;
        //how often to store data(in time step)
        if(CStep == 0) {
            CStep = (unsigned int)(TotalTimeStep/TotalTime*2.5e-10+0.5);
        }
        if (CStep > 10000)
            CStep = 100;

    } else {

        //initial wave parameters
        f = FREQUENCY;
        lamda = c / f;
        E0 = E_0;
        H0 = E0 * sqrt(eps_0 / mu_0);
        T = 1 / f;
        omega = 2 * M_PI*f;
        k = omega / c;

        //wave incident angle
        phi = INC_ANGLE ;

        //wave type
        IsTMz = _SOURCE_TMZ_;
        IsTEz = _SOURCE_TEZ_;

        //number of wavelength in domain
        NumOfWaveLength = NUMBER_OF_WAVELENGTHS_IN_DOMAIN;

        //initial FDTD grids
        NumOfCellPerWaveLen = NUMCELLPERWAVELEN;
        ds_M = dx = dy = lamda/NumOfCellPerWaveLen;
        //initial total iteration step
        TotalTimeStep = (int) (0.5 + TOTAL_TIME / dt);

        //initial fine grid
        m = FINE_GRID_SIZE;

        //set PML boundary width
        nbound = NUMBER_OF_CELLS_IN_PML_BOUND;

        // wave amplitudes
        if (IsTMz) {
            Ratio_x = sin(phi); //sin(0.5*M_PI);
            Ratio_y = -cos(phi); //cos(0.5*M_PI);

            Ex0 = Ratio_x*E0;
            Hz0 = H0;
            Ey0 = Ratio_y*E0;
        }
        if (IsTEz) {
            Ratio_x = sin(phi); //sin(0.5*M_PI);
            Ratio_y = cos(phi); //cos(0.5*M_PI);

            Hx0 = Ratio_x*H0;
            Hy0 = Ratio_y*H0;
            Ez0 = E0;
        }

        //initial  FDTD time step
        dt = 0.5 * dx / c;

        //set free space width
        scatwidth = 4; //in Yee cells
        // plot step
        PlotStep = 5;
        //how often to store data(in time step)
        CStep = 100;

        //Density Formula
        denFormula = 4;
    }
}

/*
 * @brief initial electricity density distribution
 *
 */
void InitEleDen() {


    MyDataF bndsz = m*(nbound + scatwidth+0.5);
    MyDataF temp = 2500e-12;
    unsigned int i, j;
    unsigned int mi=0,mj=0;
    MyDataF maxne=0;
    MyDataF cx = 2.25*lamda+bndsz * ds_F;
    MyDataF cy = 0.5*NumOfWaveLength*lamda+bndsz * ds_F;
    MyDataF px,py;

    for (i = 0; i < Ne.nx; i++) {
        px = (i * ds_F - cx);
        px = px*px;
        for (j = 0; j < Ne.ny; j++) {
            py = (j * ds_F - cy);
            py = py*py;
            Ne.data[i][j] = NE0 * exp(-(px+py) / temp);
            if(maxne<Ne.data[i][j])
            {
                mi = i;
                mj = j;
                maxne=Ne.data[i][j];
            }
        }
    }
    cout<<"max density position:("<<mi<<','<<mj<<')'<<endl;
}

void InitCee(){
	if (ifWithDensity)
	{
	
	unsigned int i,j,im,jm;	
	if(IsTMz){		
		for(i=0,im=m/2;i<Ceex.nx;++i,im+=m)
			for(j=0,jm=0;j<Ceex.ny;j++,jm+=m){
				Ceex.data[i][j]=-beta.data[im][jm]/(1+beta.data[im][jm])+1/(1+beta.data[im][jm]);
			}
			for(i=0,im=0;i<Ceey.ny;++i,im+=m)
				for(j=0,jm=m2;j<Ceey.ny;j++,jm+=m){
					Ceey.data[i][j]=-beta.data[im][jm]/(1+beta.data[im][jm])+1/(1+beta.data[im][jm]);
				}
	}
	if(IsTEz){
		for(i=0,im=0;i<Ceez.nx;++i,im+=m)
			for(j=0,jm=0;j<Ceez.ny;j++,jm+=m){
				Ceez.data[i][j]=-beta.data[im][jm]/(1+beta.data[im][jm])+1/(1+beta.data[im][jm]);
			}
	}
	//PrintData(Ceez);
	}else{
		if (IsTEz)
		{
			Ceez.InitStructData(0.0);
		}
		if (IsTMz)
		{
			Ceey.InitStructData(0.0);
			Ceez.InitStructData(0.0);
		}
	}
}
void InitCeh(){
	if (ifWithDensity)
	{
	
	unsigned int i,j,index=0;
	unsigned int im,jm;
	//MyDataF temp=dt/eps_0;
	if(IsTMz){		
		for(i=0,im=m2;i<Cehx.nx;++i,im+=m)
			for(j=0,jm=0;j<Cehx.ny;j++,jm+=m){
				Cehx.data[i][j]=1/(eps_0/dt+eps_0*beta.data[im][jm]/dt);
			}
			for(i=0,im=0;i<Cehy.ny;++i,im+=m)
				for(j=0,jm=m2;j<Cehy.ny;j++,jm+=m){
					Cehy.data[i][j]=1/(eps_0/dt+eps_0*beta.data[im][jm]/dt);
				}
	}
	if(IsTEz){
		for(i=0,im=0;i<Cehz.nx;++i,im+=m)
			for(j=0,jm=0;j<Cehz.ny;j++,jm+=m){
				Cehz.data[i][j]=1/(eps_0/dt+eps_0*beta.data[im][jm]/dt);
			}
	}
	//PrintData(Cehz);
	}else{
		if (IsTEz)
		{
			Cehz.InitStructData(0.0);
		}
		if (IsTMz)
		{
			Cehx.InitStructData(0.0);
			Cehy.InitStructData(0.0);
		}
	}
}
void InitCev(MyDataF alpha_ev){
	if (ifWithDensity)
	{
	
	unsigned int i,j,index=0;
	MyDataF temp;
	unsigned int im,jm;

	temp = 0.5*e*dt*(1+alpha_ev)/eps_0;
	if(IsTMz){		
		for(i=0,im=m2;i<Cevx.nx;++i,im+=m)
			for(j=0,jm=0;j<Cevx.ny;j++,jm+=m){
				Cevx.data[i][j]=temp*Ne.data[im][jm]/(1+beta.data[im][jm]);
			}
			for(i=0,im=0;i<Cevy.ny;++i,im+=m)
				for(j=0,jm=m2;j<Cevy.ny;j++,jm+=m){
					Cevy.data[i][j]=temp*Ne.data[im][jm]/(1+beta.data[im][jm]);
				}
	}
	if(IsTEz){
		for(im=0,i=0;i<Cevz.nx;++i,im+=m)
			for(j=0,jm=0;j<Cevz.ny;j++,jm+=m){
				Cevz.data[i][j]=temp*Ne.data[im][jm]/(1+beta.data[im][jm]);
			}
			//PrintData(Cevz);
	}
	}else{
		if (IsTEz)
		{
			Cevz.InitStructData(0.0);
		}
		if (IsTMz)
		{
			Cevy.InitStructData(0.0);
			Cevy.InitStructData(0.0);
		}
	}
}
void InitBeta(MyDataF gamma){
	if(ifWithDensity){
		unsigned int i,j;
		MyDataF temp=0.25*e*e*dt*dt/me/eps_0/gamma;
		for(i=0;i<beta.nx;i++)
			for(j=0;j<beta.ny;j++)
				beta.data[i][j]=temp*Ne.data[i][j];
	}		
}
void InitCoeffForVec(){

	MyDataF gamma,a=0;
	a	  = 0.5*dt*vm;
	gamma = 1+a;
	alpha = (1-a)/(1+a);
	Cve   = e*dt*0.5/(me*gamma);

	InitBeta(gamma);
	/*PrintData(beta);*/
	InitCee();
	InitCeh();
	InitCev(alpha);

}
void InitCoeff() {

	if(IsTMz){
		cexux = -dt * e / eps_0;
		cexhz = dt / dy / eps_0;

		ceyhz = -dt / eps_0 / dx;
		ceyuy = -e * dt / eps_0;

		chzex = dt / mu_0 / dy;
		chzey = -dt / mu_0 / dx;
	}
	if(IsTEz){
		
		chxez = -dt/mu_0/dy;
		chyez = dt/mu_0/dx;

		cezuz = -dt*e/eps_0;
		cezhx = -dt/dy/eps_0;
		cezhy = dt/dx/eps_0;

	}
	if (denFormula==4)
	{
		InitCoeffForVec();
	}
}

void CreateFields() {
    string path = "resdata/";
    if (IsTMz) {

        Ex.CreateStruct(nx, nyp1);
        Ey.CreateStruct(nxp1, ny);
        Hz.CreateStruct(nx, ny);

        Pex.CreateStruct(Ex);
        Pey.CreateStruct(Ey);

        Ux.CreateStruct(Ex,0);
        Uy.CreateStruct(Ey,0);

        Ex.SetName(path + "ex");
        Ey.SetName(path + "ey");
        Pex.SetName(path + "pex");
        Pey.SetName(path + "pey");
        Hz.SetName(path + "hz");

        Ux.SetName(path + "ux");
        Uy.SetName(path + "uy");
    }
    if (IsTEz) {
        Hx.CreateStruct(nxp1, ny);
        Hy.CreateStruct(nx, nyp1);
        Ez.CreateStruct(nxp1, nyp1);
		Pez.CreateStruct(Ez);
		Uz.CreateStruct(Ez,0.0);

        Hx.SetName(path + "hx");
        Hy.SetName(path + "hy");
        Pez.SetName(path + "pez");
        Ez.SetName(path + "ez");
		Uz.SetName(path + "uz");
    }
    if(ifWithDensity) {
        Ne.CreateStruct(nxp1*m, nyp1 * m);
        Pne.CreateStruct(Ne);
        Ue.CreateStruct(Ne);
        if(denFormula==2) {
            ppne.CreateStruct(Ne);
            ppne.SetName("ppne");
        }
        if(denFormula==3) {
            Deff.CreateStruct(Ne);
            Niu_i.CreateStruct(Ne);
            Niu_a.CreateStruct(Ne);
            Niu_i.SetName("Niu_i");
            Deff.SetName("deff");
            Niu_a.SetName("Niu_a");
        }
        Ne.SetName(path + "ne");
        Pne.SetName(path + "pne");
        Ue.SetName(path + "ue");
    	Em.CreateStruct(Ne);
    	Pem.CreateStruct(Em);
    	Em.SetName("Emax");
        Em.SetName(path + "Emax");
    }
}

void InitBreakDownParam() {
    //////////////////////////////////////////////
    ///END OF SET COMMON DATA ZERO
    //////////////////////////////////////////////

    mu_e = e / me / vm; //3.7e-2;
    mu_i = mu_e / MueDivMui; //mu_e/mu_i (MueDivMui)ranges from 100 to 200 ;
    miu2DivE = mu_e * 2 / e;
    De = mu_e*2*1.6021e-19/e;//8.73e-2; //
    Da = De * mu_i / mu_e;

    if (rei < 0)rei = 0;
    D_kasi_max = (Da > De ? Da : De);
}
//initial data that frequently used

void InitComData() {

    m2 = (int) floor(0.5 + m / 2.0);
    ds_F = ds_M / m;
    dt_F = T;//CFL_factor*ds_F*ds_F*0.5/D_kasi_max;// T; //
    dt_M = dt;
    half_dt = dt / 2;
    dt2 = 2 * dt;
    dt_me_e = dt * e / me;
    dt_me_e_2 = dt_me_e * 2;

    ds_Pow_2 = ds_F*ds_F;
    dt_ds2_2 = dt_F*2 / ds_F / ds_F; //2*dt_F/ds_F/ds_F
    DtfDivDsfs = dt_F/ds_F/ds_F; //dt_F/ds_F/ds_F
    eps_m_e_miu = eps_0 / (e * (mu_e + mu_i));
    CurTime = 0.0;
    pci  = (2.25-1.5)*lamda/ds_F+midi;
}

void PrintParam() {
    int width = 10;
    //Feild components
    cout << "(pci,midj)=("<<pci<<','<<midj<<')'<<endl;
    cout << "=====================< component size >===================" << endl;
    cout << "Ex:" << Ex.nx << '\t' << Ex.ny << endl;
    cout << "Ey:" << Ey.nx << '\t' << Ey.ny << endl;
    cout << "Ux:" << Ux.nx << '\t' << Ux.ny << endl;
    cout << "Uy:" << Uy.nx << '\t' << Uy.ny << endl;
    cout << "Hz:" << Hz.nx << '\t' << Hz.ny << endl;
    cout << "Ue:" << Ue.nx << '\t' << Ue.ny << endl;
    cout << "Ne:" << Ne.nx << '\t' << Ne.ny << endl;
    cout << "Pne:" << Pne.nx << '\t' << Pne.ny << endl;
    cout << "Pex:" << Pex.nx << '\t' << Pex.ny << endl;
    cout << "Pey:" << Pey.nx << '\t' << Pey.ny << endl;

    //update coefficients
    cout << "===========<Coefficients>====================" << endl;
    cout << "CHz:" << chzex << '\t' << chzey << endl; //coefficients updating Hz
    cout << "CEx:" << cexhz << '\t' << cexux << endl; //coefficients updating Ex
    cout << "CEy:" << ceyhz << '\t' << ceyuy << endl; //coefficients updating Ey

    //Fields need to interpolate
    cout << "=============<Interpolate size>================" << endl;
    width = 20;
    cout << setw(width) << "length of interp Niu_i:" << vi.n << endl;
    cout << setw(width) << "length of interp Niu_a:" << va.n << endl;
    cout << setw(width) << "length of interp Niu_c:" << vc.n << endl;
    cout << setw(width) << "length of interp lossE:" << lossE.n << endl;

    cout << "==========<Breakdown parameters>=================" << endl;
    cout << setw(width) << "De:" << De << endl;
    cout << setw(width) << "Da:" << Da << endl;
    cout << setw(width) << "D_kasi_max:" << D_kasi_max << endl;
    cout << setw(width) << "Mu_e/Mu_i:" << MueDivMui <<endl;
    cout << setw(width) << "Miu_e:" << mu_e << endl;
    cout << setw(width) << "Miu_i:" << mu_i << endl;
    cout << setw(width) << "R_ei:" << rei << endl;

    cout << "===================<Common constants>=================" << endl;
    width = 20;
    cout << setw(width) << "p:" << p << endl;
    //cout<< const unsigned Density_Time_Step<<endl;
    cout << setw(width) << "e:" << e << endl;
    //speed of light
    cout << setw(width) << "c:" << c << endl;
    cout << setw(width) << "GasDen:" << GasDen << endl;
    cout << setw(width) <<"N_air:"<<N_air<<endl;
    //electron mass
    cout << setw(width) << "Mass_e:" << me << endl;
    cout << setw(width) << "Miu_0:" << mu_0 << endl;
    cout << setw(width) << "Eps_0:" << eps_0 << endl;
    cout << setw(width) << "Niu_m:" << vm << endl;
    cout << setw(width) << "denFormula:" << denFormula << endl;


    //FDTD DATA
    cout << "==========<Frequently used variables>=========" << endl;
    cout << setw(width) << "dt:" << dt << endl;
    cout << setw(width) << "dt_F:" << dt_F << endl;
    cout << setw(width) << "dt_M:" << dt_M << endl;
    cout << setw(width) << "half_dt:" << half_dt << endl;
    cout << setw(width) << "dt2:" << dt2 << endl; //dt*2
    cout << setw(width) << "dx:" << dx << endl;
    cout << setw(width) << "dy:" << dy << endl;
    cout << setw(width) << "ds_F:" << ds_F << endl;
    cout << setw(width) << "ds_M:" << ds_M << endl;
    cout << setw(width) << "m:" << m << endl;
    cout << setw(width) << "m2:" << m2 << endl;
    cout << setw(width) << "dt_me_e:" << dt_me_e << endl;
    cout << setw(width) << "dt_me_e_2:" << dt_me_e_2 << endl;
    cout << setw(width) << "dt_ds2_2:" << dt_ds2_2 << endl; //2*dt/dx/dx
    cout << setw(width) << "miu2DivE:" << miu2DivE << endl; 
    cout << setw(width) << "DtfDivDsf :" << DtfDivDsfs << endl; //2*dt/dx/dx
    //DOMAIN DATA
    cout << setw(width) << "nx:" << nx << endl;
    cout << setw(width) << "ny:" << ny << endl;
    cout << setw(width) << "pml width:" << nbound << endl;
    cout << setw(width) << "NumOfWaveLength:" << NumOfWaveLength << endl;
    cout << "==============<Plain Wave Parameters>============" << endl;
    width = 20;
    //PARAMETERS OF INCIDENT WAVE
    cout << setw(width) << "frequency:" << f << endl; //frequency
    cout << setw(width) << "k:" << k << endl; //
    cout << setw(width) << "T:" << T << endl; //
    cout << setw(width) << "lambda:" << lamda << endl;
    cout << setw(width) << "omega:" << omega << endl;
    cout << setw(width) << "incident phrase:" << phi << endl; //incidence wave inject angle on x-axis
    cout << setw(width) << "E0:" << E0 << endl;
    cout << setw(width) << "H0:" << H0 << endl;
    cout << setw(width) << "NE0:" << NE0 << endl;
    cout << setw(width) << "Ex0:" << Ex0 << endl;
    cout << setw(width) << "Ey0:" << Ey0 << endl;
    cout << setw(width) << "Ez0:" << Ez0 << endl;
    cout << setw(width) << "Hz0:" << Hz0 << endl;
    cout << setw(width) << "Hy0:" << Hy0 << endl;
    cout << setw(width) << "Hx0:" << Hx0 << endl;

    cout << setw(width) << "IsTMz:" << IsTMz << endl;
    cout << setw(width) << "IsTEz:" << IsTEz << endl;

    cout << setw(width) << "TotalTimeStep:" << TotalTimeStep << endl;
    cout << setw(width) << "save steps:" << CStep << endl; //Capture steps
    cout << setw(width) << "scatwidth:" << scatwidth << endl;
    cout << setw(width) << "plot step:" << PlotStep << endl;
    cout << setw(width) << "NumOfCellPerWaveLen:" << NumOfCellPerWaveLen << endl;
    cout << setw(width) << "CurTime:" << CurTime << endl; //current time
    cout << setw(width) << "ifWithDensity:"<<ifWithDensity<<endl;
    cout << "======================================================" << endl;
}

void Initial(const string &fname) {
    OpenFiles();
    InitFDTDProblem(fname);
    InitDomain();
    InitBreakDownParam();

    InitComData();
    CreateFields();
    if(ifWithDensity)InitEleDen();
    InitCoeff();
    midi = Ne.nx/2;
    midj = Ne.ny/2;
    InitStoreIndex();

    PrintParam();
}
