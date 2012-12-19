
#include <iostream>
#include <math.h>

#ifdef MATLAB_SIMULATION
#include "engine.h"
#endif

#define _USE_MATH_DEFINES

#include "connectingdef.h"
#include "commondata.h"
#include "bounddata.h"

using namespace std;
//delay arrays for incident fields
MyDataF *Dhyl, *Dhxb, *Dhxt, *Dhyr, *Dezl, *Dezr, *Dezt, *Dezb;
MyDataF *Deyl, *Dexb, *Dext, *Deyr, *Dhzl, *Dhzr, *Dhzt, *Dhzb;

MyDataF *Ei, *Hi;
int eilen, hilen;
MyDataF CMur;

MyDataF Vp_ratio; //phase velocity ratio
//first point of total field in x direction
unsigned tpis;
//first point of total field in y direction
unsigned tpjs;
//last point of total field in x direction
unsigned tpie;
//last point of total field in y direction
unsigned tpje;
//wave start position in x direction
int xs;
//wave start position in y direction
int ys;
//wave end position in x direction
int xe;
//wave end position in x direction
int ye;


MyDataF chiei;
MyDataF ceihi;
MyDataF cos_phi;
MyDataF sin_phi;

MyDataF Chzex;
MyDataF Chzey;
MyDataF Ceyhz;
MyDataF Cexhz;
MyDataF Chyez;
MyDataF Chxez;

MyDataF delays(MyDataF px, MyDataF py) {
    return ((px - xs) * cos_phi + (py - ys) * sin_phi);
}
#define MAX_ITER 10000

MyDataF PhaseVelRatio(MyDataF angle) {
    MyDataF A, B, CC, S, N, kp, err;
    int i = 0;
    N = lamda / dx;
    S = c * dt / dx;
    A = 0.5 * dx * cos(angle);
    B = 0.5 * dx * sin(angle);
    CC = sin(M_PI * S / N) * sin(M_PI * S / N) / S / S;
    err = 1;
    kp = 1;

    while (fabs(err) > 1e-6 && i < MAX_ITER) {

        err = (sin(A * kp) * sin(A * kp) + sin(B * kp) * sin(B * kp) - CC) / (A * sin(2 * A * kp) + B * sin(2 * B * kp));
        kp = kp - err;
        i++;
    }
    // cout << i << '\t' << kp << endl;
    return 2 * M_PI * c / kp;
}

/**
 * CalculateDelay
 * @brief calculate delay for connecting interface
 * @param IncAngle phase angle of the incident wave
 */

void CalculateDelay(MyDataF IncAngle) {
    int i; // width, height;

    //width = (tpie - tpis + 1);
    //height = (tpje - tpjs + 1);
    //double t_per_cell = dx / c;

    cos_phi = cos(IncAngle);
    sin_phi = sin(IncAngle);

    IncAngle = IncAngle - 2 * M_PI * (floor(IncAngle / (2 * M_PI)));

    if ((IncAngle <= 0.5 * M_PI)) {
        xs = tpis;
        ys = tpjs;
        xe = tpie;
        ye = tpje;
    } else if ((IncAngle <= M_PI)) {
        xs = tpie;
        ys = tpjs;
        xe = tpis;
        ye = tpje;
    } else if ((IncAngle <= 1.5 * M_PI)) {
        xs = tpie;
        ys = tpje;
        xe = tpis;
        ye = tpjs;
    } else {
        xs = tpis;
        ys = tpje;
        xe = tpie;
        ye = tpjs;
    }

    if (IsTEz) {
        int Max_Index_Dez_x;
        int Max_Index_Dez_y;

        int Max_Index_Dhx;
        int Max_Index_Dhy;

        Max_Index_Dhx = Max_Index_Dez_x = (int) fabs((double) (xe - xs));
        Max_Index_Dhy = Max_Index_Dez_y = (int) fabs((double) (ye - ys));

        Dhxb = (MyDataF*) malloc((Max_Index_Dhx + 1) * sizeof (MyDataF));
        Dhxt = (MyDataF*) malloc((Max_Index_Dhx + 1) * sizeof (MyDataF));
        if (Dhxb == NULL || Dhxt == NULL) {
            cerr << "Cannot allocate enough memory for Dhx!" << endl;
            exit(EXIT_FAILURE);
        }
        Dhyl = (MyDataF*) malloc((Max_Index_Dhy + 1) * sizeof (MyDataF));
        Dhyr = (MyDataF*) malloc((Max_Index_Dhy + 1) * sizeof (MyDataF));
        if (Dhyl == NULL || Dhyr == NULL) {
            cerr << "Cannot allocate enough memory for Dhy!" << endl;
            exit(EXIT_FAILURE);
        }

        Dezb = (MyDataF*) malloc((Max_Index_Dez_x + 1) * sizeof (MyDataF));
        Dezt = (MyDataF*) malloc((Max_Index_Dez_x + 1) * sizeof (MyDataF));
        Dezl = (MyDataF*) malloc((Max_Index_Dez_y + 1) * sizeof (MyDataF));
        Dezr = (MyDataF*) malloc((Max_Index_Dez_y + 1) * sizeof (MyDataF));

        if (Dezt == NULL || Dezb == NULL || Dezl == NULL || Dezr == NULL) {
            cerr << "Cannot allocate enough memory for Dez!" << endl;
            exit(EXIT_FAILURE);
        }

        for (i = 0; i <= Max_Index_Dhx; i++) {
            Dhxb[i] = delays(tpis + i, tpjs - 0.5);
            Dhxt[i] = delays(tpis + i, tpje + 0.5);
        }
        for (i = 0; i <= Max_Index_Dhy; i++) {
            Dhyl[i] = delays(tpis - 0.5, tpjs + i);
            Dhyr[i] = delays(tpie + 0.5, tpjs + i);
        }
        for (i = 0; i <= Max_Index_Dez_x; i++) {
            Dezb[i] = delays(tpis + i, tpjs);
            Dezt[i] = delays(tpis + i, tpje);
        }
        for (i = 0; i <= Max_Index_Dez_y; i++) {
            Dezl[i] = delays(tpis, tpjs + i);
            Dezr[i] = delays(tpie, tpjs + i);
        }
    }
    if (IsTMz) {
        int Max_Index_Dhz_x;
        int Max_Index_Dhz_y;

        int Max_Index_Dex;
        int Max_Index_Dey;

        Max_Index_Dex = Max_Index_Dhz_x = (int) fabs((double) (xe - xs));
        Max_Index_Dey = Max_Index_Dhz_y = (int) fabs((double) (ye - ys));

        Dexb = new MyDataF[Max_Index_Dex]();
        Dext = new MyDataF[Max_Index_Dex]();

        if (Dexb == NULL || Dext == NULL) {
            cerr << "Cannot allocate enough memory for Dex!" << endl;
            exit(EXIT_FAILURE);
        }
        Deyl = new MyDataF[Max_Index_Dey]();
        Deyr = new MyDataF[Max_Index_Dey]();
        if (Deyl == NULL || Deyr == NULL) {
            cerr << "Cannot allocate enough memory for Dey!" << endl;
            exit(EXIT_FAILURE);
        }

        Dhzb = new MyDataF[Max_Index_Dhz_x]();
        Dhzt = new MyDataF[Max_Index_Dhz_x]();
        Dhzl = new MyDataF[Max_Index_Dhz_y]();
        Dhzr = new MyDataF[Max_Index_Dhz_y]();

        if (Dhzt == NULL || Dhzb == NULL || Dhzl == NULL || Dhzr == NULL) {
            cerr << "Cannot allocate enough memory for Dhz!" << endl;
            exit(EXIT_FAILURE);
        }
        //E incident field is in total region
        for (i = 0; i < Max_Index_Dex; i++) {
            Dexb[i] = delays(tpis + i + 0.5, tpjs) + 1;
            Dext[i] = delays(tpis + i + 0.5, tpje) + 1;
        }
        for (i = 0; i < Max_Index_Dey; i++) {
            Deyl[i] = delays(tpis, tpjs + i + 0.5) + 1;
            Deyr[i] = delays(tpie, tpjs + i + 0.5) + 1;
        }
        //H incident is in scattered region
        for (i = 0; i < Max_Index_Dhz_x; i++) {
            Dhzb[i] = delays(tpis + i + 0.5, tpjs - 0.5);
            Dhzt[i] = delays(tpis + i + 0.5, tpje + 0.5);
        }
        for (i = 0; i < Max_Index_Dhz_y; i++) {
            Dhzl[i] = delays(tpis - 0.5, tpjs + i + 0.5);
            Dhzr[i] = delays(tpie + 0.5, tpjs + i + 0.5);
        }
    }

}

void InitConnectingInterface(MyDataF IncAngle) {
    //define position of the connecting interface
    tpis = pis + scatwidth;
    tpjs = pjs + scatwidth;

    tpie = pie - scatwidth;
    tpje = pje - scatwidth;

    CMur = (c * dt - dx) / (c * dt + dx);

    cos_phi = cos(IncAngle);
    sin_phi = sin(IncAngle);

    MyDataF vl = PhaseVelRatio(0) / PhaseVelRatio(IncAngle);
    //updating coefficients
    if (IsTMz) {
        Chzey = -dt / mu_0 / dx*cos_phi;
        Chzex = -dt / mu_0 / dy*sin_phi;

        Ceyhz = -dt / eps_0 / dx;
        Cexhz = dt / eps_0 / dy;

        ceihi = -dt / eps_0 / dx / vl;
        chiei = -dt / mu_0 / dx / vl;
    }

    if (IsTEz) {
        ceihi = dt / eps_0 / dx / vl;
        chiei = dt / mu_0 / dx / vl;
    }

    //Create space for Line EM fields
    eilen = (int) ceil(sqrt((double) ((tpie - tpis)*(tpie - tpis)+(tpje - tpjs)*(tpje - tpjs)))) + 2;
    hilen = eilen - 1;

    Ei = new MyDataF[eilen]();
    Hi = new MyDataF[hilen]();
    Ei[0] = E0 * Source(0);
    //Calculate Delay at connecting interface
    CalculateDelay(IncAngle);

    //#ifdef  MATLAB_SIMULATION
    //	mxArray *MyArray;
    //	MyArray=mxCreateDoubleMatrix(1,eilen,mxREAL);
    //	memcpy(mxGetPr(MyArray), Ei,eilen*sizeof(MyDataF));
    //	engPutVariable(ep,"Ei",MyArray);
    //	engEvalString(ep,"figure,eih=plot(Ei);grid on;");//%surf(data1);");//
    //	mxDestroyArray(MyArray);
    //#endif
}

void UpdateConnectingE(const MyDataF t) {
    int i;
    //int di;
    MyDataF ei_last, ei_last2;

    ei_last = Ei[eilen - 1];
    ei_last2 = Ei[eilen - 2];
    Ei[0] = E0 * Source(t);
    for (i = 1; i <= eilen - 2; i++)
        Ei[i] = Ei[i] + ceihi * (Hi[i] - Hi[i - 1]);

    //mur boundary 
    Ei[eilen - 1] = ei_last2 + CMur * (Ei[eilen - 2] - ei_last);
}

void UpdateConnectingM(const MyDataF t) {
    int i;
    //    MyDataF hi_last, hi_last2;

    //    hi_last = Hi[hilen - 1];
    //    hi_last2 = Hi[hilen - 2];
    for (i = 0; i <= hilen - 1; i++)
        Hi[i] = Hi[i] + chiei * (Ei[i + 1] - Ei[i]);
    //Mur boundary
    //Hi[hilen-1]=hi_last2+CMur*(Hi[hilen-2]-hi_last);
    //#ifdef  MATLAB_SIMULATION
    //	mxArray *MyArray;
    //	MyArray=mxCreateDoubleMatrix(1,eilen,mxREAL);
    //	memcpy(mxGetPr(MyArray), Ei,eilen*sizeof(MyDataF));
    //	engPutVariable(ep,"Ei",MyArray);
    //	engEvalString(ep,"set(eih,'ydata',Ei);");//surf(data1);");//
    //	//engEvalString(ep,"title(ind);ind=ind+1;");
    //	mxDestroyArray(MyArray);
    //#endif
}

void ApplyConnectingE(const MyDataF t) {
    MyDataF df, d;
    int di;
    int ind;
    UpdateConnectingE(t);
    if (IsTEz) {
        /*
              unsigned int index,ind;

              for(ind = tpis;ind<=tpie;ind++){
                      index = ind*Ez.ny;
                      //bottom 
                      d		=	Dhxb[ind+start_index_x]+0.5;
                      di		=	(int)floor(d);
                      df		=	d - di;
                      Ez.data[ind][tpjs]  +=	Cehz.data[ind][tpjs]*Ratio_x*(Hi[di]+df*(Hi[di+1]-Hi[di]))/dy;

                      //top
                      d		=	Dhxt[ind+start_index_x]+0.5;
                      di		=	(int)floor(d);
                      df		=	d - di;
                      Ez.data[ind][tpje] -=	Cehz.data[ind][tpje]*Ratio_x*(Hi[di]+df*(Hi[di+1]-Hi[di]))/dy;
              }
              //bound xn,xp
              for(ind = tpjs;ind<=tpje;ind++){
                      //left side
                      index	=	tpis*Ez.ny+ind;
                      d		=	Dhyl[ind+start_index_y]+0.5;
                      di		=	(int)floor(d);
                      df		=	d - di;
                      Ez.data[tpis][ind] 	-=	 Cehz.data[tpis][ind]*(-Ratio_y)*(Hi[di]+df*(Hi[di+1]-Hi[di]))/dx;

                      //right side
                      index	=	tpie*Ez.ny+ind;
                      d		=	Dhyr[ind+start_index_y]+0.5;
                      di		=	(int)floor(d);
                      df		=	d - di;
                      Ez.data[tpie][ind]	+=	Cehz.data[tpie][ind]*(-Ratio_y)*(Hi[di]+df*(Hi[di+1]-Hi[di]))/dx;
              }
         */
    }
    if (IsTMz) {
        unsigned int i, j;
        //Ceyhz = dt/eps_0/dx;
        for (ind = 0, j = tpjs; j < tpje; j++, ind++) {
            d = Dhzl[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ey.data[tpis][j] -= Ceyhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //left side

            d = Dhzr[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ey.data[tpie][j] += Ceyhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //right side

        }
        //Cexhz = dt/eps_0/dy;
        for (ind = 0, i = tpis; i < tpie; i++, ind++) {
            d = Dhzb[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ex.data[i][tpjs] -= Cexhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //bottom side

            d = Dhzt[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ex.data[i][tpje] += Cexhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //top side
        }
    }
}

void ApplyConnectingM(const MyDataF t) {

    unsigned int ind, ind1;
    int di;
    MyDataF df;
    UpdateConnectingM(t);
    if (IsTEz) {
        /*
                  for(ind=tpjs;ind<=tpje;ind++){
                          //left side
                          di = (int)floor(Dezl[ind+start_index_y]);
                          df = Dezl[ind+start_index_y] - di;
                          di = di+1;
                          Hy.data[(tpis-1)][ind]	-=	Chyez*(Ei[di]+df*(Ei[di+1]-Ei[di]));//0;//
                          //right side
                          di = (int)floor(Dezr[ind+start_index_y]);
                          df = Dezr[ind+start_index_y] - di;
                          di = di+1;
                          Hy.data[tpie][ind]		+=	Chyez*(Ei[di]+df*(Ei[di+1]-Ei[di]));
                  }
                  //Adjust Hy
                  for(ind=tpis;ind<=tpie;ind++){
                          // bottom
                          di = (int)floor(Dezb[ind+start_index_x]);
                          df = Dezb[ind+start_index_x] - di;
                          di = di+1;

              Hx.data[ind][tpjs-1]		-=		Chxez*(Ei[di]+df*(Ei[di+1]-Ei[di]));
                          //yp
                          di = (int)floor(Dezt[ind+start_index_x]);
                          df = Dezt[ind+start_index_x] - di;
                          di = di+1;
              Hx.data[ind][tpje]		+=	Chxez*(Ei[di]+df*(Ei[di+1]-Ei[di]));
                  }
         */
    }
    if (IsTMz) {
        unsigned int mtpis = tpis - 1;
        unsigned int mtpjs = tpjs - 1;
        //Chzex = dt/mu_0/dy;
        for (ind = tpis, ind1 = 0; ind < tpie; ind++, ind1++) {
            di = (int) floor(Dexb[ind1]);
            df = Dexb[ind1] - di;
            //under (i+1/2,tpjs):bottom,j0
            Hz.data[ind][mtpjs] -= Chzex * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //bottom side
            di = (int) floor(Dext[ind1]);
            df = Dext[ind1] - di;
            //top,j1
            Hz.data[ind][tpje] += Chzex * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //top side
        }
        //Chzey = dt/mu_0/dx;
        for (ind = tpjs, ind1 = 0; ind < tpje; ind++, ind1++) {
            di = (int) floor(Deyl[ind1]);
            df = Deyl[ind1] - di;
            //left,i0
            Hz.data[mtpis][ind] -= Chzey * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //left side
            di = (int) floor(Deyr[ind1]);
            df = Deyr[ind1] - di;
            //right,i1
            Hz.data[tpie][ind] += Chzey * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //right side
        }
    }
}

void FreeDelayArrays() {
    if (IsTEz) {
        delete[] Dhyl;
        delete[] Dhxb;
        delete[] Dhxt;
        delete[] Dhyr;
        delete[] Dezl;
        delete[] Dezr;
        delete[] Dezt;
        delete[] Dezb;
    }
    if (IsTMz) {
        delete[] Deyl;
        delete[] Dexb;
        delete[] Dext;
        delete[] Deyr;
        delete[] Dhzl;
        delete[] Dhzr;
        delete[] Dhzt;
        delete[] Dhzb;
    }
    delete[] Ei;
    delete[] Hi;
}

MyDataF Source(MyDataF t) {
    //MyDataF t_0;
    //MyDataF tau = 3.1E-9;
    //t_0 = 0.8 * tau; //dt*ttstep/10;//||t>1000*dt

    if (t < 0)return 0;
    return sin(omega * t); //-cos(2 * M_PI * f * t) * exp(-4 * M_PI * pow((t - t_0) / tau, 2)); //1;//

    //if (fabs(t/dt-10)<1e-3)
    //	return 1.0;
    //else return 0;
}
