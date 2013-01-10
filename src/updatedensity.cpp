/*
 *		This file defines functions that update the density.
 *
 *
 *
 *
 */
#include <iostream>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <iomanip>


#include "datastruct.h"
#include "commondata.h"
#include "updatefields.h"

using namespace std;

extern unsigned tpis, tpjs, tpie, tpje;
const MyDataF v2Div3 = 2.00000000000/3.000000000;
void updateDeff() {
    unsigned int i, j, mt = m2*tpis;
    MyDataF eps = 0;
    MyDataF neij, kasi;

    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            neij = Ne.data[i][j];
            eps = ElectricEnergy(i, j);
            Niu_i.data[i][j] = GasDen * vi.Interp(eps);
            Niu_a.data[i][j] = GasDen * va.Interp(eps);

            if (neij < 1e-5) {
                Deff.data[i][j] = De;
            } else {
                kasi = eps_m_e_miu * Niu_i.data[i][j] / neij;
                Deff.data[i][j] = (kasi * De + Da) / (kasi + 1);
            }
        }
    }
}

void StoreOpt4(unsigned i, unsigned j, MyDataF opt2, MyDataF opt4) {
    //cout<<opt4<<endl;
    if (j == midj && i == pci)
        cout << opt4 << '\t';
    for (unsigned c = 0; c < 10; c++) {
        if (i == Deff_Store_Index_x[c] && j == Deff_Store_Index_y[c]) {
            deff_file << opt2 << '\t' << opt4 << '\t';
            denfile << Ne.data[i][j] << '\t' ;
#ifdef _DEBUG
            cout << "store field ("<<c<<"):" << i << "," << j << endl;
#endif
            break;
        }
    }
}

/**
 * deff
 *
 *
 */
void UpdateDenDeff() {
    unsigned int i, j, mt = m2*tpis;
    MyDataF opt1 = 0, opt2 = 0, opt3 = 0, opt4 = 0;
    unsigned mi = 0, mj = 0;
    MyDataF maxne = 0;
    Pne = Ne;
    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {

            opt1 = 1 + dt * Niu_i.data[i][j];
            opt2 = dt_ds2_2 * Deff.data[i][j] * (Pne.data[i - 1][j] + Pne.data[i + 1][j] +
                                                 Pne.data[i][j + 1] + Pne.data[i][j - 1] - 4 * Pne.data[i][j]
                                                );
            opt3 = 1 + dt * (Niu_a.data[i][j] + rei * Ne.data[i][j]);
            opt4 = ((Deff.data[i + 1][j] - Deff.data[i][j])*(Ne.data[i + 1][j] - Ne.data[i][j])+
                    (Deff.data[i][j + 1] - Deff.data[i][j])*(Ne.data[i][j + 1] - Ne.data[i][j])) * dt_ds2_2;
            Ne.data[i][j] = (Ne.data[i][j] * opt1 + opt2 + opt4) / opt3;
            //if (j == midj || i == pci)
            //   StoreOpt4(i, j, opt2, opt4);
            if(i>=minSI&&j>=minSJ&&i<=maxSI&&j<=maxSJ)
                StoreOpt4(i, j, opt2, opt4);
            if (maxne < Ne.data[i][j]) {
                mi = i;
                mj = j;
                maxne = Ne.data[i][j];
            }
        }
    }
    DensityBound(Ne, m*tpis, 0);
    cout << maxne << '\t' << mi << '\t' << mj << '\t';
    cout << Ne.data[Deff_Store_Index_x[5]][Deff_Store_Index_y[5]] << '\t' << Ne.data[midi][pci] << '\t';
    deff_file << endl;
    denfile << endl;
#ifdef _DEBUG
    cout << "===================== end deff file ================" << endl;
#endif
}
// update density with formular from Zhao at 2012-10
// with Boltzmann equation
// 2012-10-23
// by Baofeng Shi
// email: skiloop@126.com
void UpdateDensity1210() {
    unsigned int i, j, mt = m2*tpis;
    MyDataF opt1 = 0, opt2 = 0, opt3 = 0 ;//, opt4 = 0;
    unsigned mi = 0, mj = 0;
    MyDataF maxne = 0;
    MyDataF EmDivN;
    MyDataF niu_i,niu_a,energy,Te;
    MyDataF deff,kasi;
#ifdef DEBUG
    unsigned midx=Ne.nx/2;
    unsigned midy=Ne.ny/2;
#endif
    Pne = Ne;
    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            // get EmDivN at (i,j)
            EmDivN=getEmax(i,j);
            // compute niu_i
            niu_i = ViDivN.Interp(EmDivN)*N_air;
            // compute niu_a
            niu_a = VcDivN.Interp(EmDivN)*N_air;
            // compute energy
            energy = EnergyDivN.Interp(EmDivN);

            // calculate Tempreture
            Te = v2Div3*energy;

            // calculate effective diffusion coefficient deff
	    if (Ne.data[i][j]<=0.0){
		    deff = miu2DivE*Te;//mu_e*2*Te/e;//8.73e-2; //
	    }else{
		    kasi = eps_m_e_miu * niu_i /Ne.data[i][j];
		    De = miu2DivE*Te;//mu_e*2*Te/e;//8.73e-2; //
		    //De = mu_e*2*Te/e;//8.73e-2; //
		    Da = De / MueDivMui;
		    deff = (kasi*De + Da)/(kasi + 1);
	    }


            // calculate Ne at (i,j)
            opt1 = 1+dt_F*niu_i;
            opt2 = deff*dt_F*(Pne.data[i][j+1]+Pne.data[i-1][j]+Pne.data[i+1][j]+Pne.data[i][j-1]-4*Ne.data[i][j])/ds_F/ds_F;
            opt3 = 1+dt_F*(niu_a+rei*Ne.data[i][j]);
            Ne.data[i][j] = (Ne.data[i][j]*opt1+opt2)/opt3;
#ifdef DEBUG
	    //if (i==midx&&j==midy)
	    //    cout << "Debug:" << getEmax(i,j)<<'\t'<< EmDivN <<'\t'<<niu_i<<'\t'<<niu_a<<'\t'<<energy <<endl;
	    //if (i==Deff_Store_Index_x[5]&&j==Deff_Store_Index_y[5])
	    //   cout << Ne.data[i][j] << endl;
#endif

            if (maxne < Ne.data[i][j]) {
                mi = i;
                mj = j;
                maxne = Ne.data[i][j];
            }
        }
    }
    DensityBound(Ne, m*tpis, 0);
    //cout << maxne << '\t' << mi << '\t' << mj << '\t';
    //cout << Ne.data[Deff_Store_Index_x[5]][Deff_Store_Index_y[5]] << '\t' << Ne.data[midi][pci] << '\t';
    //cout << "Nemiddle " << Ne.data[Ne.nx/2][Ne.ny/2] << '\t';
    deff_file << endl;
    denfile << endl;
#ifdef _DEBUG
    cout << "===================== end deff file ================" << endl;
#endif
}


void UpdateDensityOther() {
    unsigned int i, j, mt = m2*tpis;
    MyDataF niu_a, niu_i, eps, deff;
    MyDataF maxne = 0;
    MyDataF minne = 0;
    unsigned ci = 0, cj = 0;
    MyDataF maxvi = 0, va_maxvi = 0, ne_maxvi = 0;
    MyDataF neij, kasi, gamma1 = 0, gamma2 = 0, down = 1;
    //MyDataF mvi = 0, mva = 0, mga1 = 0, mga2 = 0;
    //MyDataF mue, mneij;
    MyDataF meps = 0;
    MyDataF opt1 = 0, opt2 = 0, opt3 = 0;

    if (denFormula == 2)
        ppne = Pne;

    Pne = Ne;
    for (i = mt; i < Ne.nx - mt; i++) {
        for (j = mt; j < Ne.ny - mt; j++) {
            neij = Ne.data[i][j];
            eps = ElectricEnergy(i, j);
            niu_i = GasDen * vi.Interp(eps);
            niu_a = GasDen * va.Interp(eps);
            //niu_i = neij * vi.Interp(eps);
            //niu_a = neij * va.Interp(eps);
            if (neij < 1e-5) {
                deff = De;
            } else {
                kasi = eps_m_e_miu * niu_i / neij;
                deff = (kasi * De + Da) / (kasi + 1);
            }
            switch (denFormula) {
            case 1:
                down = (2 - dt * (niu_i - niu_a));
                gamma1 = (2 + dt * (niu_i - niu_a)) / down;
                gamma2 = (dt2 * deff) / down;
                Ne.data[i][j] = Ne.data[i][j] * gamma1 +
                                gamma2 * (Pne.data[i - 1][j] + Pne.data[i + 1][j] +
                                          Pne.data[i][j + 1] + Pne.data[i][j - 1] - 4 * Pne.data[i][j]
                                         );
                 break;
                 case 2:
                 Ne.data[i][j] = Ne.data[i][j] * dt2 * (niu_i - niu_a) - ppne.data[i][j] +
                 dt_ds2_2 * deff * (Pne.data[i - 1][j] + Pne.data[i + 1][j] +
                 Pne.data[i][j + 1] + Pne.data[i][j - 1] - 4 * Pne.data[i][j]
                 );
                 break;
                 case 3:
                 break;
                 case 4:
                 break;
                 default:
                 //Bhaskar Chaudhury's formula
                 opt1 = 1 + dt_F*niu_i;
                 opt2 = DtfDivDsfs * deff * (Pne.data[i - 1][j] + Pne.data[i + 1][j] +
                                             Pne.data[i][j + 1] + Pne.data[i][j - 1] - 4 * Pne.data[i][j]
                                            );
                 opt3 = 1 + dt_F * (niu_a + rei * neij);
                 Ne.data[i][j] = (neij * opt1 + opt2) / opt3;
                 break;
            }
             if (Ne.data[i][j] < 0)
                 Ne.data[i][j] = 0;
            if (Ne.data[i][j] > maxne) {
                maxne = Ne.data[i][j];
            }
            if (eps > meps) {
                //mue = Ue.data[i][j];
                //mneij = neij;
                //mvi = niu_i;
                //mva = niu_a;
                //mga1 = gamma1;
                //mga2 = gamma2;
                meps = eps;
            }
            if (Ne.data[i][j] < minne)
                minne = Ne.data[i][j];
            if (maxvi < niu_i) {
                ne_maxvi = Ne.data[i][j];
                ci = i;
                cj = j;
                va_maxvi = niu_a;
                maxvi = niu_i;
            }
        }
    }

    DensityBound(Ne, m*tpis, 0);
    denfile << setiosflags(ios_base::scientific)<<maxne << '\t' << minne << '\t' << maxvi << '\t' << ne_maxvi << '\t' << va_maxvi << '\t' << ci << '\t' << cj << endl;
    cout << setiosflags(ios_base::scientific)<<maxne << '\t' << minne << '\t' << maxvi << '\t' << ne_maxvi << '\t' << va_maxvi << '\t' << ci << '\t' << cj << '\t' << meps << '\t';
    cout << "Nemiddle" << Ne.data[Ne.nx/2][Ne.ny/2] << '\t';
    /*
        cout << endl;
        cout << "meps\t" << meps << endl;
        cout << "mneij\t" << mneij << endl;
        cout << "mga2\t" << mga2 << endl;
        cout << "mga1\t" << mga1 << endl;
        cout << "mvi\t" << mvi << endl;
        cout << "mva\t" << mva << endl;
        cout << "mue\t" << mue << endl;
     */
    //PrintData(Ne);
    //if(max_dtmui>1.0)
    //cout << max_dtmui << "\t" << endl;
    //system("pause");

}

void UpdateDensity() {
    switch (denFormula) {
    case 3:
        updateDeff();
        UpdateDenDeff();
        break;
    case 4:
        //InterpEmax();
        UpdateDensity1210();
        break;
    default:
        UpdateDensityOther();
    }
}

void InterpUe(unsigned InterStep) {

    unsigned i, j;
    unsigned im, jm;
    unsigned jsw, jew;
    unsigned ilast = Ue.nx - InterStep;
    unsigned jlast = Ue.ny - InterStep;
    unsigned InterStep2 = InterStep*InterStep;
    if (InterStep == 0) {
        cerr << "Invalid interpolate step,it must be larger than 0" << endl;
        exit(-1);
    }
    unsigned ime, jme; //end of i and j
    for (im = InterStep; im < ilast; im = ime) {
        ime = im + InterStep;
        for (jm = InterStep; jm < jlast; jm = jme) {
            jme = jm + InterStep;
            for (i = im + 1; i < ime; i++) {
                unsigned iew = i - im; //weight of point ie
                unsigned isw = ime - i; //weight of point is
                for (j = jm + 1; j < jme; j++) {
                    jew = j - jm; //weight of point je
                    jsw = jme - j; //weight of point js
                    //interpolate Ue which not on coarse lines
                    Ue.data[i][j] = (isw * jsw * Ue.data[im][jm] + iew * jsw * Ue.data[ime][jm] +
                                     isw * jew * Ue.data[im][jme] + iew * jew * Ue.data[ime][jme]) / InterStep2;

                }
                //interpolate Ue at lines jm and jme
                Ue.data[i][jm] = (isw * Ue.data[im][jm] + iew * Ue.data[ime][jm]) / InterStep;
                //Ue.data[i][jme] = (isw * Ue.data[im][jme] + iew * Ue.data[ime][jme]) / InterStep;
            }
            //interpolate Ue at lines im and ime
            for (j = jm + 1; j < jme; j++) {
                jew = j - jm;
                jsw = jme - j;
                Ue.data[im][j] = (jsw * Ue.data[im][jm] + jew * Ue.data[im][jme]) / InterStep;
                //Ue.data[ime][j] = (jsw * Ue.data[ime][jm] + jew * Ue.data[ime][jme]) / InterStep;
            }
        }
    }
    //cout<<"%4.4e\t%4.4e\t%4.4e\t%4.4e\n",Ue.data[cnepx*Ue.ny+cnepy],Ue.data[cnepx*Ue.ny+cnepy+m],
    //	Ue.data[cnepx*Ue.ny+cnepy+m*Ue.ny],Ue.data[cnepx*Ue.ny+cnepy+m*Ue.ny+m]<<endl;
}

void UpdateUeField() {
    unsigned i, j, im, jm;
    unsigned ilast = Ue.nx - m;
    unsigned jlast = Ue.ny - m;
    MyDataF dtQe_8 = dt / 8.0;
    MyDataF Q = 0;
    MyDataF eps;
    //MyDataF maxue = 0;
    // <editor-fold defaultstate="collapsed" desc="Update at coarse point (i+1/2,j+1/2)">
    //    //Ue at (i+1/2,j+1/2)
    //    for (i = 0, im = m2; im < ilast; i++, im += m) {
    //        for (j = 0, jm = m2; jm < jlast; j++, jm += m) {
    //
    //            eps = ElectricEnergy(im, jm);
    //            Q = GasDen * lossE.Interp(eps);
    //            //Q = Ne.data[im][jm] * lossE.Interp(eps);
    //            if (jm <= 325 && jm >= 300)Ue.data[im][jm] = 100;
    //            else
    //                Ue.data[im][jm] += dtQe_8 * (
    //                    (Ex.data[i][j] + Ex.data[i][j + 1] + Pex.data[i][j] + Pex.data[i][j + 1])*
    //                    (Ux.data[i][j + 1] + Ux.data[i][j])
    //                    +(Ey.data[i][j] + Ey.data[i + 1][j] + Pey.data[i][j] + Pey.data[i + 1][j])*
    //                    (Uy.data[i][j] + Uy.data[i + 1][j]))
    //                - half_dt * Q * (Ne.data[im][jm] + Pne.data[im][jm]);
    //        }
    //    }// </editor-fold>

    // <editor-fold defaultstate="collapsed" desc="Update Ue at coarse point (i,j)">
    //Ue at (i,j)
    for (i = 1, im = m; im < ilast; i++, im += m) {
        for (j = 1, jm = m; jm < jlast; j++, jm += m) {

            eps = ElectricEnergy(im, jm);
            //Q = Ne.data[im][jm] * lossE.Interp(eps);
            Q = GasDen * lossE.Interp(eps);
            //            if ((im <= 325 && im >= 300) || (jm <= 325 && jm >= 200))
            //                Ue.data[im][jm] = jm;
            //            else
            Ue.data[im][jm] += dtQe_8 * (
                                   (Ex.data[i][j] + Ex.data[i - 1][j] + Pex.data[i][j] + Pex.data[i - 1][j])*
                                   (Ux.data[i - 1][j] + Ux.data[i][j])
                                   +(Ey.data[i][j - 1] + Ey.data[i][j] + Pey.data[i][j] + Pey.data[i][j - 1])*
                                   (Uy.data[i][j] + Uy.data[i][j - 1]))
                               - half_dt * Q * (Ne.data[im][jm] + Pne.data[im][jm]);
            //if (Ue.data[im][jm] > maxue)
            //maxue = Ue.data[im][jm];
        }
    }// </editor-fold>
    //cout << "Max Ue:" << maxue << endl;
    //interpolate Ue
    InterpUe(m);
}

void DensityBound(MyStruct stru, int bndwidth, const int swidth) {
    /*
    //Checked at start3.25 11:14

    //////////////////////////////////////
    //			//		//
    //	2		//		//
    //			//		//
    //////////////////////              //
    //		//	//		//
    //		//	//		//
    //		//	//	3	//
    //	1	//	//		//
    //		//	//		//
    //		//	//		//
    //		//////////////////////////
    //		//			//
    //		//		4	//
    //		//			//
    //////////////////////////////////////

     */

    int i, j; //, ind;
    int ii, i1, i2;
    double tmp;
    ii = swidth + bndwidth;
    i1 = stru.nx - swidth - bndwidth;
    i2 = stru.nx - swidth;
    //Section 1
    for (j = swidth; j < i1; j++) {
        //ind =  ii * stru.ny + j;
        tmp = 0; //2*stru.data[ind+stru.ny]-stru.data[ind+2*stru.ny];
        for (i = swidth; i <= ii; i++)
            stru.data[i][j] = tmp;
    }
    //Section 2
    for (i = 0; i < i1; i++) {
        //ind =  i * stru.ny + i1;
        tmp = 0; //2*stru.data[ind-1]-stru.data[ind-2];
        for (j = i2 - 1; j >= i1; j--)
            stru.data[i][j] = tmp;
    }
    //Section 3
    for (j = ii; j < i2; j++) {
        //ind =  i1 * stru.ny + j;
        tmp = 0; //2*stru.data[ind-stru.ny]-stru.data[ind-2*stru.ny];
        for (i = i2 - 1; i >= i1; i--)
            stru.data[i][j] = tmp;
    }

    //Section 4
    for (i = i1; i < i2; i++) {
        //ind =  ii + i * stru.ny;
        tmp = 0; //2*stru.data[ind+1]-stru.data[ind+2];
        for (j = swidth; j <= ii; j++)
            stru.data[i][j] = tmp;
    }
}

