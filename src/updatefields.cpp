/*
 *		This file defines functions that update scatter fields.
 *
 *
 *
 *
 */

#include <iostream>
#include <math.h>
#include "commondata.h"
#include "updatefields.h"

extern unsigned pis, pie, pjs, pje;

void InterpEmax() {
    unsigned i,j;
    unsigned im,jm;
    for(i=m2; i<Em.nx-m2; i++) {
        for(j=m2; j<Em.ny-m2; j++) {
            if((i+m2)%m!=0 && (j+m2)%m!=0) {
                im = ((i-m2)/m)*m+m2;
                jm = ((j-m2)/m)*m+m2;
                Em.data[i][j] = ((im+m-i)*(jm-m-j)*Em.data[im][jm]
                                 +(im+m-i)*(j-jm)*Em.data[im][jm+m]+(i-im)*(jm+m-j)*Em.data[im+m][jm]+(im+m-i)*(jm+m-j)*Em.data[im+m][jm+m])/(m*m);
            }
        }
    }

}
void CalEmax() {
    //unsigned xi,xj,yi,yj;
    unsigned i,j;
    unsigned im,jm;
    MyDataF eabs;
#ifdef DEBUG
    MyDataF maxEm=0;
    unsigned mi = 0;
    unsigned mj = 0;
#endif

    if(IsTMz) {
        MyDataF emx,emy;
        for(im=m2,i=0; im<Em.nx-m2; i++,im+=m) {
            for(jm=m2,j=0; jm<Em.ny-m2; j++,jm+=m) {
                emx = (Ex.data[i][j]+Ex.data[i][j+1])/2;
                emy = (Ey.data[i][j]+Ey.data[i+1][j])/2;
                eabs=sqrt(emx*emx+emy*emy);
                if(eabs>Em.data[im][jm]){
			Em.data[im][jm] = eabs;
#ifdef DEBUG
			if (eabs>maxEm){
				maxEm=eabs;
				mi = im;
				mj = jm;
			}
#endif
		}
            }
        }
    }
    if(IsTEz) {
        for(im=0,i=0; im<Em.nx; i++,im+=m) {
            for(jm=0,j=0; jm<Em.ny; j++,jm+=m) {
                eabs = fabs(Hz.data[i][j]);
                if(eabs>Em.data[im][jm]){
                    Em.data[im][jm] = eabs;
#ifdef DEBUG
		    if (eabs>maxEm){
			    maxEm=eabs;
			    mi = im;
			    mj = jm;
		    }
#endif
		}
            }
        }
    }
#ifdef DEBUG
    cout << "maxEm: " << maxEm << '\t' << mi << '\t' << mj << endl;
#endif
}

MyDataF getEmax(int i, int j) {
    return Em.data[i][j];
}


void CapFields(unsigned int step) {
    //std::cout << "==========Capture fields the " << step << " times ============" << std::endl;
    //Ex.CaptData(step);
    //Ey.CaptData(step);
    Hz.CaptData(step);
    if(ifWithDensity) {
        //Ue.CaptData(step);
        Ne.CaptData(step);
    }
}

void UpdateEField() {

    unsigned int i, j;

    if (IsTMz) {
        Pex = Ex;
        Pey = Ey;
        for (i = 0; i < Ex.nx; i++) {
            for (j = pjs + 1; j < pje; j++) {
                Ex.data[i][j] += cexhz * (Hz.data[i][j] - Hz.data[i][j - 1])
                                 + cexux * Ux.data[i][j];
            }
        }
        for (i = 1 + pis; i < pie; i++) {
            for (j = 0; j < Ey.ny - 1; j++) {
                Ey.data[i][j] += ceyhz * (Hz.data[i][j] - Hz.data[i - 1][j])
                                 + ceyuy * Uy.data[i][j];
            }
        }
    }

}

void UpdateMField() {

    unsigned int i, j;

    for (i = pis; i < pie; i++) {
        for (j = pjs; j < pje; j++) {
            Hz.data[i][j] +=
                chzex * (Ex.data[i][j + 1] - Ex.data[i][j]) +
                chzey * (Ey.data[i + 1][j] - Ey.data[i][j]);
        }
    }
}

void UpdateUField() {
    unsigned i, j, im, jm;
    MyDataF tvc, eps;
    MyDataF alpha, beta;

    for (i = 0, im = m2; i < Ux.nx; im += m, i++) {
        for (j = 0, jm = 0; j < Ux.ny; jm += m, j++) {
            eps = ElectricEnergy(im, jm);
            tvc = GasDen * vc.Interp(eps) * dt;
            //tvc = Ne.data[im][jm] * vc.Interp(eps);
            alpha = (2.0 - tvc) / (2.0 + tvc);
            beta = dt_me_e_2 / (2.0 + tvc);
            Ux.data[i][j] = alpha * Ux.data[i][j] + beta * Ne.data[im][jm] * Ex.data[i][j];
        }
    }
    for (i = 0, im = 0; i < Uy.nx; im += m, i++) {
        for (j = 0, jm = m2; j < Uy.ny; j++, jm += m) {
            eps = ElectricEnergy(im, jm);
            tvc = GasDen * vc.Interp(eps) * dt;
            //tvc = Ne.data[im][jm] * vc.Interp(eps);
            alpha = (2.0 - tvc) / (2.0 + tvc);
            beta = dt_me_e_2 / (2.0 + tvc);
            Uy.data[i][j] = alpha * Uy.data[i][j] + beta * Ne.data[im][jm] * Ey.data[i][j];
        }
    }

}

void SumErms() {
    /*
    unsigned i, j;
    MyDataF tEx, tEy;

    if (IsTMz) {
        for (i = 1; i < Erms.nx - 1; i++) {
            for (j = 1; j < Erms.ny - 1; j++) {
                tEx = 0.5 * (Ex.data[i - 1][j] + Ex.data[i][j]);
                tEy = 0.5 * (Ey.data[i][j] + Ey.data[i][j - 1]);

                Erms.data[i][j] += tEx * tEx + tEy*tEy;
            }
        }
    }
    if (IsTEz) {
        for (i = 1; i < Erms.nx - 1; i++) {
            for (j = 1; j < Erms.ny - 1; j++) {
                Erms.data[i][j] += Ez.data[i][j] * Ez.data[i][j];
            }
        }
    }
     */
}

void UpdateErms() {
    /*
    unsigned i, j;
    static MyDataF mdt = dt / dt_M;
    for (i = 0; i < Erms.nx; i++) {
        for (j = 0; j < Erms.ny; j++) {
            Erms.data[i][j] = sqrt(mdt * Erms.data[i][j]);
        }
    }
     */
}

void UpdateUCoeffiecients() {
    /*unsigned i, j;
    MyDataF a;
    const MyDataF tQ = e / me;
    for (i = 0; i < Cua.nx; i++) {
        for (j = 0; j < Cua.ny; j++) {
            a = Niu_c.data[i][j] * dt / 2;
            Cua.data[i][j] = (1 - a) / (1 + a);
            Cub.data[i][j] = tQ * dt / (1 + a);
        }
    }
     */
}


