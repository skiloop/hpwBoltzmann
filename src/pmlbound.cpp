#include <iostream>
#include <cstdlib>
#include "pmlbound.h"
#include "commondata.h"
#include "bounddata.h"
using namespace std;

/**
 * @brief define parameters of the pml bound
 * @param ncpml number of cells in pml boundaries
 */
void boundPMLdef(unsigned ncpml) {
    pis = ncpml;
    pjs = ncpml;
    pie = nx - ncpml;
    pje = ny - ncpml;

    IsPMLxn = 1;
    IsPMLxp = 1;
    IsPMLyn = 1;
    IsPMLyp = 1;

    PMLCellxn = ncpml;
    PMLCellxp = ncpml;
    PMLCellyn = ncpml;
    PMLCellyp = ncpml;

    IsAnySidePML = 1;
}

/**
 * Create PML data
 * @param ncx number of cells in x dimension
 * @param ncy number of cells in y dimension
 * @brief create space and initialize the data of PML arrays
 * @return tree if success,or false if invalid inputs and not enough space
 */
bool CreatePMLArrays(unsigned int ncx, unsigned int ncy) {
    if (IsAnySidePML) {
        if (IsPMLxn) {
            if (IsTEz) {
                Ezx_xn.CreateStruct(PMLCellxn, nym1);
                Ezy_xn.CreateStruct(PMLCellxn, nym1 - PMLCellyn - PMLCellyp);
            }
            if (IsTMz) {
                Hzx_xn.CreateStruct(PMLCellxn, ny);
                Hzy_xn.CreateStruct(PMLCellxn, ny - PMLCellyn - PMLCellyp);
            }
        }
        if (IsPMLxp) {
            if (IsTEz) {
                Ezy_xp.CreateStruct(PMLCellxp, nym1 - PMLCellyn - PMLCellyp);
                Ezx_xp.CreateStruct(PMLCellxp, nym1);
            }
            if (IsTMz) {
                Hzx_xp.CreateStruct(PMLCellxp, ny);
                Hzy_xp.CreateStruct(PMLCellxp, ny - PMLCellyp - PMLCellyn);
            }
        }
        if (IsPMLyn) {
            if (IsTEz) {
                Ezy_yn.CreateStruct(nxm1, PMLCellyn);
                Ezx_yn.CreateStruct(nxm1 - PMLCellxn - PMLCellxp, PMLCellyn);

            }
            if (IsTMz) {
                Hzx_yn.CreateStruct(nx - PMLCellxn - PMLCellxp, PMLCellyn);
                Hzy_yn.CreateStruct(nx, PMLCellyn);
            }
        }
        if (IsPMLyp) {
            if (IsTEz) {
                Ezy_yp.CreateStruct(nxm1, PMLCellyp);
                Ezx_yp.CreateStruct(nxm1 - PMLCellxn - PMLCellxp, PMLCellyp);
            }
            if (IsTMz) {
                Hzx_yp.CreateStruct(nx - PMLCellxp - PMLCellxn, PMLCellyp);
                Hzy_yp.CreateStruct(nx, PMLCellyp);
            }
        }
    }
    return true;
}

/**
 * Create PML Coefficients data
 * @param ncx number of cells in x dimension
 * @param ncy number of cells in y dimension
 * @brief create space and initialize the data of PML arrays
 * @return tree if success,or false if invalid inputs and not enough space
 */
bool CreatePMLCoeffients(unsigned int ncx, unsigned int ncy) {
    if (IsAnySidePML) {
        if (IsPMLxn) {
            if (IsTEz) {
                Chyh_xn.CreateStruct(PMLCellxn, nym1);
                Chyez_xn.CreateStruct(PMLCellxn, nym1);
                Cezxe_xn.CreateStruct(Ezx_xn);
                Cezxhy_xn.CreateStruct(Ezx_xn);
                Cezye_xn.CreateStruct(Ezy_xn);
                Cezyhx_xn.CreateStruct(Ezy_xn);
            }
            if (IsTMz) {
                Ceye_xn.CreateStruct(PMLCellxn, ny);
                Ceyhz_xn.CreateStruct(PMLCellxn, ny);

                Chzxh_xn.CreateStruct(Hzx_xn);
                Chzxey_xn.CreateStruct(Hzx_xn);

                Chzyh_xn.CreateStruct(Hzy_xn);
                Chzyex_xn.CreateStruct(Hzy_xn);

                //Chyh_xn.CreateStruct(PMLCellxn, ny);
                //Chyez_xn.CreateStruct(PMLCellxn, ny);
                //Cezxe_xn.CreateStruct(Hzx_xn);
                //Cezxhy_xn.CreateStruct(Hzx_xn);
                //Cezye_xn.CreateStruct(Hzy_xn);
                //Cezyhx_xn.CreateStruct(Hzy_xn);
                //Ceye_xn.CreateStruct(PMLCellxn, ny);
                //Ceyhz_xn.CreateStruct(PMLCellxn, ny);
                //Chzxh_xn.CreateStruct(Hzx_xn);
                //Chzxey_xn.CreateStruct(Hzx_xn);
                //Chzyh_xn.CreateStruct(Hzy_xn);
                //Chzyex_xn.CreateStruct(Hzy_xn);
            }
        }
        if (IsPMLxp) {
            if (IsTEz) {
                Chyh_xp.CreateStruct(PMLCellxp, nym1);
                Chyez_xp.CreateStruct(PMLCellxp, nym1);
                Cezxe_xp.CreateStruct(Ezx_xp);
                Cezxhy_xp.CreateStruct(Ezx_xp);
                Cezye_xp.CreateStruct(Ezy_xp);
                Cezyhx_xp.CreateStruct(Ezy_xp);
            }
            if (IsTMz) {
                Ceye_xp.CreateStruct(PMLCellxp, ny);
                Ceyhz_xp.CreateStruct(PMLCellxp, ny);

                Chzxh_xp.CreateStruct(Hzx_xp);
                Chzxey_xp.CreateStruct(Hzx_xp);

                Chzyh_xp.CreateStruct(Hzy_xp);
                Chzyex_xp.CreateStruct(Hzy_xp);

                //Chyh_xp.CreateStruct(PMLCellxp, ny);
                //Chyez_xp.CreateStruct(PMLCellxp, ny);
                //Cezxe_xp.CreateStruct(Hzx_xp);
                //Cezxhy_xp.CreateStruct(Hzx_xp);
                //Cezye_xp.CreateStruct(Hzy_xp);
                //Cezyhx_xp.CreateStruct(Hzy_xp);
                //Ceye_xp.CreateStruct(PMLCellxp, ny);
                //Ceyhz_xp.CreateStruct(PMLCellxp, ny);
                //Chzxh_xp.CreateStruct(Hzx_xp);
                //Chzxey_xp.CreateStruct(Hzx_xp);
                //Chzyh_xp.CreateStruct(Hzy_xp);
                //Chzyex_xp.CreateStruct(Hzy_xp);
            }
        }
        if (IsPMLyn) {
            if (IsTEz) {
                Chxh_yn.CreateStruct(nxm1, PMLCellyn);
                Chxez_yn.CreateStruct(nxm1, PMLCellyn);
                Cezye_yn.CreateStruct(Ezy_yn);
                Cezyhx_yn.CreateStruct(Ezy_yn);
                Cezxe_yn.CreateStruct(Ezx_yn);
                Cezxhy_yn.CreateStruct(Ezx_yn);
            }
            if (IsTMz) {
                Cexe_yn.CreateStruct(nx, PMLCellyn, 0.0);
                Cexhz_yn.CreateStruct(nx, PMLCellyn, 0.0);

                Chzxh_yn.CreateStruct(Hzx_yn);
                Chzxey_yn.CreateStruct(Hzx_yn);

                Chzyh_yn.CreateStruct(Hzy_yn, 0.0);
                Chzyex_yn.CreateStruct(Hzy_yn, 0.0);

                //Chxh_yn.CreateStruct(nx, PMLCellyn);
                //Chxez_yn.CreateStruct(nx, PMLCellyn);
                //Cezye_yn.CreateStruct(Hzy_yn);
                //Cezyhx_yn.CreateStruct(Hzy_yn);
                //Cezxe_yn.CreateStruct(Hzx_yn);
                //Cezxhy_yn.CreateStruct(Hzx_yn);
                //Cexe_yn.CreateStruct(nx, PMLCellyn);
                //Cexhz_yn.CreateStruct(nx, PMLCellyn);
                //Chzxh_yn.CreateStruct(Hzx_yn);
                //Chzxey_yn.CreateStruct(Hzx_yn);
                //Chzyh_yn.CreateStruct(Hzy_yn);
                //Chzyex_yn.CreateStruct(Hzy_yn);
            }
        }
        if (IsPMLyp) {
            if (IsTEz) {
                Chxh_yp.CreateStruct(nxm1, PMLCellyp);
                Chxez_yp.CreateStruct(nxm1, PMLCellyp);
                Cezye_yp.CreateStruct(Ezy_yp);
                Cezyhx_yp.CreateStruct(Ezy_yp);
                Cezxe_yp.CreateStruct(Ezx_yp);
                Cezxhy_yp.CreateStruct(Ezx_yp);
            }
            if (IsTMz) {
                Cexe_yp.CreateStruct(nx, PMLCellyp, 0.0);
                Cexhz_yp.CreateStruct(nx, PMLCellyp, 0.0);

                Chzxh_yp.CreateStruct(Hzx_yp);
                Chzxey_yp.CreateStruct(Hzx_yp);

                Chzyh_yp.CreateStruct(Hzy_yp, 0.0);
                Chzyex_yp.CreateStruct(Hzy_yp, 0.0);

                //Chxh_yp.CreateStruct(nx, PMLCellyp);
                //Chxez_yp.CreateStruct(nx, PMLCellyp);
                //Cezye_yp.CreateStruct(Hzy_yp);
                //Cezyhx_yp.CreateStruct(Hzy_yp);
                //Cezxe_yp.CreateStruct(Hzx_yp);
                //Cezxhy_yp.CreateStruct(Hzx_yp);
                //Cexe_yp.CreateStruct(nx, PMLCellyp);
                //Cexhz_yp.CreateStruct(nx, PMLCellyp);
                //Chzxh_yp.CreateStruct(Hzx_yp);
                //Chzxey_yp.CreateStruct(Hzx_yp);
                //Chzyh_yp.CreateStruct(Hzy_yp);
                //Chzyex_yp.CreateStruct(Hzy_yp);
            }
        }
    }
    return true;
}
//void InitPMLBoundTMx(){
//
//	unsigned int i,j;
//
//	MyDataF	*rho_e;
//	MyDataF	*rho_m;
//
//	MyStruct	sigma_pex_xn;
//	MyStruct	sigma_pmx_xn;
//
//	MyStruct	sigma_pey_yn;
//	MyStruct	sigma_pmy_yn;
//
//	MyStruct	sigma_pex_xp;
//	MyStruct	sigma_pmx_xp;
//
//	MyStruct	sigma_pey_yp;
//	MyStruct	sigma_pmy_yp;
//	
//
//	MyDataF	sigma_max;
//
//	printf("Initializing pml boundary conditions 2d TMx...\n");
//	
//
//
///*=====================================================================================================================*/
//	if(IsPMLxn){
//		rho_e = (MyDataF*)malloc(PMLCellxn*sizeof(MyDataF));
//		rho_m = (MyDataF*)malloc(PMLCellxn*sizeof(MyDataF));		
//
//		sigma_pex_xn.CreateStruct(PMLCellxn,ny,0.0);		
//		sigma_pmx_xn.CreateStruct(PMLCellxn,ny,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);
//
//		for(i=PMLCellxn;i>0;--i){
//			rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
//			rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
//		}
//
//		for(i=0;i<PMLCellxn;i++)
//			for(j=0;j<ny;++j){
//				sigma_pex_xn.data[i][j] = sigma_max*pow(rho_e[i],pml_order);
//				sigma_pmx_xn.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
//			}
//
//
//		for(i=0;i<PMLCellxn;++i)
//			for(j=0;j<ny;j++){
//				/* Coefficients updating Ey */
//				Ceye_xn.data[i][j] = (2*eps_0-dt*sigma_pex_xn.data[i][j])/(2*eps_0+dt*sigma_pex_xn.data[i][j]);
//				Ceyhz_xn.data[i][j] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xn.data[i][j]);
//			
//				/* Coefficients updating Hzx */
//				Chzxh_xn.data[i][j] = (2*mu_0-dt*sigma_pmx_xn.data[i][j])/(2*mu_0+dt*sigma_pmx_xn.data[i][j]);
//				Chzxey_xn.data[i][j] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xn.data[i][j]);
//
//			}
//		
//		free(rho_e);
//		free(rho_m);
//	}
//
///*====================================================================================================================*/
//	if(IsPMLxp){
//
//		rho_e = (MyDataF*)malloc(PMLCellxp*sizeof(MyDataF));
//		rho_m = (MyDataF*)malloc(PMLCellxp*sizeof(MyDataF));
//
//		sigma_pex_xp.CreateStruct(PMLCellxp,ny,0.0);		
//		sigma_pmx_xp.CreateStruct(PMLCellxp,ny,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);
//
//		for(i=0;i<PMLCellxp;++i){
//			rho_e[i] = (i+1-0.75)/PMLCellxp;
//			rho_m[i] = (i+1-0.25)/PMLCellxp;
//		}
//
//		for(i=0;i<PMLCellxp;i++)
//			for(j=0;j<ny;++j){
//				sigma_pex_xp.data[i][j] = sigma_max*pow(rho_e[i],pml_order);
//				sigma_pmx_xp.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
//			}
//
//
//		for(i=0;i<PMLCellxp;++i)
//			for(j=0;j<ny;j++){
//				/* Coefficients updating Ey */
//				Ceye_xp.data[i][j] = (2*eps_0-dt*sigma_pex_xp.data[i][j])/(2*eps_0+dt*sigma_pex_xp.data[i][j]);
//				Ceyhz_xp.data[i][j] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xp.data[i][j]);
//				
//				/* Coefficients updating Hzx */
//				Chzxh_xp.data[i][j] = (2*mu_0-dt*sigma_pmx_xp.data[i][j])/(2*mu_0+dt*sigma_pmx_xp.data[i][j]);
//				Chzxey_xp.data[i][j] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp.data[i][j]);
//
//			
//			}
//			
//		free(rho_e);
//		free(rho_m);
//		
//	}
///*=====================================================================================================================*/
//	if(IsPMLyn){
//		rho_e = (MyDataF*)malloc(PMLCellyn*sizeof(MyDataF));
//		rho_m = (MyDataF*)malloc(PMLCellyn*sizeof(MyDataF));
//
//		sigma_pey_yn.CreateStruct(nx,PMLCellyn,0.0);		
//		sigma_pmy_yn.CreateStruct(nx,PMLCellyn,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);
//
//		for(i=PMLCellyn;i>0;--i){
//			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
//			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
//		}
//
//		for(i=0;i<nx;i++)
//			for(j=0;j<PMLCellyn;++j){
//				sigma_pey_yn.data[i][j] = sigma_max*pow(rho_e[j],pml_order);
//				sigma_pmy_yn.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
//			}
//
//
//		for(i=0;i<nx;++i)
//			for(j=0;j<PMLCellyn;j++){
//				/* Coefficients updating Ey */
//				Cexe_yn.data[i][j] = (2*eps_0-dt*sigma_pey_yn.data[i][j])/(2*eps_0+dt*sigma_pey_yn.data[i][j]);
//				Cexhz_yn.data[i][j] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yn.data[i][j]);
//				
//				/* Coefficients updating Hzy */
//				Chzyh_yn.data[i][j] = (2*mu_0-dt*sigma_pmy_yn.data[i][j])/(2*mu_0+dt*sigma_pmy_yn.data[i][j]);
//				Chzyex_yn.data[i][j] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn.data[i][j]);
//			
//			}
//		
//		free(rho_e);
//		free(rho_m);
//	}
///*=====================================================================================================================*/
//	if(IsPMLyp){
//
//		rho_e = (MyDataF*)malloc(PMLCellyp*sizeof(MyDataF));
//		rho_m = (MyDataF*)malloc(PMLCellyp*sizeof(MyDataF));
//
//		sigma_pey_yp.CreateStruct(nx,PMLCellyp,0.0);
//		sigma_pmy_yp.CreateStruct(nx,PMLCellyp,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);
//
//		for(i=0;i<PMLCellyp;++i){
//			rho_e[i] = (i+1-0.75)/PMLCellyp;
//			rho_m[i] = (i+1-0.25)/PMLCellyp;
//		}
//
//		for(i=0;i<nx;i++)
//			for(j=0;j<PMLCellyp;++j){
//				sigma_pey_yp.data[i][j] = sigma_max*pow(rho_e[j],pml_order);
//				sigma_pmy_yp.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
//			}
//
//		for(i=0;i<nx;++i)
//			for(j=0;j<PMLCellyp;j++){
//				/* Coefficients updating Ey */
//				Cexe_yp.data[i][j] = (2*eps_0-dt*sigma_pey_yp.data[i][j])/(2*eps_0+dt*sigma_pey_yp.data[i][j]);
//				Cexhz_yp.data[i][j] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yp.data[i][j]);
//								
//				/* Coefficients updating Hzy */
//				Chzyh_yp.data[i][j] = (2*mu_0-dt*sigma_pmy_yp.data[i][j])/(2*mu_0+dt*sigma_pmy_yp.data[i][j]);
//				Chzyex_yp.data[i][j] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yp.data[i][j]);
//			
//			}
//	free(rho_e);
//	free(rho_m);
//	}
//	
//	printf("End of initializing pml boundary conditions 2d TMx...\n");
//}
//
//

/**
 * @brief initialize the coefficients for the PML arrays
 */
void InitPMLCoefficients() {
    double *rho_e;
    double *rho_m;
    double sigma_max;
    unsigned int i, j;
    MyStruct sigma_pex_xn, sigma_pmx_xn, sigma_pex_xp, sigma_pmx_xp;
    MyStruct sigma_pey_yn, sigma_pmy_yn, sigma_pey_yp, sigma_pmy_yp;


    if (IsTEz) {
        //		if (IsPMLxn) {
        //        rho_e = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));
        //        rho_m = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));
        //
        //        sigma_pex_xn.CreateStruct(PMLCellxn, ny);
        //        sigma_pmx_xn.CreateStruct(PMLCellxn, ny);
        //
        //        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxn);
        //
        //        for (i = PMLCellxn; i > 0; --i) {
        //            rho_e[PMLCellxn - i] = (i - 0.75) / PMLCellxn;
        //            rho_m[PMLCellxn - i] = (i - 0.25) / PMLCellxn;
        //        }
        //
        //        for (i = 0; i < PMLCellxn; i++)
        //            for (j = 0; j < ny; ++j) {
        //                sigma_pex_xn.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
        //                sigma_pmx_xn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
        //            }
        //        Chzyh_xn.InitStructData(1.0);
        //        Chzyex_xn.InitStructData(dt / (dy * mu_0));
        //
        //
        //        for (i = 0; i < PMLCellxn; ++i)
        //            for (j = 0; j < ny; j++) {
        //                /* Coefficients updating Ey */
        //                Ceye_xn.data[i][j] = (2 * eps_0 - dt * sigma_pex_xn.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);
        //                Ceyhz_xn.data[i][j] = -(2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);
        //
        //                /* Coefficients updating Hzx */
        //                Chzxh_xn.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xn.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);
        //                Chzxey_xn.data[i][j] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);
        //            }
        //
        //        free(rho_e);
        //        free(rho_m);
        //    }
        //
        //    /*====================================================================================================================*/
        //    if (IsPMLxp) {
        //
        //        rho_e = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));
        //        rho_m = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));
        //
        //        sigma_pex_xp.CreateStruct(PMLCellxp, ny);
        //        sigma_pmx_xp.CreateStruct(PMLCellxp, ny);
        //
        //        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxp);
        //
        //        for (i = 0; i < PMLCellxp; ++i) {
        //            rho_e[i] = (i + 1 - 0.75) / PMLCellxp;
        //            rho_m[i] = (i + 1 - 0.25) / PMLCellxp;
        //        }
        //
        //        for (i = 0; i < PMLCellxp; i++)
        //            for (j = 0; j < ny; ++j) {
        //                sigma_pex_xp.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
        //                sigma_pmx_xp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
        //            }
        //
        //        Chzyh_xp.InitStructData(1.0);
        //        Chzyex_xp.InitStructData(dt / (dy * mu_0));
        //
        //
        //        for (i = 0; i < PMLCellxp; ++i)
        //            for (j = 0; j < ny; j++) {
        //                /* Coefficients updating Ey */
        //                Ceye_xp.data[i][j] = (2 * eps_0 - dt * sigma_pex_xp.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);
        //                Ceyhz_xp.data[i][j] = -(2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);
        //
        //                /* Coefficients updating Hzx */
        //                Chzxh_xp.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xp.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);
        //                Chzxey_xp.data[i][j] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);
        //
        //
        //            }
        //
        //        //sigma_pmx_xp.~MyStruct();
        //        //sigma_pex_xp.~MyStruct();
        //
        //        free(rho_e);
        //        free(rho_m);
        //
        //    }
        //    /*=====================================================================================================================*/
        //    if (IsPMLyn) {
        //        rho_e = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));
        //        rho_m = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));
        //
        //        sigma_pey_yn.CreateStruct(nx, PMLCellyn);
        //        sigma_pmy_yn.CreateStruct(nx, PMLCellyn);
        //
        //        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyn);
        //
        //        for (i = PMLCellyn; i > 0; --i) {
        //            rho_e[PMLCellyn - i] = (i - 0.75) / PMLCellyn;
        //            rho_m[PMLCellyn - i] = (i - 0.25) / PMLCellyn;
        //        }
        //
        //        for (i = 0; i < nx; i++)
        //            for (j = 0; j < PMLCellyn; ++j) {
        //                sigma_pey_yn.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
        //                sigma_pmy_yn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
        //            }
        //
        //        Chzxh_yn.InitStructData(1.0);
        //        Chzxey_yn.InitStructData(-dt / (dx * mu_0));
        //
        //
        //        for (i = 0; i < nx; ++i)
        //            for (j = 0; j < PMLCellyn; j++) {
        //                /* Coefficients updating Ey */
        //                Cexe_yn.data[i][j] = (2 * eps_0 - dt * sigma_pey_yn.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);
        //                Cexhz_yn.data[i][j] = (2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);
        //
        //                /* Coefficients updating Hzy */
        //                Chzyh_yn.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yn.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);
        //                Chzyex_yn.data[i][j] = (2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);
        //
        //            }
        //
        //        //sigma_pmy_yn.~MyStruct();
        //        //sigma_pey_yn.~MyStruct();
        //
        //        free(rho_e);
        //        free(rho_m);
        //    }
        //    /*=====================================================================================================================*/
        //    if (IsPMLyp) {
        //
        //        rho_e = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));
        //        rho_m = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));
        //
        //        sigma_pey_yp.CreateStruct(nx, PMLCellyp);
        //        sigma_pmy_yp.CreateStruct(nx, PMLCellyp);
        //
        //        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyp);
        //
        //        for (i = 0; i < PMLCellyp; ++i) {
        //            rho_e[i] = (i + 1 - 0.75) / PMLCellyp;
        //            rho_m[i] = (i + 1 - 0.25) / PMLCellyp;
        //        }
        //
        //        for (i = 0; i < nx; i++)
        //            for (j = 0; j < PMLCellyp; ++j) {
        //                sigma_pey_yp.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
        //                sigma_pmy_yp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
        //            }
        //
        //        Chzxh_yp.InitStructData(1.0);
        //        Chzxey_yp.InitStructData(-dt / (dx * mu_0));
        //
        //
        //        for (i = 0; i < nx; ++i)
        //            for (j = 0; j < PMLCellyp; j++) {
        //                /* Coefficients updating Ey */
        //                Cexe_yp.data[i][j] = (2 * eps_0 - dt * sigma_pey_yp.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);
        //                Cexhz_yp.data[i][j] = (2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);
        //
        //                /* Coefficients updating Hzy */
        //                Chzyh_yp.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yp.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);
        //                Chzyex_yp.data[i][j] = (2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);
        //
        //            }
        //
        //        //sigma_pmy_yp.~MyStruct();
        //        //sigma_pey_yp.~MyStruct();
        //
        //        free(rho_e);
        //        free(rho_m);
        //    }
        ////====================================================================================================
        //		if(IsPMLxn){
        //		rho_e = (MyDataF*)malloc(PMLCellxn*sizeof(MyDataF));
        //		rho_m = (MyDataF*)malloc(PMLCellxn*sizeof(MyDataF));		
        //
        //		sigma_pex_xn.CreateStruct(PMLCellxn,ny,0.0);		
        //		sigma_pmx_xn.CreateStruct(PMLCellxn,ny,0.0);
        //		
        //		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);
        //
        //		for(i=PMLCellxn;i>0;--i){
        //			rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
        //			rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
        //		}
        //
        //		for(i=0;i<PMLCellxn;i++)
        //			for(j=0;j<ny;++j){
        //				sigma_pex_xn.data[i][j] = sigma_max*pow(rho_e[i],pml_order);
        //				sigma_pmx_xn.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
        //			}
        //
        //
        //		for(i=0;i<PMLCellxn;++i)
        //			for(j=0;j<ny;j++){
        //				/* Coefficients updating Ey */
        //				Ceye_xn.data[i][j] = (2*eps_0-dt*sigma_pex_xn.data[i][j])/(2*eps_0+dt*sigma_pex_xn.data[i][j]);
        //				Ceyhz_xn.data[i][j] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xn.data[i][j]);
        //			
        //				/* Coefficients updating Hzx */
        //				Chzxh_xn.data[i][j] = (2*mu_0-dt*sigma_pmx_xn.data[i][j])/(2*mu_0+dt*sigma_pmx_xn.data[i][j]);
        //				Chzxey_xn.data[i][j] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xn.data[i][j]);
        //
        //			}
        //		
        //		free(rho_e);
        //		free(rho_m);
        //	}
        //
        ///*====================================================================================================================*/
        //	if(IsPMLxp){
        //
        //		rho_e = (MyDataF*)malloc(PMLCellxp*sizeof(MyDataF));
        //		rho_m = (MyDataF*)malloc(PMLCellxp*sizeof(MyDataF));
        //
        //		sigma_pex_xp.CreateStruct(PMLCellxp,ny,0.0);		
        //		sigma_pmx_xp.CreateStruct(PMLCellxp,ny,0.0);
        //		
        //		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);
        //
        //		for(i=0;i<PMLCellxp;++i){
        //			rho_e[i] = (i+1-0.75)/PMLCellxp;
        //			rho_m[i] = (i+1-0.25)/PMLCellxp;
        //		}
        //
        //		for(i=0;i<PMLCellxp;i++)
        //			for(j=0;j<ny;++j){
        //				sigma_pex_xp.data[i][j] = sigma_max*pow(rho_e[i],pml_order);
        //				sigma_pmx_xp.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
        //			}
        //
        //
        //		for(i=0;i<PMLCellxp;++i)
        //			for(j=0;j<ny;j++){
        //				/* Coefficients updating Ey */
        //				Ceye_xp.data[i][j] = (2*eps_0-dt*sigma_pex_xp.data[i][j])/(2*eps_0+dt*sigma_pex_xp.data[i][j]);
        //				Ceyhz_xp.data[i][j] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xp.data[i][j]);
        //				
        //				/* Coefficients updating Hzx */
        //				Chzxh_xp.data[i][j] = (2*mu_0-dt*sigma_pmx_xp.data[i][j])/(2*mu_0+dt*sigma_pmx_xp.data[i][j]);
        //				Chzxey_xp.data[i][j] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp.data[i][j]);
        //
        //			
        //			}
        //			
        //		free(rho_e);
        //		free(rho_m);
        //		
        //	}
        ///*=====================================================================================================================*/
        //	if(IsPMLyn){
        //		rho_e = (MyDataF*)malloc(PMLCellyn*sizeof(MyDataF));
        //		rho_m = (MyDataF*)malloc(PMLCellyn*sizeof(MyDataF));
        //
        //		sigma_pey_yn.CreateStruct(nx,PMLCellyn,0.0);		
        //		sigma_pmy_yn.CreateStruct(nx,PMLCellyn,0.0);
        //		
        //		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);
        //
        //		for(i=PMLCellyn;i>0;--i){
        //			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
        //			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
        //		}
        //
        //		for(i=0;i<nx;i++)
        //			for(j=0;j<PMLCellyn;++j){
        //				sigma_pey_yn.data[i][j] = sigma_max*pow(rho_e[j],pml_order);
        //				sigma_pmy_yn.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
        //			}
        //
        //
        //		for(i=0;i<nx;++i)
        //			for(j=0;j<PMLCellyn;j++){
        //				/* Coefficients updating Ey */
        //				Cexe_yn.data[i][j] = (2*eps_0-dt*sigma_pey_yn.data[i][j])/(2*eps_0+dt*sigma_pey_yn.data[i][j]);
        //				Cexhz_yn.data[i][j] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yn.data[i][j]);
        //				
        //				/* Coefficients updating Hzy */
        //				Chzyh_yn.data[i][j] = (2*mu_0-dt*sigma_pmy_yn.data[i][j])/(2*mu_0+dt*sigma_pmy_yn.data[i][j]);
        //				Chzyex_yn.data[i][j] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn.data[i][j]);
        //			
        //			}
        //		
        //		free(rho_e);
        //		free(rho_m);
        //	}
        ///*=====================================================================================================================*/
        //	if(IsPMLyp){
        //
        //		rho_e = (MyDataF*)malloc(PMLCellyp*sizeof(MyDataF));
        //		rho_m = (MyDataF*)malloc(PMLCellyp*sizeof(MyDataF));
        //
        //		sigma_pey_yp.CreateStruct(nx,PMLCellyp,0.0);
        //		sigma_pmy_yp.CreateStruct(nx,PMLCellyp,0.0);
        //		
        //		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);
        //
        //		for(i=0;i<PMLCellyp;++i){
        //			rho_e[i] = (i+1-0.75)/PMLCellyp;
        //			rho_m[i] = (i+1-0.25)/PMLCellyp;
        //		}
        //
        //		for(i=0;i<nx;i++)
        //			for(j=0;j<PMLCellyp;++j){
        //				sigma_pey_yp.data[i][j] = sigma_max*pow(rho_e[j],pml_order);
        //				sigma_pmy_yp.data[i][j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
        //			}
        //
        //		for(i=0;i<nx;++i)
        //			for(j=0;j<PMLCellyp;j++){
        //				/* Coefficients updating Ey */
        //				Cexe_yp.data[i][j] = (2*eps_0-dt*sigma_pey_yp.data[i][j])/(2*eps_0+dt*sigma_pey_yp.data[i][j]);
        //				Cexhz_yp.data[i][j] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yp.data[i][j]);
        //								
        //				/* Coefficients updating Hzy */
        //				Chzyh_yp.data[i][j] = (2*mu_0-dt*sigma_pmy_yp.data[i][j])/(2*mu_0+dt*sigma_pmy_yp.data[i][j]);
        //				Chzyex_yp.data[i][j] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yp.data[i][j]);
        //			
        //			}
        //	free(rho_e);
        //	free(rho_m);
        //	}
        if (IsPMLxn) {

            rho_e = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));

            sigma_pex_xn.CreateStruct(PMLCellxn, nym1, 0.0);
            sigma_pmx_xn.CreateStruct(PMLCellxn, nym1, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxn);

            for (i = PMLCellxn; i > 0; --i) {
                rho_e[PMLCellxn - i] = (i - 0.75) / PMLCellxn;
                rho_m[PMLCellxn - i] = (i - 0.25) / PMLCellxn;
            }

            for (i = 0; i < PMLCellxn; i++)
                for (j = 0; j < nym1; ++j) {
                    sigma_pex_xn.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
                    sigma_pmx_xn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
                }


            Cezye_xn.InitStructData(1.0);
            Cezyhx_xn.InitStructData(-dt / (dy * eps_0));

            for (i = 0; i < PMLCellxn; ++i)
                for (j = 0; j < nym1; j++) {
                    /* Coefficients updating Hx */
                    Chyh_xn.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xn.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);
                    Chyez_xn.data[i][j] = (2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);

                    /* Coefficients updating Ezx */
                    Cezxe_xn.data[i][j] = (2 * eps_0 - dt * sigma_pex_xn.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);
                    Cezxhy_xn.data[i][j] = (2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);

                }

            free(rho_e);
            free(rho_m);
        }

        if (IsPMLyn) {

            rho_e = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));

            sigma_pey_yn.CreateStruct(nxm1, PMLCellyn, 0.0);
            sigma_pmy_yn.CreateStruct(nxm1, PMLCellyn, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyn);

            for (i = PMLCellyn; i > 0; --i) {
                rho_e[PMLCellyn - i] = (i - 0.75) / PMLCellyn;
                rho_m[PMLCellyn - i] = (i - 0.25) / PMLCellyn;
            }

            for (j = 0; j < PMLCellyn; j++)
                for (i = 0; i < nxm1; ++i) {
                    sigma_pey_yn.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
                    sigma_pmy_yn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
                }


            for (i = 0; i < nxm1; ++i)
                for (j = 0; j < PMLCellyn; j++) {
                    /* Coefficients updating Hx */
                    Chxh_yn.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yn.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);
                    Chxez_yn.data[i][j] = -(2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);

                    /* Coefficients updating Ezy */
                    Cezye_yn.data[i][j] = (2 * eps_0 - dt * sigma_pey_yn.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);
                    Cezyhx_yn.data[i][j] = -(2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);

                }


            Cezxe_yn.InitStructData(1.0);
            Cezxhy_yn.InitStructData(dt / (dx * eps_0));

            free(rho_e);
            free(rho_m);
        }

        if (IsPMLxp) {

            rho_e = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));

            sigma_pex_xp.CreateStruct(PMLCellxp, nym1, 0.0);
            sigma_pmx_xp.CreateStruct(PMLCellxp, nym1, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxp);

            for (i = 0; i < PMLCellyp; ++i) {
                rho_e[i] = (i + 1 - 0.75) / PMLCellxp;
                rho_m[i] = (i + 1 - 0.25) / PMLCellxp;
            }

            for (j = 0; j < nym1; j++)
                for (i = 0; i < PMLCellxp; ++i) {
                    sigma_pex_xp.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
                    sigma_pmx_xp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
                }
            for (i = 0; i < PMLCellxp; ++i)
                for (j = 0; j < nym1; j++) {
                    /* Coefficients updating Hy */
                    Chyh_xp.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xp.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);
                    Chyez_xp.data[i][j] = (2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);

                    /* Coefficients updating Ezx */
                    Cezxe_xp.data[i][j] = (2 * eps_0 - dt * sigma_pex_xp.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);
                    Cezxhy_xp.data[i][j] = (2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);

                }
            Cezye_xp.InitStructData(1.0);
            Cezyhx_xp.InitStructData(-dt / (dy * eps_0));

            free(rho_e);
            free(rho_m);
        }
        if (IsPMLyp) {

            rho_e = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));

            sigma_pey_yp.CreateStruct(nxm1, PMLCellyp, 0.0);
            sigma_pmy_yp.CreateStruct(nxm1, PMLCellyp, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyp);

            for (i = 0; i < PMLCellyp; ++i) {
                rho_e[i] = (i + 1 - 0.75) / PMLCellyp;
                rho_m[i] = (i + 1 - 0.25) / PMLCellyp;
            }

            for (i = 0; i < nxm1; i++)
                for (j = 0; j < PMLCellyp; ++j) {
                    sigma_pey_yp.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
                    sigma_pmy_yp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
                }

            for (i = 0; i < nxm1; ++i)
                for (j = 0; j < PMLCellyp; j++) {
                    /* Coefficients updating Hy */
                    Chxh_yp.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yp.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);
                    Chxez_yp.data[i][j] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);

                    /* Coefficients updating Ezy */
                    Cezye_yp.data[i][j] = (2 * eps_0 - dt * sigma_pey_yp.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);
                    Cezyhx_yp.data[i][j] = -(2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);

                }

            Cezxe_yp.InitStructData(1.0);
            Cezxhy_yp.InitStructData(dt / (dx * eps_0));

            free(rho_e);
            free(rho_m);
        }


    }
    if (IsTMz) {
        if (IsPMLxn) {
            rho_e = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));

            sigma_pex_xn.CreateStruct(PMLCellxn, ny, 0.0);
            sigma_pmx_xn.CreateStruct(PMLCellxn, ny, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxn);

            for (i = PMLCellxn; i > 0; --i) {
                rho_e[PMLCellxn - i] = (i - 0.75) / PMLCellxn;
                rho_m[PMLCellxn - i] = (i - 0.25) / PMLCellxn;
            }

            for (i = 0; i < PMLCellxn; i++)
                for (j = 0; j < ny; ++j) {
                    sigma_pex_xn.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
                    sigma_pmx_xn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
                }


            for (i = 0; i < PMLCellxn; ++i)
                for (j = 0; j < ny; j++) {
                    /* Coefficients updating Ey */
                    Ceye_xn.data[i][j] = (2 * eps_0 - dt * sigma_pex_xn.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);
                    Ceyhz_xn.data[i][j] = -(2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xn.data[i][j]);

                    /* Coefficients updating Hzx */
                    Chzxh_xn.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xn.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);
                    Chzxey_xn.data[i][j] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xn.data[i][j]);

                }
            Chzyh_xn.InitStructData(1.0);
            Chzyex_xn.InitStructData(dt / (dy * mu_0));

            free(rho_e);
            free(rho_m);
        }

        /*====================================================================================================================*/
        if (IsPMLxp) {

            rho_e = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));

            sigma_pex_xp.CreateStruct(PMLCellxp, ny, 0.0);
            sigma_pmx_xp.CreateStruct(PMLCellxp, ny, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxp);

            for (i = 0; i < PMLCellxp; ++i) {
                rho_e[i] = (i + 1 - 0.75) / PMLCellxp;
                rho_m[i] = (i + 1 - 0.25) / PMLCellxp;
            }

            for (i = 0; i < PMLCellxp; i++)
                for (j = 0; j < ny; ++j) {
                    sigma_pex_xp.data[i][j] = sigma_max * pow(rho_e[i], pml_order);
                    sigma_pmx_xp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
                }


            for (i = 0; i < PMLCellxp; ++i)
                for (j = 0; j < ny; j++) {
                    /* Coefficients updating Ey */
                    Ceye_xp.data[i][j] = (2 * eps_0 - dt * sigma_pex_xp.data[i][j]) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);
                    Ceyhz_xp.data[i][j] = -(2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xp.data[i][j]);

                    /* Coefficients updating Hzx */
                    Chzxh_xp.data[i][j] = (2 * mu_0 - dt * sigma_pmx_xp.data[i][j]) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);
                    Chzxey_xp.data[i][j] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xp.data[i][j]);


                }
            Chzyh_xp.InitStructData(1.0);
            Chzyex_xp.InitStructData(dt / (dy * mu_0));

            free(rho_e);
            free(rho_m);

        }
        /*=====================================================================================================================*/
        if (IsPMLyn) {
            rho_e = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));

            sigma_pey_yn.CreateStruct(nx, PMLCellyn, 0.0);
            sigma_pmy_yn.CreateStruct(nx, PMLCellyn, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyn);

            for (i = PMLCellyn; i > 0; --i) {
                rho_e[PMLCellyn - i] = (i - 0.75) / PMLCellyn;
                rho_m[PMLCellyn - i] = (i - 0.25) / PMLCellyn;
            }

            for (i = 0; i < nx; i++)
                for (j = 0; j < PMLCellyn; ++j) {
                    sigma_pey_yn.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
                    sigma_pmy_yn.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
                }


            for (i = 0; i < nx; ++i)
                for (j = 0; j < PMLCellyn; j++) {
                    /* Coefficients updating Ey */
                    Cexe_yn.data[i][j] = (2 * eps_0 - dt * sigma_pey_yn.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);
                    Cexhz_yn.data[i][j] = (2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yn.data[i][j]);

                    /* Coefficients updating Hzy */
                    Chzyh_yn.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yn.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);
                    Chzyex_yn.data[i][j] = (2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yn.data[i][j]);

                }
            Chzxh_yn.InitStructData(1.0);
            Chzxey_yn.InitStructData(-dt / (dx * mu_0));

            free(rho_e);
            free(rho_m);
        }
        /*=====================================================================================================================*/
        if (IsPMLyp) {

            rho_e = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));
            rho_m = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));

            sigma_pey_yp.CreateStruct(nx, PMLCellyp, 0.0);
            sigma_pmy_yp.CreateStruct(nx, PMLCellyp, 0.0);

            sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyp);

            for (i = 0; i < PMLCellyp; ++i) {
                rho_e[i] = (i + 1 - 0.75) / PMLCellyp;
                rho_m[i] = (i + 1 - 0.25) / PMLCellyp;
            }

            for (i = 0; i < nx; i++)
                for (j = 0; j < PMLCellyp; ++j) {
                    sigma_pey_yp.data[i][j] = sigma_max * pow(rho_e[j], pml_order);
                    sigma_pmy_yp.data[i][j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
                }

            for (i = 0; i < nx; ++i)
                for (j = 0; j < PMLCellyp; j++) {
                    /* Coefficients updating Ey */
                    Cexe_yp.data[i][j] = (2 * eps_0 - dt * sigma_pey_yp.data[i][j]) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);
                    Cexhz_yp.data[i][j] = (2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yp.data[i][j]);

                    /* Coefficients updating Hzy */
                    Chzyh_yp.data[i][j] = (2 * mu_0 - dt * sigma_pmy_yp.data[i][j]) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);
                    Chzyex_yp.data[i][j] = (2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yp.data[i][j]);

                }
            Chzxh_yp.InitStructData(1.0);
            Chzxey_yp.InitStructData(-dt / (dx * mu_0));

            free(rho_e);
            free(rho_m);
        }
    }
    cout << "End of initializing pml boundary conditions 2d TMx..." << endl;
}

/**
 * @brief initial PML 
 * @param ncpml Number of cells in PML boundary
 * @param ncx size of domain in x direction
 * @param ncy size of domain in y direction
 */
void InitialPML(unsigned ncpml, unsigned ncx, unsigned ncy) {

    boundPMLdef(ncpml);
    CreatePMLArrays(ncx, ncy);
    CreatePMLCoeffients(ncx, ncy);
    InitPMLCoefficients();

}

void FreePMLSpace() {

}

void UpdateMFieldForPML(MyStruct &hx, MyStruct &hy, MyStruct &hz, const MyStruct &ex, const MyStruct &ey, const MyStruct &ez) {
    /* update magnetic fields at PML regions */

    if (IsTEz) {
        UpdMagFldForPML_TEz(hx, hy, ez);
    }
    if (IsTMz) {
        UpdMagFldForPML_TMz(hz, ex, ey);
    }
}

void UpdEltFldForPML_TEz(MyStruct &ez, const MyStruct &hx, const MyStruct &hy) {
    unsigned int i, j;
    unsigned int im, jm;
    unsigned km = 0;
    if (IsPMLxn) {
        for (i = 0; i < Ezx_xn.nx; ++i) {
            km = i + 1;
            for (j = 0, jm = j + 1; j < Ezx_xn.ny; j++, jm++) {
                Ezx_xn.data[i][j] = Cezxe_xn.data[i][j] * Ezx_xn.data[i][j]
                        + Cezxhy_xn.data[i][j]*(hy.data[km][jm] - hy.data[i][jm]);
            }
        }
        for (j = 0, jm = pjs; j < Ezy_xn.ny; j++, jm++) {
            km = jm + 1;
            for (i = 0, im = i + 1; i < Ezy_xn.nx; ++i, im++) {
                Ezy_xn.data[i][j] = Cezye_xn.data[i][j] * Ezy_xn.data[i][j]
                        + Cezyhx_xn.data[i][j]*(hx.data[im][km] - hx.data[im][jm]);
            }
        }
    }

    if (IsPMLxp) {
        for (i = 0, im = pie; i < Ezx_xp.nx; ++i, im++) {
            km = im - 1;
            for (j = 0, jm = j + 1; j < Ezx_xp.ny; j++, jm++) {
                Ezx_xp.data[i][j] = Cezxe_xp.data[i][j] * Ezx_xp.data[i][j]
                        + Cezxhy_xp.data[i][j]*(hy.data[im][jm] - hy.data[km][jm]);
            }
        }
        for (j = 0, jm = pjs; j < Ezy_xp.ny; j++, jm++) {
            km = jm + 1;
            for (i = 0, im = pie; i < Ezy_xp.nx; ++i, im++) {
                Ezy_xp.data[i][j] = Cezye_xp.data[i][j] * Ezy_xp.data[i][j]
                        + Cezyhx_xp.data[i][j]*(hx.data[im][km] - hx.data[im][jm]);
            }
        }
    }

    if (IsPMLyn) {
        for (i = 0, im = pis; i < Ezx_yn.nx; ++i, im++) {
            km = im + 1;
            for (j = 0, jm = j + 1; j < Ezx_yn.ny; j++, jm++) {
                Ezx_yn.data[i][j] = Cezxe_yn.data[i][j] * Ezx_yn.data[i][j]
                        + Cezxhy_yn.data[i][j]*(hy.data[km][jm] - hy.data[im][jm]);
            }
        }
        for (i = 0, im = i + 1; i < Ezy_yn.nx; ++i, im++)
            for (j = 0; j < Ezy_yn.ny; j++) {
                Ezy_yn.data[i][j] = Cezye_yn.data[i][j] * Ezy_yn.data[i][j]
                        + Cezyhx_yn.data[i][j]*(hx.data[im][j + 1] - hx.data[im][j]);

            }

    }
    if (IsPMLyp) {
        for (i = 0, im = pis; i < Ezx_yp.nx; ++i, im++) {
            km = im + 1;
            for (j = 0, jm = pje; j < Ezx_yp.ny; j++, jm++) {
                Ezx_yp.data[i][j] = Cezxe_yp.data[i][j] * Ezx_yp.data[i][j]
                        + Cezxhy_yp.data[i][j]*(hy.data[km][jm] - hy.data[im][jm]);
            }
        }
        for (j = 0, jm = pje; j < Ezy_yp.ny; j++, jm++) {
            km = jm - 1;
            for (i = 0, im = i + 1; i < Ezy_yp.nx; ++i, im++) {
                Ezy_yp.data[i][j] = Cezye_yp.data[i][j] * Ezy_yp.data[i][j]
                        + Cezyhx_yp.data[i][j]*(hx.data[im][jm] - hx.data[im][km]);
            }

        }
    }
    // 1
    for (i = 0; i < pis; i++) {
        km = i + 1;
        for (j = 0; j < pjs; j++) {
            ez.data[km][j + 1] = Ezx_xn.data[i][j] + Ezy_yn.data[i][j];
        }
    }
    // 2
    for (i = 0; i < pis; i++) {
        km = i + 1;
        for (j = 0, jm = pje; j < PMLCellyp; j++, jm++) {
            ez.data[km][jm] = Ezx_xn.data[i][jm - 1] + Ezy_yp.data[i][j];
        }
    }
    // 3
    for (i = pie; i < nx; i++) {
        im = i - pie;
        km = i - 1;
        for (j = pje; j < ny; j++) {
            ez.data[i][j] = Ezx_xp.data[im][j - 1] + Ezy_yp.data[km][j - pje];
        }
    }
    // 4
    for (i = pie; i < nx; i++) {
        im = i - pie;
        km = i - 1;
        for (j = 0; j < pjs; j++) {
            ez.data[i][j + 1] = Ezx_xp.data[im][j] + Ezy_yn.data[km][j];
        }
    }
    // 5
    jm = pie - 1;
    for (i = pis; i < jm; i++) {
        im = i + 1;
        km = i - pis;
        for (j = 0; j < pjs; j++) {
            ez.data[im][j + 1] = Ezx_yn.data[km][j] + Ezy_yn.data[i][j];
        }
    }
    // 6
    jm = pie - 1;
    for (i = pis; i < jm; i++) {
        im = i + 1;
        km = i - pis;
        for (j = 0; j < PMLCellyp; j++) {
            ez.data[im][j + pje] = Ezx_yp.data[km][j] + Ezy_yp.data[i][j];
        }
    }
    // 7
    im = pje - 1;
    for (j = pjs; j < im; j++) {
        jm = j + 1;
        km = j - pjs;
        for (i = 0; i < pis; i++) {
            ez.data[i + 1][jm] = Ezx_xn.data[i][j] + Ezy_xn.data[i][km];
        }
    }
    // 8
    for (j = pjs; j < pje - 1; j++) {
        jm = j + 1;
        km = j - pjs;
        for (i = 0; i < PMLCellxp; i++) {
            ez.data[i + pie][jm] = Ezx_xp.data[i][j] + Ezy_xp.data[i][km];
        }
    }
}

void UpdEltFldForPML_TMz(MyStruct &ex, MyStruct &ey, const MyStruct &hz) {
    unsigned int i, j;
    unsigned int im, jm;
    unsigned km = 0;
    if (IsPMLxn) {
        for (i = 1, im = 0; i <= pis; i++, im++) {
            for (j = 0; j < ey.ny; j++) {
                ey.data[i][j] = Ceye_xn.data[im][j] * ey.data[i][j]
                        + Ceyhz_xn.data[im][j]*(hz.data[i][j] - hz.data[im][j]);
            }
        }
    }
    if (IsPMLxp) {
        for (i = pie, im = 0; i < nx; i++, im++) {
            km = i - 1;
            for (j = 0; j < ny; j++) {
                ey.data[i][j] = Ceye_xp.data[im][j] * ey.data[i][j]
                        + Ceyhz_xp.data[im][j]*(hz.data[i][j] - hz.data[km][j]);

            }
        }
    }
    if (IsPMLyn) {
        for (j = 1, jm = 0; j <= pjs; j++, jm++) {
            for (i = 0; i <= nxm1; i++)
                ex.data[i][j] = Cexe_yn.data[i][jm] * ex.data[i][j]
                    + Cexhz_yn.data[i][jm]*(hz.data[i][j] - hz.data[i][jm]);

        }
    }
    if (IsPMLyp) {
        for (j = pje, jm = 0; j < ny; j++, jm++) {
            km = j - 1;
            for (i = 0; i <= nxm1; i++)
                ex.data[i][j] = Cexe_yp.data[i][jm] * ex.data[i][j]
                    + Cexhz_yp.data[i][jm]*(hz.data[i][j] - hz.data[i][km]);
        }
    }
}

void UpdateEFieldForPML(MyStruct &ex, MyStruct &ey, MyStruct &ez, const MyStruct &hx, const MyStruct &hy, const MyStruct &hz) {
    /* update electric fields at PML regions */

    if (IsTEz) {
        UpdEltFldForPML_TEz(ez, hx, hy);
    }
    if (IsTMz) {
        UpdEltFldForPML_TMz(ex, ey, hz);
    }

}

void UpdMagFldForPML_TEz(MyStruct &hx, MyStruct &hy, const MyStruct &ez) {
    unsigned int i, j;
    unsigned int im, jm;
    if (IsPMLxn) {
        for (j = 1, jm = 0; j <= nym1; j++, jm++) {
            for (i = 0; i < pis; i++) {
                hy.data[i][j] = Chyh_xn.data[i][jm] * hy.data[i][j]
                        + Chyez_xn.data[i][jm]*(ez.data[(i + 1)][j] - ez.data[i][j]);
                /*******************************************************************************************************************/
            }
        }
    }
    if (IsPMLxp) {
        for (j = 1, jm = 0; j <= nym1; j++, jm++) {
            for (i = pie, im = 0; i <= nxm1; i++, im++) {
                hy.data[i][j] = Chyh_xp.data[im][jm] * hy.data[i][j]
                        + Chyez_xp.data[im][jm]*(ez.data[(i + 1)][j] - ez.data[i][j]);
                /*******************************************************************************************************************/
            }
        }
    }

    if (IsPMLyn) {
        for (i = 1, im = 0; i <= nxm1; im++, i++) {
            for (j = 0; j < pjs; j++) {
                hx.data[i][j] = Chxh_yn.data[im][j] * hx.data[i][j]
                        + Chxez_yn.data[im][j]*(ez.data[i][j + 1] - ez.data[i][j]);
            }
        }
    }

    if (IsPMLyp) {
        for (i = 1, im = 0; i <= nxm1; i++, im++) {
            for (j = pje, jm = 0; j <= nym1; j++, jm++) {
                hx.data[i][j] = Chxh_yp.data[im][jm] * hx.data[i][j]
                        + Chxez_yp.data[im][jm]*(ez.data[i][j + 1] - ez.data[i][j]);
            }
        }
    }
}

void UpdMagFldForPML_TMz(MyStruct &hz, const MyStruct &ex, const MyStruct &ey) {
    unsigned int i, j;
    unsigned int im, jm;
    unsigned km;
    if (IsPMLxn) {
        for (i = 0, im = i + 1; i < Hzx_xn.nx; ++i, im++)
            for (j = 0; j < Hzx_xn.ny; j++)
                Hzx_xn.data[i][j] = Chzxh_xn.data[i][j] * Hzx_xn.data[i][j]
                    + Chzxey_xn.data[i][j]*(ey.data[im][j] - ey.data[i][j]);
        for (i = 0; i < Hzy_xn.nx; ++i)
            for (j = 0, jm = pjs; j < Hzy_xn.ny; j++, jm++) {
                Hzy_xn.data[i][j] = Chzyh_xn.data[i][j] * Hzy_xn.data[i][j]
                        + Chzyex_xn.data[i][j]*(ex.data[i][jm + 1] - ex.data[i][jm]);
            }
    }
    if (IsPMLxp) {
        for (i = 0, im = pie; i < Hzx_xp.nx; ++i, im++) {
            km = im + 1;
            for (j = 0; j < Hzx_xp.ny; j++) {
                Hzx_xp.data[i][j] = Chzxh_xp.data[i][j] * Hzx_xp.data[i][j]
                        + Chzxey_xp.data[i][j]*(ey.data[km][j] - ey.data[im][j]);
            }
            //if(fabs(Hzx_xp.data[i][j])>(1.8e-2)*eps)
            //	puts("Hello!");
        }
        for (i = 0, im = pie; i < Hzy_xp.nx; ++i, im++) {
            for (j = 0, jm = pjs; j < Hzy_xp.ny; j++, jm++) {
                Hzy_xp.data[i][j] = Chzyh_xp.data[i][j] * Hzy_xp.data[i][j]
                        + Chzyex_xp.data[i][j]*(ex.data[im][jm + 1] - ex.data[im][jm]);
                //if(fabs(Hzx_xp.data[i][j])>(1.8e-2)*eps)
                //	puts("Hello!");
            }
        }
    }
    /*===============================================================================================================*/
    if (IsPMLyn) {
        for (i = 0, im = pis; i < Hzx_yn.nx; ++i, im++) {
            jm = im + 1;
            for (j = 0; j < Hzx_yn.ny; j++) {
                Hzx_yn.data[i][j] = Chzxh_yn.data[i][j] * Hzx_yn.data[i][j]
                        + Chzxey_yn.data[i][j]*(ey.data[jm][j] - ey.data[im][j]);
            }
        }
        for (j = 0, jm = j + 1; j < Hzy_yn.ny; j++, jm++) {
            for (i = 0; i < Hzy_yn.nx; ++i) {
                Hzy_yn.data[i][j] = Chzyh_yn.data[i][j] * Hzy_yn.data[i][j]
                        + Chzyex_yn.data[i][j]*(ex.data[i][jm] - ex.data[i][j]);
            }
        }

    }
    /*===============================================================================================================*/
    if (IsPMLyp) {
        for (i = 0, im = pis; i < Hzx_yp.nx; im++, ++i) {
            km = im + 1;
            for (j = 0, jm = pje; j < Hzx_yp.ny; j++, jm++) {
                Hzx_yp.data[i][j] = Chzxh_yp.data[i][j] * Hzx_yp.data[i][j]
                        + Chzxey_yp.data[i][j]*(ey.data[km][jm] - ey.data[im][jm]);
                //if(fabs(Hzx_yp.data[i][j])>((1e-2)*eps))
                //	i=i;
            }
        }
        for (j = 0, jm = pje; j < Hzy_yp.ny; j++, jm++) {
            km = jm + 1;
            for (i = 0; i < Hzy_yp.nx; ++i) {
                Hzy_yp.data[i][j] = Chzyh_yp.data[i][j] * Hzy_yp.data[i][j]
                        + Chzyex_yp.data[i][j]*(ex.data[i][km] - ex.data[i][jm]);
                //if(fabs(Hzy_yp.data[i][j])>((1e-2)*eps))
                //	i=i;
            }
        }
    }
    /*===================================================================================================================*/
    for (i = 0; i < pis; i++) {
        for (j = 0; j < pjs; j++) {
            hz.data[i][j] = Hzx_xn.data[i][j] + Hzy_yn.data[i][j];
        }
    }
    for (j = pje; j < ny; j++) {
        jm = j - pje;
        for (i = 0; i < pis; i++) {
            hz.data[i][j] = Hzx_xn.data[i][j] + Hzy_yp.data[i][jm];
        }
    }
    for (i = pie; i < nx; i++) {
        im = i - pie;
        for (j = 0; j < pjs; j++) {
            hz.data[i][j] = Hzx_xp.data[im][j] + Hzy_yn.data[i][j];
        }
    }
    for (i = pie; i < nx; i++) {
        im = i - pie;
        for (j = pje; j < ny; j++) {
            hz.data[i][j] = Hzx_xp.data[im][j] + Hzy_yp.data[i][j - pje];
        }
    }
    /*=======================================================================================================*/
    for (j = pjs; j < pje; j++) {
        jm = j - pjs;
        for (i = 0; i < pis; i++) {
            hz.data[i][j] = Hzx_xn.data[i][j] + Hzy_xn.data[i][jm];
        }
    }
    for (i = 0; i < PMLCellxp; i++) {
        im = i + pie;
        for (j = pjs; j < pje; j++) {
            hz.data[im][j] = Hzx_xp.data[i][j] + Hzy_xp.data[i][j - pjs];
        }
    }
    for (i = pis; i < pie; i++) {
        im = i - pis;
        for (j = 0; j < pjs; j++) {
            hz.data[i][j] = Hzx_yn.data[im][j] + Hzy_yn.data[i][j];
        }
    }
    for (i = pis; i < pie; i++) {
        im = i - pis;
        for (j = 0; j < PMLCellyp; j++) {
            hz.data[i][j + pje] = Hzx_yp.data[im][j] + Hzy_yp.data[i][j];
        }
    }
}


