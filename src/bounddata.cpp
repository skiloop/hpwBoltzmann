
#ifndef BOUNDDATA_H
#define BOUNDDATA_H

#include "datastruct.h"
//PML BOUND

int IsPMLxn = 0;
int IsPMLxp = 0;
int IsPMLyn = 0;
int IsPMLyp = 0;

int IsAnySidePML = 0;

unsigned int PMLCellxn = 0;
unsigned int PMLCellxp = 0;
unsigned int PMLCellyn = 0;
unsigned int PMLCellyp = 0;

int pml_order = 2;
MyDataF pml_R_0 = 1e-8;

unsigned int pis, pie;
unsigned int pjs, pje;

#ifndef COORDINATES
#define COORDINATES MyStruct
#endif
//PML fields
COORDINATES Hzx_xn;
COORDINATES Hzx_xp;

COORDINATES Hzx_yn;
COORDINATES Hzx_yp;

COORDINATES Hzy_xn;
COORDINATES Hzy_xp;

COORDINATES Hzy_yn;
COORDINATES Hzy_yp;


COORDINATES Ceye_xn;
COORDINATES Ceyhz_xn;
COORDINATES Chzxh_xn;
COORDINATES Chzyh_xn;
COORDINATES Chzyex_xn;
COORDINATES Chzxey_xn;

COORDINATES Ceye_xp;
COORDINATES Ceyhz_xp;
COORDINATES Chzxh_xp;
COORDINATES Chzyh_xp;
COORDINATES Chzyex_xp;
COORDINATES Chzxey_xp;

COORDINATES Cexe_yp;
COORDINATES Cexhz_yp;
COORDINATES Chzxh_yp;
COORDINATES Chzyh_yp;
COORDINATES Chzxey_yp;
COORDINATES Chzyex_yp;

COORDINATES Cexe_yn;
COORDINATES Cexhz_yn;
COORDINATES Chzxh_yn;
COORDINATES Chzyh_yn;
COORDINATES Chzxey_yn;
COORDINATES Chzyex_yn;


COORDINATES Ezx_xn;
COORDINATES Ezx_xp;

COORDINATES Ezx_yn;
COORDINATES Ezx_yp;

COORDINATES Ezy_xn;
COORDINATES Ezy_xp;

COORDINATES Ezy_yn;
COORDINATES Ezy_yp;

COORDINATES Cezxe_xn;
COORDINATES Cezxhy_xn;
COORDINATES Chyh_xn;
COORDINATES Chyez_xn;
COORDINATES Cezye_xn;
COORDINATES Cezyhx_xn;

COORDINATES Cezye_yn;
COORDINATES Cezyhx_yn;
COORDINATES Cezxe_yn;
COORDINATES Cezxhy_yn;
COORDINATES Chxh_yn;
COORDINATES Chxez_yn;


COORDINATES Cezxe_xp;
COORDINATES Cezxhy_xp;
COORDINATES Cezye_xp;
COORDINATES Cezyhx_xp;
COORDINATES Chyh_xp;
COORDINATES Chyez_xp;

COORDINATES Cezxe_yp;
COORDINATES Cezxhy_yp;
COORDINATES Cezye_yp;
COORDINATES Cezyhx_yp;
COORDINATES Chxh_yp;
COORDINATES Chxez_yp;

#endif
