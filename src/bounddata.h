
#ifndef BOUNDDATA_H
#define BOUNDDATA_H

#include "datastruct.h"
//PML BOUND

extern int IsPMLxn;
extern int IsPMLxp;
extern int IsPMLyn;
extern int IsPMLyp;

extern int IsAnySidePML;

extern unsigned int PMLCellxn;
extern unsigned int PMLCellxp;
extern unsigned int PMLCellyn;
extern unsigned int PMLCellyp;

extern int pml_order;
extern MyDataF pml_R_0;

extern unsigned int pis, pie;
extern unsigned int pjs, pje;

#ifndef COORDINATES
#define COORDINATES MyStruct
#endif
//PML fields
extern COORDINATES Hzx_xn;
extern COORDINATES Hzx_xp;

extern COORDINATES Hzx_yn;
extern COORDINATES Hzx_yp;

extern COORDINATES Hzy_xn;
extern COORDINATES Hzy_xp;

extern COORDINATES Hzy_yn;
extern COORDINATES Hzy_yp;


extern COORDINATES Ceye_xn;
extern COORDINATES Ceyhz_xn;
extern COORDINATES Chzxh_xn;
extern COORDINATES Chzyh_xn;
extern COORDINATES Chzyex_xn;
extern COORDINATES Chzxey_xn;

extern COORDINATES Ceye_xp;
extern COORDINATES Ceyhz_xp;
extern COORDINATES Chzxh_xp;
extern COORDINATES Chzyh_xp;
extern COORDINATES Chzyex_xp;
extern COORDINATES Chzxey_xp;

extern COORDINATES Cexe_yp;
extern COORDINATES Cexhz_yp;
extern COORDINATES Chzxh_yp;
extern COORDINATES Chzyh_yp;
extern COORDINATES Chzxey_yp;
extern COORDINATES Chzyex_yp;

extern COORDINATES Cexe_yn;
extern COORDINATES Cexhz_yn;
extern COORDINATES Chzxh_yn;
extern COORDINATES Chzyh_yn;
extern COORDINATES Chzxey_yn;
extern COORDINATES Chzyex_yn;




extern COORDINATES Ezx_xn;
extern COORDINATES Ezx_xp;

extern COORDINATES Ezx_yn;
extern COORDINATES Ezx_yp;

extern COORDINATES Ezy_xn;
extern COORDINATES Ezy_xp;

extern COORDINATES Ezy_yn;
extern COORDINATES Ezy_yp;

extern COORDINATES Cezxe_xn;
extern COORDINATES Cezxhy_xn;
extern COORDINATES Chyh_xn;
extern COORDINATES Chyez_xn;
extern COORDINATES Cezye_xn;
extern COORDINATES Cezyhx_xn;

extern COORDINATES Cezye_yn;
extern COORDINATES Cezyhx_yn;
extern COORDINATES Cezxe_yn;
extern COORDINATES Cezxhy_yn;
extern COORDINATES Chxh_yn;
extern COORDINATES Chxez_yn;


extern COORDINATES Cezxe_xp;
extern COORDINATES Cezxhy_xp;
extern COORDINATES Cezye_xp;
extern COORDINATES Cezyhx_xp;
extern COORDINATES Chyh_xp;
extern COORDINATES Chyez_xp;

extern COORDINATES Cezxe_yp;
extern COORDINATES Cezxhy_yp;
extern COORDINATES Cezye_yp;
extern COORDINATES Cezyhx_yp;
extern COORDINATES Chxh_yp;
extern COORDINATES Chxez_yp;

#endif
