/*
 *		This file defines functions that update scatter fields.
 *
 *
 *
 *
 */

#ifndef UPDATEFIELDS_H
#define UPDATEFIELDS_H

#include "datastruct.h"
#include "commondata.h"
///////////////////////////////
//formation 3
///////////////////////////////

void SumErms();
void UpdateErms();
void UpdateEField();
void UpdateMField();
void UpdateUField();
void UpdateUeField();
void UpdateUCoeffiecients();
void CapFields(unsigned int step);
void DensityBound(MyStruct stru, int bndwidth, const int);
void UpdateDensity();

inline MyDataF ElectricEnergy(unsigned im, unsigned jm) {
    if (Ne.data[im][jm] < 1e-5)
        return 0;
    else
        return Ue.data[im][jm] / Ne.data[im][jm];
}

void CalEmax();
void InterpEmax();
MyDataF getEmax(int i, int j) ;

#endif
