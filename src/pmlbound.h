
#ifndef PMLBOUND_H
#define PMLBOUND_H
#include "datastruct.h"

/**
 * InitialPML
 * @brief Define PML,then create and initial PML arrays
 * @param ncpml
 * @param ncx
 * @param ncy
 */
void InitialPML(unsigned ncpml, unsigned ncx, unsigned ncy);

/**
 * UpdMagFldForPML_TEz
 * update magnetic fields at PML regions
 * TEz
 * @param hx component Hx
 * @param hy component Hy
 * @param ez component Ez
 */
void UpdMagFldForPML_TEz(MyStruct &hx, MyStruct &hy, const MyStruct &ez);

/**
 *UpdMagFldForPML_TMz
 *update magnetic fields at PML regions
 *TMz
 */
void UpdMagFldForPML_TMz(MyStruct &hz, const MyStruct &ex, const MyStruct &ey);
/**
 *UpdMEltFldForPML_TMz
 *update electric fields at PML regions
 *TEz
 */
void UpdEltFldForPML_TMz(MyStruct &ex, MyStruct &ey, const MyStruct &hz);
/**
 *UpdMEltFldForPML_TEz
 *update electric fields at PML regions
 *TMz
 */
void UpdEltFldForPML_TEz(MyStruct &ez, const MyStruct &hx, const MyStruct &hy);
void UpdateMFieldForPML(MyStruct &hx, MyStruct &hy, MyStruct &hz, const MyStruct &ex, const MyStruct &ey, const MyStruct &ez);
void UpdateEFieldForPML(MyStruct &ex, MyStruct &ey, MyStruct &ez, const MyStruct &hx, const MyStruct &hy, const MyStruct &hz);
void FreePMLSpace();

#endif
