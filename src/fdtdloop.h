#ifndef FDTDLOOP_H
#define FDTDLOOP_H



//#define SD(l,i,j) l.data[(i)*l.ny+j]
//#define DD(s,p,i,j) (s.data[(i)*s.ny+j]+p.data[(i)*p.ny+j])

////////////////////////////////////////////////
// FUNCTION PREDEFINE
///////////////////////////////////////////////

void UpdateCoeff();
void AdBound();
void InitSources();
void fdtdloop();


#undef SD

#endif
