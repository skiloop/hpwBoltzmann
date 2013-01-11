
#ifndef MICRODEF_H
#define MICRODEF_H
#define _USE_MATH_DEFINES
#include <math.h>

/////////////////////////////////
//FDTD DEFINITION
/////////////////////////////////
#define MAX_LIMIT 1e80
#define NE_MAX_LIMIT 1E307
#define COURANT_FACTOR 0.5
#define CFL_FACTOR 0.4
#define NUMCELLPERWAVELEN 50
#define MAXWELL_MESH_SIZE (1.0/NUMCELLPERWAVELEN)
#define NUMBER_OF_CELLS_IN_PML_BOUND 20
#define SCATTER_FIELD_DOMAIN_BND_SIZE 5
#define NUMBER_OF_WAVELENGTHS_IN_DOMAIN 3
#define DEN_TIME_STEPS 2000000

#define TOTAL_TIME 1e-9 // in second
#define SAVE_LEAP 30	//sample Ez and ne every SAVE_LEAP steps
#define SAVE_ERMS_LEAP 55
#define PULSE_LENGGTH_IN_TIME
#define FINE_GRID_SIZE 4
#define SCATTERED_FORMATION


/////////////////////////
// plasma default
//////////////////////////
#define NE0 1e13

#define DEFAULT_REI -1.0
#define DEFAULT_MIU_DIV 100.0
#define DEFAULT_AIR_PRESSURE 760.0
#define IF_WITH_DENSITY 0


////////////////////////
// wave type default
///////////////////////
//#define _SOURCE_TMZ_ 1
#ifndef _SOURCE_TEZ_
// if not define _SOURCE_TMZ_ then define _SOURCE_TEZ_=1
#ifndef _SOURCE_TMZ_
#define _SOURCE_TEZ_ 1
#define _SOURCE_TMZ_ 0
#elif (_SOURCE_TMZ_!=0)
#define _SOURCE_TEZ_ 0
#else
#define _SOURCE_TEZ_ 1
#endif
#endif // 

#define E_0 6e6

#define FREQUENCY 110E9
#define MAX_NE 1E29
#define INC_ANGLE 0.0*M_PI



// option
#define PROJECT_2012_10

////////////////////////////////////
//SAMPLE DEFINITION
///////////////////////////////////
#define LEAPSTEP_OF_DISPLAY 1
#define DTF_VI_LIMIT 2.0e10
#define LEAPSTEP_OF_CAPTURE 40

///////////////////////////
//SOURCES POSITION
///////////////////////////
#define CELLS_INSIDE_I 0
#define CELLS_INSIDE_S 0
#define SOURCE_POS_IN_Y 0
#define SOURCE_POS_IN_X 0
#define SOURCES_SIZE 1

/// Matlab Simulation
//#pragma comment( lib, "libmx.lib"  )
////#pragma comment( lib, "libmex.lib" )
//#pragma comment( lib, "libeng.lib" )

#endif
