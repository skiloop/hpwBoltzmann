/**
 * This file declare parameters of the connecting interface.
 * 
 * Create at 
 * by skiloop@126.com
 */
#ifndef __CONNECTING_H__
#define __CONNECTING_H__

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>

#include "datastruct.h"

#define INC_SIZE 0
/**
 * @brief Initialize the connecting interface:
 *  main tasks: 1 create space for connecting interfaces
 * 		2 calculate the
 */
void InitConnectingInterface(MyDataF IncAngle);

/**
 * @brief update E field along lines
 *  
 */
void UpdateConnectingE(const MyDataF t);
/**
 * @brief update M field along lines to apply connecting interface 
 */
void UpdateConnectingM(const MyDataF t);
/**
 * @brief apply E field connecting interface 
 */
void ApplyConnectingE(const MyDataF t);
/**
 * @brief apply M field connecting interface
 */
void ApplyConnectingM(const MyDataF t);

/**
 * @brief get the source wave value at time t
 */
MyDataF Source(MyDataF t);


/**
 * Calculate the delay of (@c px, @c py) from (xs,ys) the first point in total field 
 */
MyDataF delays(MyDataF px, MyDataF py);

MyDataF PhaseVelRatio(MyDataF angle);

void FreeDelayArrays();

#endif
