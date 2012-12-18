/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  <copyright holder> <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file datastruct.h
 * @version 0.0.0
 * @author skiloop ( skiloop@126.com )
 * @date 31/08/2011 0.0.0 created, by skiloop
 */
#pragma once
#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#ifdef MATLAB_SIMULATION

#include <engine.h>
#include <mex.h>
#ifdef printf
#undef printf
#endif

#endif// end MATLAB_SIMULATION

#include <string>
#include "microdef.h"


typedef double MyDataF;
typedef double* pMyDataF;

/**
 * nx:point count in x direction
 * ny:point count in x direction
 * data:pointer to store data;
 */
class MyStruct {
public:
    static int cnt;
    static std::string tail;
#ifdef MATLAB_SIMULATION
    static Engine *ep;
#endif

public:
    unsigned int nx;
    unsigned int ny;
    MyDataF** data;
    std::string name;

private:
    int Number;
#ifdef MATLAB_SIMULATION
    mxArray *num;
    mxArray *MyArray;
#endif
public:
    /**
     * if cx==0 then data = NULL
     * else data has cx pointers;
     * 
     * if cy == 0 then all cx pointers of data is NULL
     * else data[i] has cy pointers with i from 0 to cx-1;
     * 
     * when cannot create space for data and data[i],exit program;
     */
    MyStruct(unsigned int cx = 0, unsigned int cy = 0);
    MyStruct(const MyStruct &obj);
    ~MyStruct();
    /**
     * print data in struct MyStruct
     */
    void PrintData();
    /**
     * free space created for MyStruct @c mst
     */
    void FreeStructData();

    /**
     * Set Data to val
     */
    int ResetStructData(MyDataF val = 0);

    /**
     * check data of MyStruct @c mst is valid
     * if data is not NULL and none of its subpointers,then 
     * return true,otherwise false
     */
    bool CheckStruct();

    /**
     * Create Space for struct MyStruct and initialize its @c nx and @c ny
     */
    int CreateStruct(unsigned nnx, unsigned nny);
    /**
     * Create Space for struct MyStruct and initialize its @c nx and @c ny
     */
    int CreateStruct(unsigned nnx, unsigned nny, MyDataF initVal);

    /**
     * Copy all data in st to stpre
     * Dimensions of @c st and that of @c pstruct must macth,and both with valid 
     * data
     */
    int BackupMyStruct(const MyStruct &mstru);
   
    /**
     * @brief Save data of  MyStruct data skipping p rows and p columns 
     * 
     */
    void CaptData(const unsigned num, unsigned p=0);

    void operator=(MyStruct const &other);
    void InitStructData(MyDataF initVal = 0);
    void SaveData(unsigned leap=0);

    /**
     * @brief Create a MyStruct with the same size;
     * @param stru the source MyStruct to be copied.
     * @return
     */
    int CreateStruct(const MyStruct &stru);
    /**
     * @brief Create a MyStruct with the same size as @c stru and initial all var to @c initVal;
     * @param stru the source MyStruct to be copied.
     */
    int CreateStruct(const MyStruct &stru, MyDataF initVal);
    
    //set name
    void SetName(const std::string &sn){name = sn;}
public:
    void ClearSim();
    void PlotArrays();
    void InitPlot();
public:
    static int InitMatlabEngine();
    static int CloseEngine();

};
#endif
