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
 * @file datastruct.cpp
 * @version 0.0.0
 * @author skiloop ( skiloop@126.com )
 * @date 31/08/2011 0.0.0 created, by skiloop
 */


#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "datastruct.h"

using namespace std;

int MyStruct::cnt = 0;
string MyStruct::tail = ".dat";

#ifdef MATLAB_SIMULATION
Engine* MyStruct::ep = NULL;
#endif

MyStruct::MyStruct(unsigned cx, unsigned cy)
: nx(cx), ny(cy) {
    unsigned i;
    if (cx == 0) {
        data = NULL;
        return;
    }
    data = new pMyDataF[cx];
    if (data == NULL) {
        printf("Failed to create space for MyStruct!\n");
        exit(0);
    }
    for (i = 0; i < cx; i++) {
        if (cy == 0)data[i] = NULL;
        else {
            data[i] = new MyDataF[cy];
            if (data[i] == NULL) {
                printf("Failed to create space for MyStruct!\n");
                exit(0);
            }
        }
    }
}

MyStruct::MyStruct(const MyStruct& obj) : data(NULL) {
    CreateStruct(obj);
    try {
        for (unsigned i = 0; i < obj.nx; i++) {
            memcpy(data[i], obj.data[i], obj.ny * sizeof (MyDataF));
        }
    } catch (exception & e) {
        cerr << e.what() << endl;
        exit(-1);
    }
}

MyStruct::~MyStruct() {

    if (data == NULL)
        return;
    unsigned i;
    for (i = 0; i < nx; i++) {
        if (data[i] != NULL) {
            delete [] data[i];
        }
    }
    //delete []p;
    delete []data;
}

int MyStruct::CreateStruct(unsigned nnx, unsigned nny) {
    unsigned i;
    data = new pMyDataF[nnx];
    if (data == NULL) {
        printf("Failed to create space for MyStruct!\n");
        return -1;
    }
    for (i = 0; i < nnx; i++) {
        data[i] = new MyDataF[nny]();
        if (data[i] == NULL) {
            printf("Failed to create space for MyStruct!\n");
            return -2;
        }
    }
    nx = nnx;
    ny = nny;

    return 0;
}

int MyStruct::CreateStruct(unsigned nnx, unsigned nny, MyDataF initVal) {
    unsigned i;
    data = new pMyDataF[nnx];
    if (data == NULL) {
        printf("Failed to create space for MyStruct!\n");
        return -1;
    }
    for (i = 0; i < nnx; i++) {
        data[i] = new MyDataF[nny]();
        if (data[i] == NULL) {
            printf("Failed to create space for MyStruct!\n");
            return -2;
        }
    }
    nx = nnx;
    ny = nny;
    return ResetStructData(initVal);
}

int MyStruct::ResetStructData(MyDataF Val) {
    unsigned i, j;
    if (!CheckStruct())
        return -1;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            data[i][j] = Val;
        }
    }
    return 0;
}

void MyStruct::PrintData() {
    unsigned i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++)
            printf("%2.3e ", data[i][j]);
        printf("\n");
    }
    //system("pause");
}

int MyStruct::BackupMyStruct(const MyStruct &obj) {

    if (data == NULL || obj.data == NULL) {
        printf("One of the fields is empty!\n");
        return 0;
    }
    if (data == obj.data || this == &obj) {
        printf("Cannot Copy to Itself!\n");
        return 0;
    }

    if (ny != obj.ny || nx != obj.nx) {
        printf("Fields are not the same size!\n");
        return 0;
    }
    if (nx == 0 || ny == 0) {
        printf("Empty space to be copied from!\n");
        return 0;
    }
    for (unsigned i = 0; i < nx; i++)
        memcpy(obj.data[i], data[i], ny * sizeof (MyDataF));
    return 1;

}

bool MyStruct::CheckStruct() {
    if (ny <= 0 || nx <= 0 || data == NULL) {
        return false;
    }
    return true;
}

void MyStruct::CaptData(const unsigned num, unsigned p) {
    unsigned i, j;

    stringstream ss;
    ss << name << num << tail;
    ofstream ofile(ss.str().c_str());

    if (!ofile.is_open()) {
        cerr << "Cannot open" << ss << endl;
        return;
    }
    //check p
    if (p > nx || p > ny)p = 0;

    p += 1;
    for (i = 0; i < nx; i += p) {
        for (j = 0; j < ny; j += p) {
            ofile << data[i][j] << '\t';
        }
        ofile << endl;
    }
    ofile.close();
}

void MyStruct::FreeStructData() {
    unsigned i;
    if (data == NULL)
        return;
    for (i = 0; i < nx; i++)
        delete [] data[i];
    delete [] data;
}

void MyStruct::operator =(MyStruct const &other) {
    unsigned i, j;
    if (&other == this || data == NULL || other.data == NULL || other.nx != nx || other.ny != ny)
        return;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            data[i][j] = other.data[i][j];
        }
    }
}

void MyStruct::InitStructData(MyDataF initVal) {
    unsigned i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            data[i][j] = initVal;
        }
    }
}

void MyStruct::SaveData(unsigned leap) {
    string fname = name + tail;
    leap += 1;
    if (leap >= nx || leap >= ny) {
        cerr << "Invalid leap for saving data!" << endl;
        return;
    }
    ofstream out(fname.c_str(), ios_base::binary);
    if (!out.is_open()) {
        cerr << "File " << fname << "cannot be opened!" << endl;
        return;
    }
    unsigned i, j;
    for (i = 0; i < nx; i += leap) {
        for (j = 0; j < ny; j += leap)
            out << data[i][j] << '\t';
        out << endl;
    }
}

int MyStruct::InitMatlabEngine() {
#ifdef MATLAB_SIMULATION
    if (ep != NULL)return -2;
    if ((ep = engOpen(NULL)) == NULL) {
        cerr << "Can't start matlab engine!" << endl;
        exit(-1);
    }
    engSetVisible(ep, 1);
#endif
    return 0;

}

int MyStruct::CloseEngine() {
#ifdef MATLAB_SIMULATION
    engEvalString(ep, "close all;clear;");
    engClose(ep);
#endif
    return 0;
}

int MyStruct::CreateStruct(const MyStruct &stru) {
    return CreateStruct(stru.nx, stru.ny);
}

int MyStruct::CreateStruct(const MyStruct &stru, MyDataF initVal) {
    return CreateStruct(stru.nx, stru.ny, initVal);
}

void MyStruct::ClearSim() {
#ifdef MATLAB_SIMULATION
    mxDestroyArray(MyArray);
    mxDestroyArray(num);
#endif
}

void MyStruct::PlotArrays() {

#ifdef MATLAB_SIMULATION
    MyDataF *pData = (MyDataF*) malloc(nx * ny * sizeof (MyDataF));
    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            pData[i * ny + j] = data[i][j];
        }
    }
    engPutVariable(ep, "ind", num);
    engEvalString(ep, "ind=int32(ind);");
    memcpy(mxGetPr(MyArray), pData, nx * ny * sizeof (MyDataF));
    engPutVariable(ep, "array", MyArray);
    engEvalString(ep, "obj(ind).array=array;clear array;");
    engEvalString(ep, "set(obj(ind).img,'CData',obj(ind).array);drawnow;");
    free(pData);
#endif
}

void MyStruct::InitPlot() {
    MyStruct::cnt++;
    string filename = name + tail;
    Number = MyStruct::cnt;
#ifdef MATLAB_SIMULATION
    mxArray *mxStr = mxCreateString(filename.c_str());
    MyDataF *pData = (MyDataF*) malloc(nx * ny * sizeof (MyDataF));
    if (pData == NULL)return;
    MyArray = mxCreateDoubleMatrix(ny, nx, mxREAL);
    num = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    memcpy(mxGetPr(num), &Number, sizeof (unsigned));
    engPutVariable(ep, "ind", num);
    engEvalString(ep, "ind=int32(ind);");
    engEvalString(ep, "obj(ind).fig=figure('NumberTitle','OFF');");

    engPutVariable(ep, "name", mxStr);
    engEvalString(ep, "obj(ind).name=name;");

    for (unsigned i = 0; i < nx; i++) {
        for (unsigned j = 0; j < ny; j++) {
            pData[i * ny + j] = data[i][j];
        }
    }
    memcpy(mxGetPr(MyArray), pData, nx * ny * sizeof (MyDataF));
    engPutVariable(ep, "array", MyArray);

    engEvalString(ep, "obj(ind).array=array;clear array;");
    engEvalString(ep, "obj(ind).img=imagesc(obj(ind).array);obj(ind).ax=gca;title(obj(ind).ax,obj(ind).name);drawnow;");
    engEvalString(ep, "set(gca,'YDir','Normal');colorbar;");
    free(pData);
    mxDestroyArray(mxStr);
#endif
}

