/* 
 * File:   FileInterp.cpp
 * Author: skiloop
 * 
 * Created on April 19, 2012, 3:45 PM
 */
#include <fstream>
#include <iostream>
#include <exception>
#include <cstring>
#include <cstdlib>

#include "FileInterp.h"

using namespace std;

FileInterp::FileInterp(int nv, const string &filex, const string &filey)
: n(nv)
, x(NULL)
, y(NULL) {
    //files
    ifstream infx, infy;

    //open file for array x
    try
    {
        infx.open(filex.c_str());
        if (!infx.is_open()) {
            cerr << "Cannot open file:" << filex << endl;
            exit(-1);
        }
#ifdef _DEBUG
        cout << "open file :" << filex << endl;
#endif
    }

    catch(exception & e) {
        cerr << e.what() << endl;
        exit(-1);
    }
    //open file for array y
    try
    {
        infy.open(filey.c_str());
        if (!infy.is_open()) {
            cerr << "Cannot open file:" << filey << endl;
            exit(-1);
        }
#ifdef _DEBUG
        cout << "open file" << filey << endl;
#endif
    }
    catch(exception & e) {
        cerr << e.what() << endl;
        exit(-1);
    }

    //Create space for Value Arrays:x and y
    if (n <= 0) {
        cerr << "Invalid n for FileInterp" << endl;
        exit(-1);
    }
    try
    {
        x = new double[n];
        y = new double[n];
    }

    catch(exception & e) {
        cerr << e.what() << endl;
        if (x != NULL)delete[]x;
        if (y != NULL)delete[]y;
        exit(-1);
    }

    //read files
    try
    {
        int i = 0;
        while (i < n && infx.good() && infy.good()) {
            infx >> x[i];
            infy >> y[i];
            i++;
        }
        n = i;
    }

    catch(exception & e) {
        cerr << e.what() << endl;
        exit(-1);
    }
    //close files
    infx.close();
    infy.close();
}

FileInterp::FileInterp(const FileInterp& orig) {

    n = orig.n;
    //create array for x and  y
    if (n <= 0) {
        cerr << "Invalid n for FileInterp" << endl;
        exit(-1);
    }
    try
    {
        x = new double[n];
        y = new double[n];
    }

    catch(exception & e) {
        cerr << e.what() << endl;
        if (x != NULL)delete[]x;
        if (y != NULL)delete[]y;
        exit(-1);
    }

    //copy values to x and y
    memcpy(x, orig.x, n * sizeof (double));
    memcpy(y, orig.y, n * sizeof (double));
}

FileInterp::~FileInterp() {
    delete []x;
    delete []y;
}

