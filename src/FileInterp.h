/* 
 * File:   FileInterp.h
 * Author: skiloop
 *
 * Created on April 19, 2012, 3:45 PM
 */

#ifndef FILEINTERP_H
#define	FILEINTERP_H

#include <string>

#include "interp.h"

using std::string;

class FileInterp : public Interpolation {
public:
    /*
     * @param len length of the arrays
     * @param filex file name for @c x array
     * @param filey file name for @c y array
     */
    FileInterp(int len, const string &filex, const string &filey);

    FileInterp(const FileInterp& orig);
    virtual ~FileInterp();

    double Interp(double vx) {
        return Linear_interpolation(x, y, n, vx);
    }
    int n;
private:
    double *x;
    double *y;
};

#endif	/* FILEINTERP_H */

