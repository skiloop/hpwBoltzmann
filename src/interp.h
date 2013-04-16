
#pragma once


class Interpolation //七次拉格朗日插值
{
public:
    static double lgr(double const *lgrx, double const *lgry, int lgrn, double lgrt);
    static double lg3(double const *x, double const *y, int n, double t);
    static double Linear_interpolation(double const *x, double const*y, int n, double t);
    //   double Linear_interpolation1(double *x, double *y, int n, double t);
};


