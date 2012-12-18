
#include <math.h>

#include "interp.h"


//-------------------------------------------------------------------------------

double Interpolation::lgr(double const *lgrx, double const *lgry, int lgrn, double lgrt) {
    int lgri, lgrj, lgrk, lgrm;
    double lgrz, lgrs;
    lgrz = 0.0;
    if (lgrn < 1) return (lgrz);
    if (lgrn == 1) {
        lgrz = lgry[0];
        return (lgrz);
    }
    if (lgrn == 2) {
        lgrz = (lgry[0]*(lgrt - lgrx[1]) - lgry[1]*(lgrt - lgrx[0])) / (lgrx[0] - lgrx[1]);
        return (lgrz);
    }
    lgri = 0; //上限
    while ((lgrx[lgri] < lgrt) && (lgri < lgrn)) lgri = lgri + 1;
    lgrk = lgri - 4; //下限
    if (lgrk < 0) lgrk = 0;
    lgrm = lgri + 3;
    if (lgrm > (lgrn - 1)) lgrm = lgrn - 1;
    for (lgri = lgrk; lgri <= lgrm; lgri++) {
        lgrs = 1.0;
        for (lgrj = lgrk; lgrj <= lgrm; lgrj++)
            if (lgrj != lgri) lgrs = lgrs * (lgrt - lgrx[lgrj]) / (lgrx[lgri] - lgrx[lgrj]);
        lgrz = lgrz + lgrs * lgry[lgri];
    }
    return (lgrz);
}
//---------------------------------------------------------------------------------

double Interpolation::lg3(double const *x, double const *y, int n, double t) {
    int i, j, k, m;
    double z, s;
    z = 0.0;
    if (n < 1) return (z);
    if (n == 1) {
        z = y[0];
        return (z);
    }
    if (n == 2) {
        z = (y[0]*(t - x[1]) - y[1]*(t - x[0])) / (x[0] - x[1]);
        return (z);
    }
    if (t <= x[1]) {
        k = 0;
        m = 2;
    } else if (t >= x[n - 2]) {
        k = n - 3;
        m = n - 1;
    } else {
        k = 1;
        m = n;
        while (m - k != 1) {
            i = (k + m) / 2;
            if (t < x[i - 1]) m = i;
            else k = i;
        }
        k = k - 1;
        m = m - 1;
        if (fabs(t - x[k]) < fabs(t - x[m])) k = k - 1;
        else m = m + 1;
    }
    z = 0.0;
    for (i = k; i <= m; i++) {
        s = 1.0;
        for (j = k; j <= m; j++)
            if (j != i) s = s * (t - x[j]) / (x[i] - x[j]);
        z = z + s * y[i];
    }
    return (z);
}
//-----------------------Linear_interpolation----------------------------------------------------
/*

double Interpolation::Linear_interpolation(double *x, double *y, int n, double t)
{
int up,down,i;
double y0,y1,yt;
up=0;
down=0;



if(t<x[0]) yt=y[0];
else if(t>=x[n-1]) yt=y[n-1];
else {
   //判断上下限
  for(i=0;i<=(n-2);i++)
  {
   if(t>=x[i]) { down=i;up=i+1;}
   else break;
  }


//  cout<<"普通法"<<"i="<<i<<endl;

//   cout<<"t="<<t<<"  "<<"up="<<x[up]<<"  "<<"down="<<x[down]<<endl;

  y0=y[down];
  y1=y[up];
  yt=(t-x[up])/(x[down]-x[up])*y0+(t-x[down])/(x[up]-x[down])*y1;
 }
//   cout<<"t="<<t<<"  "<<"up"<<y[up]<<"  "<<yt<<"   "<<"down="<<y[down]<<endl;
return yt;
}
 */
//---------------------------------------------------------------------------------------

double Interpolation::Linear_interpolation(double const *x, double const *y, int n, double t) {
    int up, down, i;
    int d1, d2, mid; //查找的大小区间，中间
    double y0, y1, yt;
    up = 0;
    down = 0;


    //   ofstream time ("E:\\普通time2.dat");

    if (t <= x[0]) yt = y[0];
    else if (t >= x[n - 1]) yt = y[n - 1];
    else {
        d1 = 0;
        d2 = n - 1;
        //判断上下限
        for (i = 0; i <= (n - 2); i++) {
            mid = (d1 + d2) / 2;
            if (t >= x[mid]) {
                if (t < x[mid + 1]) {
                    down = mid;
                    up = mid + 1;
                    break;
                } else {
                    d1 = mid;
                }
            } else {
                if (t >= x[mid - 1]) {
                    down = mid - 1;
                    up = mid;
                    break;
                } else {
                    d2 = mid;
                }
            }
            // else break;
        }



        //   cout<<"对分法"<<"i="<<i<<endl;
        //     cout<<"t="<<t<<"  "<<"up"<<x[up]<<"  "<<"down="<<x[down]<<endl;
        y0 = y[down];
        y1 = y[up];
        yt = (t - x[up]) / (x[down] - x[up]) * y0 + (t - x[down]) / (x[up] - x[down]) * y1;
    }
    //         cout<<"t="<<t<<"  "<<"up"<<y[up]<<"  "<<yt<<"   "<<"down="<<y[down]<<endl;
    return yt;
}
