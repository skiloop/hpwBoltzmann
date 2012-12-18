#include <iostream>

#include "datastruct.h"
#include "initials.h"
#include "pmlbound.h"
#include "commondata.h"
#include "connectingdef.h"
#include "fdtdloop.h"
#include "freedata.h"

using namespace std;

int main(int argc, char **argv) {
    //system("dir");
    if(argc<=1)
        Initial();
    else
        Initial(argv[1]);
    //system("pause");
    InitialPML(nbound,nx,ny);
    InitConnectingInterface(phi);
    //CalculateDelay(phi);
    fdtdloop();
    FreeDelayArrays();
    FreePMLSpace();
    savedata();
    freedata();
    //system("pause");
    return 0;
}
