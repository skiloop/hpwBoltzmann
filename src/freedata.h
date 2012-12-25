#ifndef FREEDATA_H
#define FREEDATA_H
#include "commondata.h"
void closeFiles(){
	denfile.close();
	deff_file.close();
}
void freedata()
{
    if(IsTMz)
    {

    }
    if(IsTEz)
    {
    }
	
}


void savedata()
{
	closeFiles();
    if(IsTMz)
    {
        Ex.SaveData();
        Ey.SaveData();
        Hz.SaveData();
        Ux.SaveData();
        Uy.SaveData();
    }
    if(IsTEz)
    {
        Ex.SaveData();
        Ey.SaveData();
        Hz.SaveData();
        Ux.SaveData();
        
    }
    if(ifWithDensity){
	    Ne.SaveData();
	    Ue.SaveData();	
    }
}

#endif
