#include "def.h"

void ReadGlobalParams( const char filename[30] = "output/GlobalBWparams.txt" )
{
    ifstream txtFile;
    txtFile.open(filename);

    int charge;
    while (true)
    {
        for (int centr: CENTR)
        {
        txtFile >> charge >> centr 
                >> paramsGlobal[charge][centr][0] >> paramsGlobal[charge][centr][1]
                >> paramsGlobal[charge][centr][2] >> paramsGlobal[charge][centr][3] >> paramsGlobal[charge][centr][4];
        }
        
        if( txtFile.eof() ) break;
    }
    
    
    txtFile.close();
}

void ReadParams( const char filename[30] = "output/BWparams.txt" )
{
    ifstream f;
    f.open(filename);

    int p, c;

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            f >> p >> c >> constPar[p][c] >> Tpar[p][c] >> Tpar_err[p][c] >> utPar[p][c] >> utPar_err[p][c];
        }
    }
    f.close();
}
