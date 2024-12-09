#include "def.h"

 void ReadGlobalParams( const char filename[30] = "output/GlobalBWparams.txt" )
{
    ifstream txtFile;
    txtFile.open(filename);

    for (int centr: CENTR)
    {
       txtFile >> centr 
               >> paramsGlobal[centr][0] >> paramsGlobal[centr][1]
               >> paramsGlobal[centr][2] >> paramsGlobal[centr][3] >> paramsGlobal[centr][4];
    }
    
    txtFile.close();
}
