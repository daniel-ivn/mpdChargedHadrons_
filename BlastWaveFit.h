#include "def.h"
#include "WriteReadFiles.h"

using namespace std;

void setParamsForSys ( int systematicType, double parResults[4] )
{
    switch(systematicType)
    {
        case 1: 
            parResults[1] *= 0.9;
        case 2: 
            parResults[1] *= 1.1;
        case 3:
            parResults[2] *= 0.9;
        case 4:
            parResults[2] *= 1.1;
        case 5:
            parResults[0] *= 0.9;
        case 6:
            parResults[0] *= 1.1;
        // case 7:
        //     xmin[part] *= 0.9;
        //     xmax[part] *= 0.9;
        // case 8:
        //     xmin[part] *= 1.1;
        //     xmax[part] *= 1.1;
    }
}

void GetContourPlots( int part, int centr )
{
    for (int s = 1; s < N_SIGMA; s++)
    {
        gMinuit->SetErrorDef(s * 8); //note 4 and not 2!
        contour[part][centr][s] = (TGraph*)gMinuit->Contour(40, 2, 1);
        if (contour[part][centr][s]) 
        {
            contour[part][centr][s]->SetLineColor(partColors[part] + centr);
            contour[part][centr][s]->SetFillStyle(0);
        }
        else 
            cout << "ERROR " << part << " " <<centr << " " << s << endl;
        
    }
}

class BlastWaveFit {
public:

    bool isContour = false;
    bool isDraw = true;
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];

    void Fit( int initParamsType = 0, int systematicType = 0  )
    {    
        // ++++++ Read data +++++++++++++++++++++++++++++++++++++

        string inputFileName = "postprocess_mpdpid10";
        SetSpectra(inputFileName, "mt");

        // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

        TVirtualFitter::SetDefaultFitter("Minuit");  
        // TMinuit* minuit = new TMinuit(5); 
        auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
        funcx->SetParameters(2,1);
        funcx->SetParNames("constant", "T", "beta", "mass", "pt");
        MyIntegFunc *integ = new MyIntegFunc(funcx);

        for (int part: PARTS)
        {  
            for (int centr: CENTR)
            {
                string ifuncxName = "MyIntegFunc_" + to_string(part) + "_" + to_string(centr);
                ifuncx[part][centr] = new TF1("ifuncx", integ, xmin[part], xmax[part], 4, ifuncxName.c_str());

                switch(initParamsType)
                {
                    case 0: { /* DEFAULT */
                        // ================== version1 Params from Global fit ============================
                        double parResults[5];
                        ReadGlobalParams(paramsGlobal);
                        getGlobalParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        parResults[2] = parResults[2] > 0.8 ? 0.8 : parResults[2];
                        parResults[2] = (parResults[2] > 0.6) ? 0.6 : parResults[2];
                        ifuncx[part][centr]->SetParameters(parResults);
                        for (int par = 0; par < 3; par++)
                        {
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.6, parResults[par] * 1.3);
                            if (part <= 1 && par == 2) 
                                ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par]);
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================

                    case 1: {
                        // ================== version2 Params from individual fit results =====================
                        double parResults[4];
                        ReadParams(part, centr, parResults);
                        if (parResults[0] == 0)
                            continue;
                            
                        ifuncx[part][centr]->SetParameters(parResults);
                        for (int par = 0; par < 3; par++)
                        {
                            ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.6, parResults[par] * 1.5);
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================

                    case 2: {
                        // ================= version 2 Params with limits =================================
                        double customParams[4] = {con[part], 0.09, 0.75, masses[part]};
                        ifuncx[part][centr]->SetParameters(customParams);
                        ifuncx[part][centr]->SetParLimits(0, conmin[part], conmax[part]);
                        ifuncx[part][centr]->SetParLimits(1, 0.8, 0.1);	
                        ifuncx[part][centr]->SetParLimits(2, 0.5, 0.8);	
                        ifuncx[part][centr]->FixParameter(3, masses[part]);	//	mass
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }   
                        // ====================================================================================

                    case 3: {
                        // ================= version 3 hand Params without Fit =============================
                        double handParams[4] = {handConst[part][centr], handT[centr], handBeta[centr], masses[part]};
                        ifuncx[part][centr]->SetParameters(handParams);
                        break;
                    }
                        // ====================================================================================
                }
            
                ifuncx[part][centr]->SetLineColor(centrColors[centr]);
                
                double *params = ifuncx[part][centr]->GetParameters();
                const double *paramsErr = ifuncx[part][centr]->GetParErrors();

                std::copy(params, params + 4, outParams[part][centr]);
                std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);
            }
        }
    }
};



