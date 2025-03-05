#ifndef __DEF_H_
#define __DEF_H_

#include "FormatOfEverything.h"
#include "BlastWave.h"

double gAvgT;
double gAvgUt;
TH1D *hSpectra[6][12];
TGraphErrors *grSpectra[6][12];

// Названия систем столкновений
TString systNamesT[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
string systNames[5] = {"AuAu", "pAl", "HeAu", "CuAu", "UU"};
int Ncentr[5] = {12, 4, 5, 5, 4};

const int MAX_CENTR = 20;                   // Максимальное число центральностей
const int MAX_PARTS = 6;                    // Максимальное число частиц
const int N_CENTR = 12;                      // Текущее число центральностей
const int N_PARTS = 6;                      // Текущее число частиц
const int N_SIGMA = 7;                      // Число сигм для контуров
const int PARTS[] = {0, 1, 2, 3, 4, 5};     // Индексы всех частиц
const int PARTS_POS[] = {0, 2, 4};          // Положительные частицы (π⁺, K⁺, p)
const int PARTS_NEG[] = {1, 3, 5};          // Отрицательные частицы (π⁻, K⁻, анти-p)
const int PARTS_ALL[] = {0, 1, 2, 3, 4, 5};

const int CENTR[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};          // Индексы центральностей
// const int CENTR[] = {0, 1, 2, 3}; 
const int N_CENTR_SYST[5] = {12, 4, 5, 5, 4}; // Число центральностей для каждой системы 
const int CENTR_SYST[5][20] = {               // Конкретные центральности для систем
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, // AuAu
        {0, 1, 2, 3},                           // pAl
        {0, 1, 2, 3, 4},                        // HeAu
        {0, 1, 2, 3, 4},                        // CuAu
        {0, 1, 2, 3},                           // UU
};

string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
string partTitles[6] = {"#pi^{+}", "#pi^{#minus}", "K^{+}", "K^{#minus}", "p", "#bar{p}"};
double masses[6] = {0.13957061, 0.13957061, 0.493667, 0.493667, 0.938272, 0.938272};
Color_t centrColors[11] = {kRed, kBlue, kGreen + 2, kBlack, kMagenta, kBlue+2, kBlack, kBlack, kBlack, kBlack, kBlack};
Color_t partColors[6] = {kRed, kRed, kBlue, kBlue, kGreen + 2, kGreen + 2};
string centrTitles[10] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-60%", "60-80%"};
double centrX[10] = {5, 15, 25, 35, 45, 70};

string centrTitlesAuAu[MAX_CENTR] = {"Minimum bias", "0-5%", "5-10%", "10-15%", "15-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-92%"};
double TAuAu[MAX_CENTR] = {132, 107.8, 109.8, 113.3, 116.5, 123, 132, 142, 153, 163, 168, 179};


TString NpartStr[5][MAX_CENTR] =
    {
        {"3.1", "4.35", "3.3", "2.7"},
        {"11.34", "21.84", "15.38", "9.51", "4.87"},
        {"57.0", "154.8", "80.4", "34.9", "7.5"},
        {"330", "159.0", "61.6", "17.8"},
        {"109.1", "351.4", "325.2", "299.0", 
         "253.9", "234.6", "215.3", "166.6", 
         "114.2", "74.4 ", "45.5 ", "25.7 ", 
         "19.5 ", "14.5 ", "13.4 ", "9.5 ", "6.3 ", "14.5 "}
    };
TString centrStr[5][MAX_CENTR] =
    {
        {"0-72%", "0-20%", "20-40%", "40-72%"},
        {"0-88%", "0-20%", "20-40%", "40-60%", "60-88%"},
        {"0-80%", "0-20%", "20-40%", "40-60%", "60-80%"},
        {"0-80%", "0-20%", "20-40%", "40-60%", "60-80%"},
        {"MB", "0-10", "5-10", "10-15", 
         "10-20", "15-20", "20-30", "30-40", 
         "40-50", "50-60", "60-70", "60-80", 
         "60-92", "70-80", "70-92", "80-92", "60-92"}
    };



// =============== Для  BlastWave ======================
double xmin[] = {0.5, 0.5, 0.12, 0.4, 0.2, 0.12};
double xmax[] = {1., 1., 1, 1, 1, 1};

// double xmin[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
// double xmax[] = {4, 4, 4, 4, 4, 4, };

TGraph *contour[MAX_PARTS][N_CENTR][N_SIGMA];
TF1 *ifuncx[MAX_PARTS][N_CENTR], *ifuncxGlobal[MAX_PARTS][N_CENTR];

double paramsGlobal[2][N_CENTR][5]; // [2] - charge, 5 - количество параметров 1) T 2) ut 3) const pi 4) const K 5) const p

// По центральностям
double T[] = {0.132,  0.1078, 0.1098, 0.1133, 0.1165, 0.123, 0.132, 0.142, 0.153, 0.163, 0,168, 0.179};
double Tmin[] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
double Tmax[] = {0.22, 0.22, 0.22, 0.22, 0.22, 0.22};

double betaAuAu[] = {0.71, 0.773, 0.769, 0.763, 0.754, 0.738, 0.71, 0.67, 0.614, 0.555, 0.497, 0.399};
double beta_[MAX_CENTR] = {0.673, 0.769, 0.614, 0.497, 0.399};
double betamin[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
double betamax[] = {0.9, 0.9, 0.8, 0.9, 0.9, 0.9};

// По частицам
// AuAu
// double con[6] = {100.0, 100.0, 50.0, 5000, 5000, 5000};
// double conmin[] = {0, 0, 0, 3000, 3000, 3000};
// double conmax[] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000};
// pAl, HeAu, CuAu
double con[]    = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000};  
double conmin[] = {0, 0, 0, 0, 0, 0};
double conmax[] = {30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000};

// Для GlobalFit
// AuAu
// double conGlobal[]    = {0, 0, 0, 0, 0, 0};  
// double conminGlobal[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
// double conmaxGlobal[] = {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000};

// pAl, HeAu, CuAu
double conGlobal[]    = {100, 100, 120, 60, 0.01, 0.01};  
double conminGlobal[] = {0, 0, 0, 0, 0, 0};
double conmaxGlobal[] = {5000, 5000, 5000, 5000, 5000, 5000};

// Поиск констант вручную
double handT[] = {132, 107.8, 109.8, 113.3, 116.5, 123, 132, 142, 153, 163, 168, 179};
double handBeta[] = {0.71, 0.773, 0.769, 0.763, 0.754, 0.738, 0.71, 0.67, 0.614, 0.555, 0.497, 0.399};
double handConst[MAX_PARTS][MAX_CENTR] = 
    {
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
        {10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
    };


// Извлечение глобальных параметры модели 
void getGlobalParams( int part, int centr, double parResults[4] )
{
    int charge = part % 2; 

    parResults[0] = paramsGlobal[charge][centr][2 + part / 2];
    parResults[1] = paramsGlobal[charge][centr][0]; 
    parResults[2] = paramsGlobal[charge][centr][1]; 
    parResults[3] = masses[part];     
}

double Npart[5][MAX_CENTR] =
    {   // systN, centr
        {109.1, 351.4, 299.0, 253.9, 
         215.3, 166.6, 114.2, 74.4 , 
         45.5 , 25.7 , 13.4 , 6.3  },      // AuAu
        {3.1, 4.35, 3.3, 2.7, 0},          // pAl
        {11.34, 21.84, 15.38, 9.51, 4.87}, // HeAu
        {57.0, 154.8, 80.4, 34.9, 7.5},    // CuAu
        {330, 159.0, 61.6, -17.8},         // UU
    };

double constPar[N_PARTS][N_CENTR],
    Tpar[N_PARTS][N_CENTR], Tpar_err[N_PARTS][N_CENTR], 
    utPar[N_PARTS][N_CENTR], utPar_err[N_PARTS][N_CENTR];

double paramsGlobalAllParts[N_CENTR][5];

Color_t systColors[6] = {kBlack, kBlue, kGreen + 2, kRed + 2, kMagenta};

// Вычисление поперечной массы частицы
double GetMt( int part, double pT )
{
    return sqrt(pT * pT  + masses[part] * masses[part]) - masses[part];
}


//=======================================================

void SetSpectra(string inputFileName = "postprocess_mpdpid10", string type = "pt")
{
    TFile *f = new TFile(("input/" + inputFileName + ".root").c_str());
    TDirectory *fd;

    for (int i = 0; i < 6; i++)
    {
        fd = (TDirectory*)f->Get(particles[i].c_str());
        fd->cd();
        for (int centr = 0; centr < N_CENTR; centr++)
        {
            string name = "h__pt_" + particles[i] +"_centrality" + to_string(centr) + "_mc_y-0.5_0.5";
            hSpectra[i][centr] = (TH1D *)fd->Get(name.c_str());    
            
            if (!hSpectra[i][centr]) continue;
            const int N_BINS = hSpectra[i][centr]->GetNbinsX();
            double mT[N_BINS], pT[N_BINS], sp[N_BINS], sp_err[N_BINS], xerr[N_BINS];
            for (int bin = 1; bin < N_BINS; bin++)
            {
                sp[bin - 1] = hSpectra[i][centr]->GetBinContent(bin);
                sp_err[bin - 1] = hSpectra[i][centr]->GetBinError(bin);
                xerr[bin - 1] = 0.;
                pT[bin - 1] = hSpectra[i][centr]->GetBinCenter(bin);
                mT[bin - 1] = sqrt(pT[bin - 1] * pT[bin - 1] + masses[i] * masses[i]) - masses[i];
                // cout << i << "  " << pT[bin - 1] << "  " << sp_err[bin - 1] << endl;
            }
    
            if (type == "mt") 
                grSpectra[i][centr] = new TGraphErrors(N_BINS - 1, mT, sp, xerr, sp_err);
            else if (type == "pt")
                grSpectra[i][centr] = new TGraphErrors(hSpectra[i][centr]);

            grSpectra[i][centr]->SetLineColor(centrColors[centr]);
        }
    }
}

#endif /* __DEF_H_ */