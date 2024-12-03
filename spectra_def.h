#ifndef __DEF_H_
#define __DEF_H_

#include "input/FormatOfEverything.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <TH1.h>
#include <TFile.h>
#include "TLegend.h"


const int N_CENTR = 6;
const int N_PARTS = 6;
const int PARTS[N_PARTS] = {0, 1, 2, 3, 4, 5};
const int CENTR[N_CENTR] = {0, 1, 2, 3, 4, 5};

TH1D *hSpectra[6][12];
TGraphErrors *grSpectra[6][12];

string particles[6] = {"pip", "pim", "kp", "km", "p", "ap"};
string partTitles[6] = { "#pi^{+}","#pi^{#minus}","K^{+}","K^{#minus}","p", "#bar{p}"};
double masses[6] = {0.13957061, 0.13957061, 0.493667, 0.493667, 0.938272, 0.938272};
Color_t centrColors[11] = {kRed, kBlue, kGreen + 2, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack};
Color_t partColors[6] = {kRed, kRed, kBlue, kBlue, kGreen + 2, kGreen + 2};
string centrTitles[10] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-60%", "60-80%"};
// string centrTitles[10] = {"0-20%", "20-40%", "40-60%", "60-80%", "40-50%", "50-60%", "60-70%", "70-80%"};

#endif /* __DEF_H_ */