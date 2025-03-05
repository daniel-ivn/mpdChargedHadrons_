#include "input/def.h"
#include "input/WriteReadFiles.h"

TGraphErrors *gr[2][5]; // [charge][system]
int SYSTS[] = {0, 1, 2, 3, 4};

// Определяем двумерный массив цветов: первая ось – заряд, вторая – система
int chargeSystColors[2][5] = {
    {kRed, kGreen, kBlue, kMagenta, kOrange},   // для charge = 0
    {kCyan, kBlack, kYellow, kPink, kViolet}       // для charge = 1
};

int markerStyles[2][5] = {
    {20, 21, 22, 23, 24},  // для charge = 0
    {25, 26, 27, 28, 29}   // для charge = 1
};

void DrawParam(TString paramName = "T")
{
    TCanvas *c2 = new TCanvas("c2", "c2", 30, 30, 1200, 1000);
    c2->cd();
    c2->SetGrid();
    c2->SetLogx();
    c2->SetLeftMargin(0.2);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.15);

    double ll = 1, rl = 1000., pad_min = 0., pad_max = (paramName == "T") ? 0.3 : 1., 
            pad_offset_x = 1.2, pad_offset_y = 1.5, 
            pad_tsize = 0.055, pad_lsize = 0.055;
    // TString pad_title_y = "T [GeV]";
    TString pad_title_x = "N_{part}";
    TString pad_title_y;
    if (paramName == "T") 
    {
        pad_title_y = "T [GeV]";
    } else if (paramName == "ut") 
    {
        pad_title_y = "u_{t} [GeV]";
    } 

    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    // Положение легенды в зависимости от параметра
    double leg_x1, leg_y1, leg_x2, leg_y2;
    if (paramName == "T") 
    {
        leg_x1 = 0.60; leg_y1 = 0.60; 
        leg_x2 = 0.90; leg_y2 = 0.90; 
    } else if (paramName == "ut") 
    {
        leg_x1 = 0.60; leg_y1 = 0.20; 
        leg_x2 = 0.90; leg_y2 = 0.50; 
    }

    TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);

    for (int systN: SYSTS)
    {
        for (int charge: {0, 1})
        {
            gr[charge][systN]->Draw("P SAME");
            legend->AddEntry(gr[charge][systN], systNamesT[systN], "P");
        }
    }

    if (paramName == "T" || paramName == "ut") 
    {
        const int nPoints = 12; 
    
        double yValues[nPoints];
    
        // Выбираем правильные значения
        if (paramName == "T") {
            for (int i = 0; i < nPoints; i++) yValues[i] = T[i];
        } else if (paramName == "ut") {
            for (int i = 0; i < nPoints; i++) yValues[i] = betaAuAu[i];
        }
    
        TGraph *theoryGraph = new TGraph(nPoints, Npart[0], yValues);
        theoryGraph->SetMarkerStyle(43); 
        theoryGraph->SetMarkerColor(kMagenta);
        theoryGraph->SetMarkerSize(1.5);
        theoryGraph->SetLineColor(kMagenta);
        theoryGraph->SetLineStyle(2); 
        theoryGraph->Draw("P SAME");
    
        legend->AddEntry(theoryGraph, "AuAu_{Th}", "P");
    }
    legend->Draw();

    double avgValue = (paramName == "T") ? gAvgT : gAvgUt;
    
    TLine *avgLine = new TLine(ll, avgValue, rl, avgValue);
    avgLine->SetLineColor(kBlack);
    avgLine->SetLineStyle(9); 
    avgLine->SetLineWidth(2);
    avgLine->Draw("SAME");

    double avgLeg_x1 = 0.23, avgLeg_y1 = 0.17; 
    double avgLeg_x2 = 0.43, avgLeg_y2 = 0.27;

    TLegend *avgLegend = new TLegend(avgLeg_x1, avgLeg_y1, avgLeg_x2, avgLeg_y2);
    avgLegend->SetBorderSize(0);
    avgLegend->SetFillStyle(0);
    avgLegend->SetTextSize(0.05);

    TString avgLabel = (paramName == "T") 
        ? Form("T_{av} = %.3f GeV", avgValue) 
        : Form("u_{t,av} = %.3f GeV", avgValue);

    avgLegend->AddEntry(avgLine, avgLabel, "L");
    avgLegend->Draw();

    cout << "test" << endl;
    TString name = "output/pics/BWparamGlobal_" + paramName + ".png";
    c2->SaveAs(name);
}


void SetGraphs( int systN, TString paramName )
{
    cout << systNamesT[systN] << endl;
    TString filename = "output/parameters/GlobalBWparams_" + systNamesT[systN] + ".txt";
    ReadGlobalParams(paramsGlobal, "output/parameters/GlobalBWparams_" + systNamesT[systN] + ".txt");
    cout << filename << "  " <<  N_CENTR_SYST[systN] << endl;
    
    for (int charge: {0, 1})
    {
        for (int centr = 0; centr < N_CENTR_SYST[systN]; centr++)
        {
            // cout << charge << " " << centr << " " << Tpar[charge][centr] << endl;

            Tpar[charge][centr] = paramsGlobal[charge][centr][0];
            utPar[charge][centr] = paramsGlobal[charge][centr][1];
        }
    }
    cout << "DEFINE GRAPHS " << N_CENTR_SYST[systN] << endl;
    
    double xerr[MAX_CENTR];
    for (int i: CENTR_SYST[systN]) xerr[i] = 0.;

    for (int charge: {0, 1})
    {
        if (paramName == "T")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], Tpar[charge], xerr, xerr); // Tpar_err[part]);
        else if (paramName == "ut")
            gr[charge][systN] = new TGraphErrors(N_CENTR_SYST[systN], Npart[systN], utPar[charge], xerr, xerr); // utPar_err[part]);
        
        // gr[charge][systN]->SetMarkerStyle(8);
        gr[charge][systN]->SetMarkerStyle(markerStyles[charge][systN]);
        gr[charge][systN]->SetMarkerSize(2);
        gr[charge][systN]->SetMarkerColor(systColors[systN]);
        // gr[charge][systN]->SetMarkerColor(chargeSystColors[charge][systN]);
    }

}

void CalculateAverage(TString paramName) {
    double sum = 0.0;
    int count = 0;

    for (int systN : SYSTS) {
        for (int charge : {0, 1}) {
            TGraphErrors* graph = gr[charge][systN];
            if (!graph) continue;
            int nPoints = graph->GetN();
            double* yValues = graph->GetY();
            
            for (int i = 0; i < nPoints; ++i) {
                sum += yValues[i];
                count++;
            }
        }
    }

    if (count > 0) {
        double average = sum / count;
        if (paramName == "T") {
            gAvgT = average;
        } else if (paramName == "ut") {
            gAvgUt = average;
        }
        cout << "Average " << paramName << " = " << average 
             << " (from " << count << " points)" << endl;
    } else {
        cout << "No data for " << paramName << endl;
    }
}

void GlobalDrawParams (void) {
    // Для параметра T
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesT[systN] << endl;
        SetGraphs(systN, "T");
    }
    CalculateAverage("T"); 
    DrawParam("T");
    
    // Для параметра ut
    for (int systN: SYSTS) 
    {
        cout << systN << " " << systNamesT[systN] << endl;
        SetGraphs(systN, "ut");
    }
    CalculateAverage("ut"); 
    DrawParam("ut");

    WriteAveragesToFile();
}
