#include "def.h"
#include "WriteReadFiles.h"
#include "BlastWaveFit.h"

// Hello
using namespace std;

void DrawSpectraPart( TString partName, int part1, int part2 )
{
    TCanvas *c4 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    int NpadX = 1, NpadY = 2;
    Format_Canvas(c4, NpadX, NpadY, 0);

    int padN = 1;
    
    for (int i: {part1, part2})
    {
        c4->SetLogy();
        c4->cd(padN);  

        double shiftX = (i % NpadX == 0) ? 0 : 0.1;
        double texScale = (padN == 1 ) ? 1 : 0.9;
        TLegend *legend = new TLegend(0.5 - shiftX, 0.65, 0.9 - shiftX, 0.95); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(2);
        legend->SetTextSize(0.073 * texScale);

        TLatex *titleTex = new TLatex(0.4, 2000, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.09);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int centr: CENTR)
        {
            if (!ifuncx[i][centr]) continue;

            grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            grSpectra[i][centr]->SetMarkerSize(1);
            grSpectra[i][centr]->SetMarkerStyle(8);
            grSpectra[i][centr]->Draw("P SAME");      
            ifuncx[i][centr]->Draw("SAME");
            legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "p");        
        }
        legend->Draw();
        titleTex->Draw();    
        padN++;
    }

    c4->SaveAs("output/BlastWave_" + partName + ".pdf");
}

void BlastWave( void )
{
    bool isContour = false;
    bool isDraw = true;

    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    string inputFileName = "postprocess_mpdpid10";
    SetSpectra(inputFileName, "mt");

    // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++
    BlastWaveFit *bwFit = new BlastWaveFit();
    bwFit->Fit(0);
    WriteParams(bwFit->outParams, bwFit->outParamsErr);
    WriteParams(bwFit->outParams, bwFit->outParamsErr, false);
   
    if (!isDraw)
        return;

    // ++++++ Draw spectra All +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    Format_Canvas(c2, 1, 2, 0);

    for (int i: {0, 1})
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.9;
        TLegend *legend = new TLegend(0.5 - shiftX, 0.6, 0.95 - shiftX, 0.9); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(2);
        legend->SetTextSize(0.075 * texScale);

        TLatex *titleTex = new TLatex(0.4, 500, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.09);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int centr: CENTR)
        {
            if (!ifuncx[i][centr]) continue;

            grSpectra[i][centr]->SetMarkerColor(centrColors[centr]);
            grSpectra[i][centr]->SetMarkerSize(1);
            grSpectra[i][centr]->SetMarkerStyle(8);
            grSpectra[i][centr]->Draw("P SAME");      
            ifuncx[i][centr]->Draw("SAME");
            legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "p");        
        }
        legend->Draw();
        titleTex->Draw();    
    }

    c2->SaveAs("output/BlastWave.pdf");
    delete c2;

    DrawSpectraPart("pi", 0, 1);
    DrawSpectraPart("K", 2, 3);
    DrawSpectraPart("p", 4, 5);

    if (!isContour) {
        gROOT->ProcessLine(".q");
        return;
    }

    //++++++++ Draw Contour plots ++++++++++++++++++++++++++++++

    TCanvas *c3 = new TCanvas("c3", "c3", 29, 30, 1200, 1200);
    c3->cd();
    c3->SetGrid();
    double ll = 0.00001, rl = 2, pad_min = 0.00001, pad_max = 0.2, 
        pad_offset_x = 1., pad_offset_y = 1., 
        pad_tsize = 0.05, pad_lsize=0.05;
    TString pad_title_y = "T";
    TString pad_title_x = "#beta";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "", 8);        

    TLegend *legendContour = new TLegend(0.6, 0.35, 0.85, 0.85);
    legendContour->SetBorderSize(0);
    legendContour->SetFillStyle(0);
    legendContour->SetTextSize(0.04);

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            string legendText = partTitles[part] + ", " + centrTitles[centr];
            legendContour->AddEntry(contour[part][centr][1], legendText.c_str(), "l"); 

            for (int nsigma = 1; nsigma < N_SIGMA; nsigma++)
            {
                if (contour[part][centr][nsigma])
                    contour[part][centr][nsigma]->Draw("lf");
                else 
                    cout << "ERROR " << part << " " <<centr << " " << nsigma << endl;
            }
        }
    }
    legendContour->Draw();

    c3->SaveAs("output/BlastWave_contour.pdf");
    gROOT->ProcessLine(".q");
}

