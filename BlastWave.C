#include "def.h"
#include "WriteReadFiles.h"
#include "BlastWaveFit.h"

using namespace std;

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
   
    if (!isDraw)
        return;

    // ++++++ Draw spectra +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    Format_Canvas(c2, 2, 3, 0);

    for (int i: PARTS)
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.9;
        TLegend *legend = new TLegend(0.5 - shiftX, 0.7, 0.95 - shiftX, 0.9); //1 column
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetNColumns(2);
        legend->SetTextSize(0.072 * texScale);

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

