

#include "spectra_def.h"

using namespace std;
void FormatSpectraPad( double texScale = 1 )
{
    double ll = 0.01, rl = 2.49, pad_min = 0.0011, pad_max = 3000, 
            pad_offset_x = 1., pad_offset_y = 1., 
            pad_tsize = 0.09 * texScale, pad_lsize=0.08 * texScale;
    TString pad_title_y = "d^{2}N/(p_{T}dydp_{T})";
    TString pad_title_x = "p_{T} [GeV/c]";
    Format_Pad(ll, rl, pad_min, pad_max, pad_title_x, pad_title_y, pad_offset_x, pad_offset_y, pad_tsize, pad_lsize, "");        
}

void spectra( int systN = 0 )
{
    // ++++++ Read data +++++++++++++++++++++++++++++++++++++

    string inputFileName = "postprocess_mpdpid10";
    TFile *f = new TFile(("input/" + inputFileName + ".root").c_str());
    TDirectory *fd;

    for (int i = 0; i < 6; i++)
    {
        fd = (TDirectory*)f->Get(particles[i].c_str());
        fd->cd();
        for (int centr = 0; centr < N_CENTR; centr++)
        {
            string name = "h__pt_" + particles[i] +"_centrality" + to_string(centr) + "_mc_y-0.5_0.5";
            cout << name << endl;
            hSpectra[i][centr] = (TH1D *)fd->Get(name.c_str());    
            grSpectra[i][centr] = new TGraphErrors(hSpectra[i][centr]);

            grSpectra[i][centr]->SetLineColor(centrColors[centr]);

        }
    }

    // ++++++ Draw spectra +++++++++++++++++++++++++++++++++++++

    TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1200, 1200);
    Format_Canvas(c2, 2, 3, 0);

    for (int i: PARTS)
    {
        c2->SetLogy();
        c2->cd(i + 1);  
        
        double shiftX = (i % 2 == 0) ? 0 : 0.1;
        double texScale = (i < 3) ? 1 : 0.9;

        TLegend *legend = new TLegend(0.55 - shiftX, 0.7, 0.98 - shiftX, 0.9); //1 column
        legend->SetNColumns(2);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.07 * texScale);

        TLatex *titleTex = new TLatex(0.6, 500, partTitles[i].c_str());
        titleTex->SetTextFont(42);
        titleTex->SetTextSize(0.08);
        titleTex->SetLineWidth(2 * texScale);

        FormatSpectraPad(texScale);
        for (int centr: CENTR)
        {
            grSpectra[i][centr]->Draw("SAME");      
            legend->AddEntry(grSpectra[i][centr], centrTitles[centr].c_str(), "l");        
        }
        legend->Draw();
        titleTex->Draw();    
    }

    c2->SaveAs(("output/spectra_" + inputFileName + ".pdf").c_str());
    delete c2;
}

