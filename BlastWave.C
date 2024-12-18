#include "def.h"
#include "ReadFiles.h"

using namespace std;

void GetContourPlots( int part, int centr )
{
    for (int s = 1; s < N_SIGMA; s++)
    {
        gMinuit->SetErrorDef(s * 4); //note 4 and not 2!
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

void WriteParams( const char filename[30] = "output/BWparams.txt" )
{
    ofstream txtFile;
    txtFile.open(filename);

    for (int part: PARTS)
    {
        for (int centr: CENTR)
        {
            txtFile << part << "  " << centr << "  " << constPar[part][centr] << "  "
                    << Tpar[part][centr] << "   " << Tpar_err[part][centr] << "   "
                    << utPar[part][centr] << "   " << utPar_err[part][centr] << endl;
        }
    }
    
    txtFile.close();
}

void BlastWave( void )
{
    bool isContour = false;
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

            // ================== version1 Params from Global fit ============================
            // double parResults[5];
            // getGlobalParams(part, centr, parResults);
            // if (parResults[0] == 0)
            //     continue;
                
            // ifuncx[part][centr]->SetParameters(parResults);
            // for (int par = 0; par < 3; par++)
            // {
            //     ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.9, parResults[par] * 1.1);
            // }

            // ================= version 2 Params with limits =================================
            // double customParams[4] = {con[part], 0.09, 0.75, masses[part]};
            // ifuncx[part][centr]->SetParameters(customParams);
            // ifuncx[part][centr]->SetParLimits(0, conmin[part], conmax[part]);
            // ifuncx[part][centr]->SetParLimits(1, 0.8, 0.1);	
            // ifuncx[part][centr]->SetParLimits(2, 0.5, 0.8);	
            // ifuncx[part][centr]->FixParameter(3, masses[part]);	//	mass

            // ================= version 3 hand Params without Fit =============================
            double handParams[4] = {handConst[part][centr], handT[centr], handBeta[centr], masses[part]};
            ifuncx[part][centr]->SetParameters(handParams);

            ifuncx[part][centr]->SetLineColor(centrColors[centr]);
            // grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[centr], xmax[centr]);

            // constPar[part][centr] = ifuncx[part][centr]->GetParameter(0);
            // Tpar[part][centr] = ifuncx[part][centr]->GetParameter(1);
            // utPar[part][centr] = ifuncx[part][centr]->GetParameter(2);
            // Tpar_err[part][centr] = ifuncx[part][centr]->GetParError(1);
            // utPar_err[part][centr] = ifuncx[part][centr]->GetParError(2);

            if (isContour) GetContourPlots(part,  centr);    
        }
    }

    // WriteParams();

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

    if (!isContour) return;

    // //++++++++ Draw Contour plots ++++++++++++++++++++++++++++++

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

}

