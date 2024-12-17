 
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "def.h"
#include "ReadFiles.h"

int ipar0[3] = {2, 0, 1};
int ipar2[3] = {3, 0, 1};
int ipar4[3] = {4, 0, 1};
 
bool isParamsFileExist = false;

void WriteParams( int charge, double paramsGlobal[2][N_CENTR][5], const char filename[30] = "output/GlobalBWparams.txt" )
{
   cout << " WriteParams " << endl;
   ofstream txtFile;
   if (isParamsFileExist) 
   {
      txtFile.open(filename, ios::app);
      txtFile << endl;
   }
   else txtFile.open(filename);

   for (int centr: CENTR)
   {
      txtFile  << charge << "  "  << centr << "  " 
               << paramsGlobal[charge][centr][0]  << "  " << paramsGlobal[charge][centr][1] << "  "
               << paramsGlobal[charge][centr][2] << "   " << paramsGlobal[charge][centr][3] << "   " << paramsGlobal[charge][centr][4] << endl;
   }
    
   txtFile.close();
   isParamsFileExist = true;
}

// Create the GlobalCHi2 structure
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;

   // par[0] = constant;
   // par[1] = Tf;
   // par[2] = beta;
   // par[3] = mass	pi
   // par[4] = mass	K
   // par[5] = mass	p
   double operator() (const double *par) const 
   {
      const int Nparams = 4, Nfunc = 3;
      double p[Nfunc][Nparams];
      int commonParams = 3;
      for (int i = 0; i < Nfunc; i++)
      {
         p[i][0] = par[2 + i];
         p[i][1] = par[0]; 
         p[i][2] = par[1]; 
         p[i][3] = masses[2 * i];
      }
 
      return (*fChi2_1)(p[0]) + (*fChi2_2)(p[1]) + (*fChi2_3)(p[2]);
   }
};
 
void GlobalFitCentr( double paramsInit[5], double paramsGlobal[2][N_CENTR][5], int centr, int charge = 0 ) 
{
   cout << " GlobalFitCentr " << endl;
   double xmin = 0.3, xmax = 1.2;

   // perform now global fit
   ROOT::Math::WrappedMultiTF1 wf0(*ifuncxGlobal[0 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf2(*ifuncxGlobal[2 + charge][centr], 1);
   ROOT::Math::WrappedMultiTF1 wf4(*ifuncxGlobal[4 + charge][centr], 1);
 
   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange range0, range2, range4;
   // set the data range
   
   range0.SetRange(xmin, xmax);
   ROOT::Fit::BinData data0(opt, range0);
   ROOT::Fit::FillData(data0, grSpectra[0 + charge][centr]);
 
   range2.SetRange(xmin, xmax);
   ROOT::Fit::BinData data2(opt, range2);
   ROOT::Fit::FillData(data2, grSpectra[2 + charge][centr]);
 
   range4.SetRange(xmin, xmax);
   ROOT::Fit::BinData data4(opt, range4);
   ROOT::Fit::FillData(data4, grSpectra[4 + charge][centr]);
   
   ROOT::Fit::Chi2Function chi2_0(data0, wf0);
   ROOT::Fit::Chi2Function chi2_2(data2, wf2);
   ROOT::Fit::Chi2Function chi2_4(data4, wf4);
 
   GlobalChi2 globalChi2(chi2_0, chi2_2, chi2_4);
 
   ROOT::Fit::Fitter fitter;

   const int Npar = 5;
  
   // create before the parameter settings in order to fix or set range on them
   //fitter.Config().SetParamsSettings(Npar, par0);
   fitter.Config().SetParamsSettings(Npar, paramsInit);

   fitter.Config().ParSettings(0).SetLimits(0.06, 0.13);
   fitter.Config().ParSettings(1).SetLimits(0.3, 0.95);
   fitter.Config().ParSettings(2).SetLimits(paramsInit[2 + 0 + charge] * 0.5, paramsInit[2 + 0 + charge] * 1.5);
   fitter.Config().ParSettings(3).SetLimits(paramsInit[2 + 2 + charge] * 0.5, paramsInit[2 + 2 + charge] * 1.5);
   fitter.Config().ParSettings(4).SetLimits(paramsInit[2 + 4 + charge] * 0.5, paramsInit[2 + 4 + charge] * 1.5);

   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer("Minuit2","Migrad");
 
   // fit FCN function directly
   // (specify optionally data size and flag to indicate that is a chi2 fit)
   fitter.FitFCN(5, globalChi2, 0, data0.Size() + data2.Size() + data4.Size(), true);
   ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);

   const double *fitResults = result.GetParams();
   for (int i = 0; i < 5; i++ ) paramsGlobal[charge][centr][i] = fitResults[i];

   string chargeFlag = (charge == 0) ? "pos" : "neg";
   cout<< "Result " << paramsGlobal[charge][centr][0] << "  " << paramsGlobal[charge][centr][1] << "  " << paramsGlobal[charge][centr][2] << "  " << paramsGlobal[charge][centr][3] << "  " << paramsGlobal[charge][centr][4] << endl;
}

void DrawFitSpectra( string chargeFlag = "all" )
{
   TCanvas *c2 = new TCanvas("c2", "c2", 29, 30, 1100, 1200);
   Format_Canvas(c2, 2, 3, 0);

   int padN = 1;   
   for (int part: PARTS_ALL)
   {
      c2->cd(padN++);
      FormatSpectraPad(1);

      if (chargeFlag == "pos" && part % 2 == 1) continue;
      if (chargeFlag == "neg" && part % 2 == 0) continue;

      double shiftX = (part % 2 == 0) ? 0 : 0.1;
      double texScale = (part < 3) ? 1 : 0.9;

      TLegend *legend = new TLegend(0.55 - shiftX, 0.7, 0.98 - shiftX, 0.9); //1 column
      legend->SetNColumns(2);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetTextSize(0.07 * texScale);

      TLatex *titleTex = new TLatex(0.6, 500, partTitles[part].c_str());
      titleTex->SetTextFont(42);
      titleTex->SetTextSize(0.08);
      titleTex->SetLineWidth(2 * texScale);

      double paramsGlobal[2][N_CENTR][5]; // For initial values
      ReadGlobalParams(paramsGlobal);
      
      for (int centr: CENTR)
      {   
         double parResults[5];
         getGlobalParams(part, centr, paramsGlobal, parResults);

         ifuncxGlobal[part][centr]->SetParameters(parResults);
         ifuncxGlobal[part][centr]->SetLineColor(centrColors[centr]);
         ifuncxGlobal[part][centr]->Draw("SAME");
         grSpectra[part][centr]->GetListOfFunctions()->Add(ifuncxGlobal[part][centr]);
         grSpectra[part][centr]->SetMarkerStyle(8);
         grSpectra[part][centr]->SetMarkerSize(1);
         grSpectra[part][centr]->Draw("P SAME");

         legend->AddEntry(ifuncxGlobal[part][centr], centrTitles[centr].c_str(), "l");        
      }

      legend->Draw();
      titleTex->Draw(); 
   }
 
   c2->SaveAs("output/BlastWaveGlobalFit.pdf");
   gROOT->ProcessLine(".q");
}

void BlastWaveGlobal_1Fit(string chargeFlag = "all", double setParams[5] = nullptr) 
{
   // ++++++ Read data ++++++++++++++++++++++++++++++++++++

    string inputFileName = "postprocess_mpdpid10";
    SetSpectra(inputFileName, "mt");

   // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

   double paramsGlobal[2][N_CENTR][5]; // Global params to store results

   double paramsInit[2][N_CENTR][5]; // For initial values
   ReadGlobalParams(paramsInit);

   for (int centr: CENTR)
   {
      TVirtualFitter::SetDefaultFitter("Minuit");  
      // TMinuit* minuit = new TMinuit(5); 
      TF1 *funcx;
      MyIntegFunc *integ;
      double xmin = 0.3, xmax = 1.2;

      funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
      funcx->SetParameters(2,1);
      funcx->SetParNames("constant", "T", "beta", "mass", "pt");
      integ = new MyIntegFunc(funcx);

      cout << paramsInit[0][1][0] << endl;

      for (int part: PARTS_ALL)
      {
         int charge = part % 2;
         string ifuncxName = "BW_" + to_string(part);
         ifuncxGlobal[part][centr] = new TF1("ifuncx", integ, xmin, xmax, 4, ifuncxName.c_str());
         ifuncxGlobal[part][centr]->FixParameter(3, masses[part]); 	    
         
         if (setParams)
            ifuncxGlobal[part][centr]->SetParameters(setParams);
         else
            ifuncxGlobal[part][centr]->SetParameters(paramsInit[charge][centr]);

         for (int par = 0; par < 3; par++)
         {
            if (setParams)
               ifuncxGlobal[part][centr]->SetParLimits(par, paramsInit[charge][centr][par] * 0.9, paramsInit[charge][centr][par] * 1.1);
            else 
               ifuncxGlobal[part][centr]->SetParLimits(par, paramsInit[charge][centr][par] * 0.9, paramsInit[charge][centr][par] * 1.1);
         }

         // double conMin = paramsInit[charge][centr][0] * 0.5;
         // double conMax = paramsInit[charge][centr][0] * 1.5;
         // ifuncxGlobal[part][centr]->SetParLimits(0, conMin, conMax);
         // ifuncxGlobal[part][centr]->SetParLimits(1, 0.06, 0.13);          	 //	temp.
         // ifuncxGlobal[part][centr]->SetParLimits(2, 0.3, 0.95);	          //	beta
      }

      if (chargeFlag != "neg") GlobalFitCentr(paramsInit[0][centr], paramsGlobal, centr, 0); // positive charged
      if (chargeFlag != "pos") GlobalFitCentr(paramsInit[1][centr], paramsGlobal, centr, 1); // negative charged
   }

   if (chargeFlag != "neg") WriteParams(0, paramsGlobal);
   if (chargeFlag != "pos") WriteParams(1, paramsGlobal);

   DrawFitSpectra();
}

void BlastWaveGlobal(string chargeFlag = "all")
{
   double setParams[5] = {0.1, 0.6, 300, 100, 0.1};
   for (int i = 0; i < 3; i++)
   {
      BlastWaveGlobal_1Fit(chargeFlag, setParams);   
   }
   gROOT->ProcessLine(".q");
}
