#include "def.h"
#include "WriteReadFiles.h"

#include "TMinuit.h"

using namespace std;


void GetContourPlots(int part, int centr) {

    TMinuit* minuit = (TMinuit*)gMinuit; 
    
    if (!minuit) {
        cerr << "Minuit not initialized!" << endl;
        return;
    }

    // Параметры T (индекс 1) и beta (индекс 2)
    const int parT = 1; 
    const int parBeta = 2;

    Double_t T_val, T_err, beta_val, beta_err;
    minuit->GetParameter(parT, T_val, T_err);
    minuit->GetParameter(parBeta, beta_val, beta_err);

    // // Явный расчет ошибок, если они нулевые
    // if (T_err < 1e-10 || beta_err < 1e-10) {
    //     int status;
    //     Double_t args[1] = {0}; // Аргументы не требуются
    //     minuit->mnexcm("HESSE", args, 0, status); // Явный расчет матрицы ошибок
    //     minuit->GetParameter(parT, T_val, T_err);
    //     minuit->GetParameter(parBeta, beta_val, beta_err);
    // }

    if (T_err < 1e-10 || beta_err < 1e-10) {
        cerr << "Contour skipped: Uncertainties too small." << endl;
        return;
    }

    for (int s = 1; s < N_SIGMA; s++) {
        minuit->SetErrorDef(pow(s, 2));
        contour[part][centr][s] = (TGraph*)minuit->Contour(200, parT, parBeta);
        
        if (contour[part][centr][s]) {
            contour[part][centr][s]->SetLineColor(centrColors[centr]);
            contour[part][centr][s]->SetLineStyle(s);
            contour[part][centr][s]->SetTitle(Form("%s %s", partTitles[part].c_str(), centrTitles[centr].c_str()));
        } else {
            cerr << "Contour failed: part=" << part << " centr=" << centr << " s=" << s << endl;
        }
    }
}


class BlastWaveFit {
public:

    bool isContour = false;
    bool isDraw = true;
    double outParams[N_PARTS][N_CENTR][4];
    double outParamsErr[N_PARTS][N_CENTR][4];
    double paramsSystematics[N_PARTS][N_CENTR][4];
    double lLimitMult = 0.5, rLimitMult = 1.5; // for parLimits in case 4 (Systematic)
    double lLimitMultPi = 0.5, rLimitMultPi = 1.; // for parLimits in case 4 (Systematic Pi meson)
    
    BlastWaveFit() {
        for (int p = 0; p < N_PARTS; p++) {
            for (int c = 0; c < N_CENTR; c++) {
                for (int s = 0; s < N_SIGMA; s++) {
                    contour[p][c][s] = nullptr;
                }
            }
        }
    }



    void Fit( int initParamsType = 0 )
    {    
        // ++++++ Read data +++++++++++++++++++++++++++++++++++++

        // Чтение данных в зависимости от системы
        if (systN == 0) 
            ReadFromFileAuAu(); // Для системы AuAu
        else                    // Для других систем
            for (int part: PARTS) ReadFromFile(part, systN);

        // +++++++++ Fit +++++++++++++++++++++++++++++++++++++++

        // TVirtualFitter::SetDefaultFitter("Minuit");  
        // TMinuit* minuit = new TMinuit(5); 
        gMinuit = new TMinuit(5);  // Инициализация глобального Minuit
        gMinuit->SetPrintLevel(1); // Включить отладочный вывод
        // auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);
        auto funcx = new TF1("funcx", bwfitfunc, 0.01, 10, 5);

        // gMinuit->SetMaxIterations(1000); // Увеличьте число итераций
        // gMinuit->SetPrecision(1e-5);     // Повысьте точность

        funcx->SetParameters(2,1);
        funcx->SetParNames("constant", "T", "beta", "mass", "pt");
        MyIntegFunc *integ = new MyIntegFunc(funcx);

        for (int part: PARTS)
        {  
            for (int j = 0; j < N_CENTR_SYST[systN]; j++) {
                int centr = CENTR_SYST[systN][j];
   
                // cout << "PART: " << part << "   CENTR: " << centr << endl;
                string ifuncxName = "MyIntegFunc_" + to_string(part) + "_" + to_string(centr);
                ifuncx[part][centr] = new TF1("ifuncx", integ, xmin[part], xmax[part], 4, ifuncxName.c_str());

                switch(initParamsType)
                {
                    case 0: { /* DEFAULT */
                        // ================== version1 Params from Global fit ============================
                        double parResults[5];

                        std::string filename = "output/parameters/ALL_GlobalBWparams_" + std::string(systNamesT[systN]) + ".txt";
                        // std::string filename = "output/parameters/ALL_FinalBWparams_" + std::string(systNamesT[systN]) + ".txt";
                        ReadGlobalParams(systN, paramsGlobal, filename.c_str());
                        getGlobalParams(part, centr, parResults);
                        if (parResults[0] == 0) continue;
                            
                        // parResults[2] = (parResults[2] > 0.9) ? 0.9 : parResults[2];
                        // ifuncx[part][centr]->SetParameters(parResults);
                        // for (int par = 0; par < 3; par++)
                        // {
                        //     ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par] * 1.5);
                        //     if (part <= 1 && par == 2) 
                        //         ifuncx[part][centr]->SetParLimits(par, parResults[par] * 0.5, parResults[par] * 1.5);
                        // }

                        // Установка начальных параметров
                        ifuncx[part][centr]->SetParameters(parResults);
                        
                        // Фиксируем параметры T (индекс 1) и beta (индекс 2)
                        // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                        // ifuncx[part][centr]->FixParameter(2, parResults[2]); // Фиксируем beta
                        if (systN == 0) {
                            ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            ifuncx[part][centr]->FixParameter(2, parResults[2]); // Фиксируем beta

                            // ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.9, parResults[1] * 1.1);
                            // ifuncx[part][centr]->SetParLimits(2, parResults[2] * 0.9, parResults[2] * 1.1);

                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 1000); // Широкие границы для константы
                        } else if (systN == 1) {
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 1.09, parResults[1] * 1.3); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.09, parResults[2] * 1.3); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 300); // Широкие границы для константы
                        } else if (systN == 2) {
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 1.09, parResults[1] * 1.3); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 1.09, parResults[2] * 1.3); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 300); // Широкие границы для константы
                        } else if (systN == 3) {
                            // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.99, parResults[1] * 1.3); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 0.99, parResults[2] * 1.3); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100); // Широкие границы для константы
                        } else if (systN == 4) {
                            // ifuncx[part][centr]->FixParameter(1, parResults[1]); // Фиксируем T
                            ifuncx[part][centr]->SetParLimits(1, parResults[1] * 0.99, parResults[1] * 1.3); // T ±1%
                            ifuncx[part][centr]->SetParLimits(2, parResults[2] * 0.99, parResults[2] * 1.3); // beta ±1%
                            ifuncx[part][centr]->SetParLimits(0, parResults[0] * 0, parResults[0] * 100); // Широкие границы для константы
                        }
                        

                        // Лимиты ТОЛЬКО для константы (параметр 0)
                        

                        // Фиксируем массу (параметр 3)
                        ifuncx[part][centr]->FixParameter(3, masses[part]);

                        // Выполняем фит и получаем результат
                        TFitResultPtr fitResult = grSpectra[part][centr]->Fit(ifuncx[part][centr], "QR+S", "", xmin[part], xmax[part]);

                        // Проверяем валидность результата
                        if (fitResult->IsValid()) {
                            
                            // gMinuit->mnmigr(); // Минимизация
                            // gMinuit->mnhesse(); // Расчет матрицы ошибок
                            // gMinuit->mnerrs(parT, T_err, T_err_neg, T_err_parab, T_globalcc); // Для T
                            // gMinuit->mnerrs(parBeta, beta_err, beta_err_neg, beta_err_parab, beta_globalcc); // Для beta
                            // cout << "T_err = " << T_err << ", beta_err = " << beta_err << endl;

                            double chi2 = fitResult->Chi2();
                            int ndf = fitResult->Ndf();
                            double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;

                            std::cout << part 
                                      << centr 
                                      << " Chi2/NDF = " << chi2_ndf 
                                      << " (Chi2 = " << chi2 
                                      << ", NDF = " << ndf << ")\n" 
                                      << std::endl;
                            
                            if (isContour) { 
                              GetContourPlots(part, centr);
                            }
                        } else {
                            std::cerr << "Fit failed for part " << part << " centr " << centr << std::endl;
}
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
                        ifuncx[part][centr]->SetParLimits(1, 0.8, 0.14);	
                        ifuncx[part][centr]->SetParLimits(2, 0.4, 0.8);	
                        ifuncx[part][centr]->FixParameter(3, masses[part]);	//	mass
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }   
                        // ====================================================================================

                    case 3: {
                        // ================= version 3 hand Params without Fit =============================+==
                        double handParams[4] = {handConst[part][centr], TCuAu[centr], betaCuAu[centr], masses[part]};
                        ifuncx[part][centr]->SetParameters(handParams);
                        break;
                    }
                        // ====================================================================================
                    
                    case 4: {
                        // ================= version 4 Params For Systematics =================================    
                        paramsSystematics[part][centr][2] = (paramsSystematics[part][centr][2] > 0.95) ? 0.95 : paramsSystematics[part][centr][2];
                        ifuncx[part][centr]->SetParameters(paramsSystematics[part][centr]);
                        for (int par = 0; par < 3; par++)
                        {
                            if (paramsSystematics[part][centr][par] * rLimitMult > 0.95) rLimitMult = 0.95 / paramsSystematics[part][centr][par];
                            ifuncx[part][centr]->SetParLimits(par, paramsSystematics[part][centr][par] * lLimitMult, paramsSystematics[part][centr][par] * rLimitMult);
                            cout << paramsSystematics[part][centr][par] * lLimitMult << "   " << paramsSystematics[part][centr][par] * rLimitMult << endl;
                            // if (part <= 1 && par == 2) {
                            //     double rl = paramsSystematics[part][centr][par] * rLimitMultPi;

                            //     if (rl > 0.95) rl = 0.95;
                            //     ifuncx[part][centr]->SetParLimits(par, paramsSystematics[part][centr][par] * lLimitMultPi, rl);
                            // }
                                
                        }

                        ifuncx[part][centr]->FixParameter(3, masses[part]);
                        grSpectra[part][centr]->Fit(ifuncx[part][centr],  "QR+", "", xmin[part], xmax[part]);
                        break;
                    }
                        // ====================================================================================
                }
                
                ifuncx[part][centr]->SetLineColor(centrColors[centr]);
                
                
        // +++++++++ Metrics ++++++++++++++++++++++++++++++++++++

                double *params = ifuncx[part][centr]->GetParameters();
                const double *paramsErr = ifuncx[part][centr]->GetParErrors();
                std::copy(params, params + 4, outParams[part][centr]);
                std::copy(paramsErr, paramsErr + 4, outParamsErr[part][centr]);

                // NormalizeErrors on chi2/NDF
                double chi2 = ifuncx[part][centr]->GetChisquare();
                double ndf = ifuncx[part][centr]->GetNDF();
                double chi2Ndf = chi2 / ndf;
                for (int i = 0; i < 3; i++ )
                {
                    outParamsErr[part][centr][i] *= sqrt(chi2Ndf);
                }

                int N = 0, fitN = 0;
                double *x, *y, d = 0;
                x = grSpectra[part][centr]->GetX();
                y = grSpectra[part][centr]->GetY();
                N = grSpectra[part][centr]->GetN();

                for (int i = 0; i < N; i++)
                {
                    if (x[i] >= xmin[part] && x[i] <= xmax[part])
                    {
                        d += pow((y[i] - ifuncx[part][centr]->Eval(x[i])) / y[i], 2);
                        fitN++;
                    }
                }
                d = sqrt(d) / fitN;

                cout << part << " " << centr << "  " << d << " " << chi2Ndf << endl;
            }
        }
    }
};



