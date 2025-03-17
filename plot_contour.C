#include "input/headers/FormatOfEverything.h"

// plot_contour.C
void plot_contour() {
    // Настройки стиля
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(99);
    gStyle->SetPalette(kTemperatureMap);

    // Создаем канвас с двумя панелями
    TCanvas *c = new TCanvas("c", "Contour Plots", 1000, 800);
    c->Divide(1, 2); // 2 строки, 1 столбец

    // ================= Верхняя панель (0-5% центральность) =================
    c->cd(1);
    gPad->SetRightMargin(0.15);
    
    TH2F *h_central = new TH2F("h_central", "0-5% Central Collisions;#beta [GeV];T [GeV]",
                              100, 0.4, 0.85, 100, 0.11, 0.18);

    // Чтение данных для центральных столкновений (centrality = 0)
    std::ifstream file("output/pics/ALL_FinalBWparams_AuAu.txt");
    Double_t particle, centrality, constant, T, T_err, beta, beta_err;
    
    while (file >> particle >> centrality >> constant >> T >> T_err >> beta >> beta_err) {
        if (centrality == 0) { // 0 = 0-5% centrality
            h_central->Fill(beta, T);
        }
    }
    file.close();

    // Настройка контуров
    Double_t contours[3] = {0.68, 0.95, 0.997};
    h_central->SetContour(3, contours);
    h_central->Draw("cont3");
    
    // Лучшая точка (примерные координаты)
    TGraph *best_point = new TGraph(1);
    best_point->SetPoint(0, 0.713, 0.122); // Замените на реальные значения
    best_point->SetMarkerStyle(29);
    best_point->SetMarkerSize(2);
    best_point->Draw("P same");

    // ================= Нижняя панель (Периферийные 60-80%) =================
    c->cd(2);
    gPad->SetRightMargin(0.15);
    
    TH2F *h_peripheral = new TH2F("h_peripheral", "60-80% Peripheral Collisions;#beta [GeV];T [GeV]",
                                 100, 0.4, 0.85, 100, 0.11, 0.18);

    // Переоткрываем файл для периферийных данных
    std::ifstream file2("output/pics/ALL_FinalBWparams_AuAu.txt");
    while (file2 >> particle >> centrality >> constant >> T >> T_err >> beta >> beta_err) {
        if (centrality == 10) { // 10 = 60-80% centrality
            h_peripheral->Fill(beta, T);
        }
    }
    file2.close();

    h_peripheral->SetContour(3, contours);
    h_peripheral->Draw("cont3");
    
    // Легенда
    TLegend *leg = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg->AddEntry(best_point, "Best fit", "P");
    leg->Draw();

    // Сохранение
    c->SaveAs("contour_comparison.png");
}