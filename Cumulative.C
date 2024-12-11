#include "math.h"
#include "TF1.h"

const double mN = 0.931494;
const double mPI = 0.13957061;
const double PI = acos(-1);

double mI = mN, mII = mN;  // - сталкиваются нуклон ядра и нуклон мишени
double m1 = mPI, m2 = 2 * mN; // - рождается 1 пи-мезон
double Minit = mI + mII; // масса покоя начальной системы сталкивающихся нуклонов
    
double maxE(double *x, double *params)
{
    double Etot, W, E1, E2, A1, pI, theta1;

    W = params[0]; // - кинетическая энергия частицы I
    //theta1 = params[1]; // угол вылета пи-мезона
    theta1 = x[0] * PI / 180.; 
    // Etot - полная энергия
 
    Etot = mI + mII + W;
    pI = sqrt(W * W + 2 * W * mI);

    A1 = Minit * Minit + 2 * mII * W + m1 * m1 - m2 * m2;
    
    double sq = A1 * A1 - 4 * m1 * m1 * (Etot * Etot - pow( pI * cos(theta1), 2));
    double E1_nom = A1 * Etot + pI * cos(theta1) * sqrt(sq);
    double E1_denom = 2 * (Etot * Etot - pow( pI * cos(theta1), 2));

    E1 = E1_nom / E1_denom;

    return E1;
}

double calculateW( double sNN )
{
    return (sNN - mI * mI - mII * mII) / 2 * mII - mI;
}

void Cumulative ( void ) 
{
    double sNN = 9.2 * 9.2;
    double W = calculateW(sNN);

    auto Efunc = new TF1("funcx", maxE, 0.1, 179, 1);
    Efunc->SetParameters(2,1);
    Efunc->FixParameter(0, W);
    Efunc->Draw();

    double theta1 = 5.7; // быстрота = 3 
    double Emax = maxE(&theta1, &W);
    double pmax = sqrt(Emax * Emax - m1 * m1);

    cout << pmax << endl;
}