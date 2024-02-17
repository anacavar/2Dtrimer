// Varijacijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nw 10      // broj šetača
#define Nt 10      // broj koraka
#define Nb 110     // broj blokova
#define NbSkip 0   // broj prvih blokova koje preskačemo
#define sigma 4    // angstrema
#define epsilon 12 // dubina jame, u kelvinima preko boltzmannove konstante
#define L0 100.    // angstrema
#define alpha 4.16 // čega?
#define gamma 2.82 // čega?
#define s 0.0027   // čega?
#define mass 4.    // u
#define hbar 1.    // koliko?
#define k_B 1.     // koliko?

double E_kin_L(double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double); // kinetički dio lokalne energije
double U_LJ(double);
double f_ddr(double), f_dr(double);

int main(void)
{
#pragma region // KONSTANTE & VARIJABLE
    long idum = -1234;
    int ib, it, iw, k;
    double x[4][Nw + 1], y[4][Nw + 1]; // indeks čestice, indeks šetača
    double x_old[4], y_old[4];
    double dx, dy;            // promjene koordinata čestica
    double dxyMax = L0 / 100; // maksimalne promjene x,y,z koordinata
    double random_number;
    double Psi(double);            // probna valna funkcija
    double psi_initial, psi_final; // ukupna inicijalna i finalna probna valna funkcija
    double P[Nw + 1];              // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
    double E_L[Nw + 1];            // lokalna energija zadnjeg koraka svakog od šetača
    double r12, r23, r13;          // udaljenosti između čestica
    double T;                      // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    double ratio;
    double SwE;        // = suma(srednjih E) po setacima
    double StE;        // = suma (srednjih E) po koracima
    double SbE;        // = suma (srednjih E) po blokovima
    double SbE2;       // = suma (srednjih E) po blokovima
    double AE, sigmaE; // srednja vrijednost i standardna devijacija
    int NbEff;         // efektivni indeks bloka
    int itmp;          // postotak prihvaćanja
#pragma endregion

    FILE *data;
    data = fopen("data.txt", "w");

    // inicijalizacija koordinata čestica gdje je gustoća Psi*Psi znacajna
    for (iw = 1; iw <= Nw; iw++) // po šetačima
    {
        for (k = 1; k <= 3; k++) // po česticama
        {
            x[k][iw] = (2 * ran1(&idum) - 1) * L0;
            y[k][iw] = (2 * ran1(&idum) - 1) * L0;
        }
        r12 = sqrt(pow((x[1][iw] - x[2][iw]), 2) + pow((y[1][iw] - y[2][iw]), 2));
        r23 = sqrt(pow((x[2][iw] - x[3][iw]), 2) + pow((y[2][iw] - y[3][iw]), 2));
        r13 = sqrt(pow((x[1][iw] - x[3][iw]), 2) + pow((y[1][iw] - y[3][iw]), 2));
        P[iw] = Psi(r12) * Psi(r23) * Psi(r13);
    }

    SbE = 0.;
    SbE2 = 0;
    for (ib = 1; ib <= Nb; ib++)
    { // blokovi
        StE = 0;
        NbEff = ib - NbSkip;
        for (it = 1; it <= Nt; it++) // po koracima
        {
            SwE = 0;
            for (iw = 1; iw <= Nw; iw++) // po šetačima
            {
                for (k = 1; k <= 3; k++) // po česticama
                {
                    dx = (ran1(&idum) * 2 - 1) * dxyMax;
                    dy = (ran1(&idum) * 2 - 1) * dxyMax;
                    x_old[k] = x[k][iw];
                    y_old[k] = y[k][iw];
                    x[k][iw] += dx;
                    y[k][iw] += dy;
                }
                // METROPOLIS ALGORITAM
                // T = |P(R_final)|^2/|P(R_initial)|^2; za svakog šetača posebno
                psi_initial = P[iw];
                r12 = sqrt(pow((x[1][iw] - x[2][iw]), 2) + pow((y[1][iw] - y[2][iw]), 2));
                r23 = sqrt(pow((x[2][iw] - x[3][iw]), 2) + pow((y[2][iw] - y[3][iw]), 2));
                r13 = sqrt(pow((x[1][iw] - x[3][iw]), 2) + pow((y[1][iw] - y[3][iw]), 2));
                psi_final = Psi(r12) * Psi(r23) * Psi(r13);
                T = psi_final * psi_final / (psi_initial * psi_initial);
                // printf("Psi_initial: %f\tPsi_final: %f => T=psi_f**2/psi_i**2= %f\n", psi_initial, psi_final, T);
                if (T > 1) // prihvaćamo pomak
                {
                    accepted++;
                    P[iw] = psi_final;
                    // printf("accepted, T: %f\n", T);
                }
                else
                {
                    random_number = ran1(&idum);
                    if (random_number <= T) // prihvaćamo pomak
                    {
                        accepted++;
                        P[iw] = psi_final;
                        // printf("accepted, T: %f\tRand: %f\n", T, random_number);
                    }
                    else // odbacujemo pomak
                    {
                        rejected++;
                        // printf("rejected, T: %f\tRand: %f\n", T, random_number);
                        for (k = 1; k <= 3; k++) // po česticama
                        {
                            x[k][iw] = x_old[k];
                            y[k][iw] = y_old[k];
                        }
                    }
                }
                // računamo Lokalnu energiju - za svakog šetača, nakon što smo prihvatili (ili odbacili) korak
                E_L[iw] = E_kin_L(r12, r13, r23) + E_pot_L(r12, r13, r23); // kinetički dio + potencijalni dio
                SwE = SwE + E_L[iw];
            } // kraj petlje šetača
            // akumulacija podataka nakon stabilizacije
            if (ib > NbSkip)
            {
                StE += SwE / Nw;
            }
        } // kraj petlje koraka
        // nakon svakog bloka podeđavamo maksimalnu duljinu koraka kako bi prihvaćanje bilo oko 50%
        ratio = (double)accepted / (double)(accepted + rejected);
        // printf("Prihvacenih: %d\nOdbijenih: %d\nOmjer: %f\n", accepted, rejected, ratio);

        if (ratio > 0.5)
            dxyMax = dxyMax * 1.05;
        if (ratio < 0.5)
            dxyMax = dxyMax * 0.95;

        // printf("dxyMax = %f\n", dxyMax);

        if (ib > NbSkip)
        {
            SbE += StE / Nt;
            SbE2 += StE * StE / (Nt * Nt);
            fprintf(data, "%d\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff); // indeks bloka, srednji E po koracima, Srednji E po blokovima
            itmp = (int)(round(ratio * 100.));
            printf("%6d. blok:  %d%% prihvacenih,  Eb = %10.2e\n", NbEff, itmp, StE / Nt);
        }
    } // kraj petlje blokova

    AE = SbE / NbEff;
    sigmaE = sqrt((SbE2 / NbEff - AE * AE) / (NbEff - 1.));
    printf("\n konacni max. korak: %6.2e\n", dxyMax);
    printf("\n E = %8.5e +- %6.2e \n\n", AE, sigmaE);
    fclose(data);
    return 0;
}

// probna valna funkcija
double Psi(double r)
{
    return exp(-pow(alpha / r, gamma) - s * r) / sqrt(r);
}

double E_pot_L(double r12, double r13, double r23)
{
    return U_LJ(r12) + U_LJ(r13) + U_LJ(r23);
}

double U_LJ(double r)
{

    return 4 * epsilon * (pow((sigma / r), 12) - pow((sigma / r), 6));
}

double E_kin_L(double r12, double r13, double r23)
{
    double D = pow(hbar, 2) / (2 * mass);
    return -D * (f_ddr(r12) + f_ddr(r13) + pow((f_dr(r12) + f_dr(r13)), 2) + f_ddr(r12) + f_ddr(r23) + pow((f_dr(r12) + f_dr(r23)), 2) + f_ddr(r13) + f_ddr(r23) + pow((f_dr(r13) + f_dr(r23)), 2));
}

double f_dr(double r)
{
    return 1 / pow(r, 2) * (gamma * pow((alpha / r), gamma) - s - 1 / 2);
}

double f_ddr(double r)
{
    return 1 / pow(r, 2) * (gamma * (gamma + 1) * pow((alpha / r), gamma) - s - 1 / 2);
}
