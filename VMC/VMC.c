// Varijacijski Monte Carlo
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

// konstante u zadanim mjernim jedinicama i njihov preračun u SI
#define k_B 1.                                         // boltzmannova konstanta
#define k_B_si 1.3806 * pow(10, -23)                   // J/K = m2 kg s-2 K-1
#define A 1.                                           // angstrem
#define A_si 1.0 * pow(10, -10)                        // m
#define K 1.                                           // kelvin
#define K_si 1.                                        // kelvin
#define u 1.                                           // jedinica mase
#define u_si 1.6605 * pow(10, -27)                     // kg
#define hbar2 4.851159708 * 10 * k_B *K *pow(A, 2) * u // (reducirana planckova konstanta)^2
#define hbar2_si 1.112121717 * pow(10, -68)            // (Js)^2 = (m2kg/s)^2
// početne vrijednosti
#define Nt 100              // broj koraka
#define Nw 100              // broj šetača
#define Nb 250              // broj blokova
#define NbSkip 0            // broj prvih blokova koje preskačemo
#define sigma 4 * A         // angstrema
#define epsilon 12 * k_B *K // dubina jame, u kelvinima preko boltzmannove konstante
#define L0 30. * A          // angstrema
#define alpha 4.16 * A      // angstrema
#define gamma 2.82          // eksponent u probnoj valnoj funkciji
#define s 0.0027            // eksponent u probnoj valnoj funkciji
#define mass 4. * u         // u

// deklaracija funkcija
double Psi(double);                                                                     // probna valna funkcija (korelacijska funkcija)
double U_LJ(double);                                                                    // Lennard-Jonesov potencijal
double E_kin_L(double, double, double, double, double, double, double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double);                                                 // potencijalni dio lokalne energije
double f_ddr(double), f_dr(double);

// OČEKIVANA ENERGIJA BI TREBALA BITI OKO 5 KELVINA
int main(void)
{
#pragma region // VARIJABLE
    long idum = -1234;
    int ib, it, iw, k;
    double x[4][Nw + 1], y[4][Nw + 1]; // indeks čestice, indeks šetača
    double x_old[4], y_old[4];
    double dx, dy;            // promjene koordinata čestica
    double dxyMax = L0 / 100; // maksimalne promjene x,y,z koordinata
    double random_number;
    double psi_initial, psi_final; // ukupna inicijalna i finalna probna valna funkcija
    double P[Nw + 1];              // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
    double E_L[Nw + 1];            // lokalna energija zadnjeg koraka svakog od šetača
    double r12, r23, r13;          // udaljenosti između čestica
    double x1, x2, x3, y1, y2, y3; // trenutne koordinate svih čestica
    double T;                      // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    double ratio;
    double SwE;        // = suma(srednjih E) po setacima
    double StE;        // = suma (srednjih E) po koracima
    double SbE;        // = suma (srednjih E) po blokovima
    double SbE2;       // = suma (srednjih E^2) po blokovima
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
    for (ib = 1; ib <= Nb; ib++) // po blokovima
    {
        StE = 0;
        NbEff = ib - NbSkip;
        // resetirat accepted, rejected ovdje?
        accepted = 0;
        rejected = 0;
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
                r12 = sqrt(pow((x[1][iw] - x[2][iw]), 2) + pow((y[1][iw] - y[2][iw]), 2));
                r23 = sqrt(pow((x[2][iw] - x[3][iw]), 2) + pow((y[2][iw] - y[3][iw]), 2));
                r13 = sqrt(pow((x[1][iw] - x[3][iw]), 2) + pow((y[1][iw] - y[3][iw]), 2));
                x1 = x[1][iw];
                x2 = x[2][iw];
                x3 = x[3][iw];
                y1 = y[1][iw];
                y2 = y[2][iw];
                y3 = y[3][iw];
                E_L[iw] = E_kin_L(r12, r13, r23, x1, x2, x3, y1, y2, y3) + E_pot_L(r12, r13, r23); // kinetički dio + potencijalni dio
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
            fprintf(data, "%d\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff); // indeks bloka, srednji E po koracima (po jednom bloku), Srednji E po blokovima (od početka simulacije)
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

double E_kin_L(double r12, double r13, double r23, double x1, double x2, double x3, double y1, double y2, double y3)
{
    double D = hbar2 / (2 * mass);
    double dio1 = f_ddr(r12) + f_ddr(r13) + pow(f_dr(r12) * (x2 - x1) + f_dr(r13) * (x3 - x1), 2) + pow(f_dr(r12) * (y2 - y1) + f_dr(r13) * (y3 - y1), 2);
    double dio2 = f_ddr(r12) + f_ddr(r23) + pow(f_dr(r12) * (x1 - x2) + f_dr(r23) * (x3 - x2), 2) + pow(f_dr(r12) * (y1 - y2) + f_dr(r23) * (y3 - y2), 2);
    double dio3 = f_ddr(r13) + f_ddr(r23) + pow(f_dr(r13) * (x1 - x3) + f_dr(r23) * (x2 - x3), 2) + pow(f_dr(r13) * (y1 - y3) + f_dr(r23) * (y2 - y3), 2);
    return -D * (dio1 + dio2 + dio3);
}

double f_dr(double r)
{
    double fdr = 1 / pow(r, 2) * (gamma * pow((alpha / r), gamma) - s * r - 1 / 2);
    return fdr;
}

double f_ddr(double r)
{
    double fddr = -1 / pow(r, 2) * (pow(gamma, 2) * pow((alpha / r), gamma) + s * r);
    return fddr;
}
