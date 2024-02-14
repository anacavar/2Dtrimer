// Varijacijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nw 1     // broj šetača
#define Nk 10    // broj koraka
#define Nb 50    // broj blokova
#define NbSkip 0 // broj prvih blokova koje preskačemo

double E_kin_L(double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double); // kinetički dio lokalne energije
double U_LJ(double);
double f_ddr(double), f_dr(double);

int main(void)
{
#pragma region // KONSTANTE & VARIJABLE
    long idum = -1234;
    int i, j, k, ib;
    double x[3][Nw];         // indeks čestice, indeks šetača
    double y[3][Nw];         // indeks čestice, indeks šetača
    double L0 = 10;          // angstromi (10^(−10) m))
    double dx, dy;           // promjene koordinata čestica
    double dxyMax = L0 / 10; // maksimalne promjene x,y,z koordinata
    double x_old[3], y_old[3];
    double random_number;
    double Psi(double);            // probna valna funkcija
    double P[Nw];                  // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
    double psi_initial, psi_final; // ukupna inicijalna i finalna probna valna funkcija
    double r12, r23, r13;          // udaljenosti između čestica
    double T;                      // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    double ratio;
    double SwE;  // = suma(srednjih E) po setacima
    double SkE;  // = suma (srednjih E) po koracima
    double SbE;  // = suma (srednjih E) po blokovima
    double SbE2; // = suma (srednjih E) po blokovima
    int NbEff;
    int itmp;       // ŠTO JE OVO???
    double E_L[Nw]; // lokalna energija zadnjeg koraka svakog od šetača

#pragma endregion

    FILE *data;
    data = fopen("data.txt", "w");
    // inicijalizacija koordinata čestica gdje je gustoća Psi*Psi znacajna
    for (i = 0; i < Nw; i++) // po šetačima
    {
        for (j = 0; j < 3; j++) // po česticama
        {
            x[j][i] = (2 * ran1(&idum) - 1) * L0;
            y[j][i] = (2 * ran1(&idum) - 1) * L0;
        }
        r12 = sqrt(pow((x[0][i] - x[1][i]), 2) + pow((y[0][i] - y[1][i]), 2));
        r23 = sqrt(pow((x[1][i] - x[2][i]), 2) + pow((y[1][i] - y[2][i]), 2));
        r13 = sqrt(pow((x[0][i] - x[2][i]), 2) + pow((y[0][i] - y[2][i]), 2));
        P[i] = Psi(r12) * Psi(r23) * Psi(r13);
    }

    SbE = 0.;
    SbE2 = 0;
    for (ib = 0; ib < Nb; ib++)
    { // blokovi
        SkE = 0;
        NbEff = ib - NbSkip;
        for (i = 0; i < Nk; i++) // po koracima
        {
            SwE = 0;
            for (j = 0; j < Nw; j++) // po šetačima
            {
                for (k = 0; k < 3; k++) // po česticama
                {
                    dx = (ran1(&idum) * 2 - 1) * dxyMax;
                    dy = (ran1(&idum) * 2 - 1) * dxyMax;
                    x_old[k] = x[k][j];
                    y_old[k] = y[k][j];
                    x[k][j] += dx;
                    y[k][j] += dy;
                }
                // METROPOLIS ALGORITAM
                // T = |P(Rf)|^2/|P(Ri)|^2; za svakog šetača posebno
                psi_initial = P[j];
                r12 = sqrt(pow((x[0][j] - x[1][j]), 2) + pow((y[0][j] - y[1][j]), 2));
                r23 = sqrt(pow((x[1][j] - x[2][j]), 2) + pow((y[1][j] - y[2][j]), 2));
                r13 = sqrt(pow((x[0][j] - x[2][j]), 2) + pow((y[0][j] - y[2][j]), 2));
                psi_final = Psi(r12) * Psi(r23) * Psi(r13);
                T = psi_final * psi_final / (psi_initial * psi_initial);
                // T = psi_final * psi_fin / psi_initial * psi_initial;
                // printf("r12: %f\tr23: %f\tr13: %f\n", r12, r23, r13);
                // printf("Psi(R12): %f\n", Psi(r12));
                printf("Psi_initial: %f\tPsi_final: %f => T=psi_f2/psi_i2= %f\n", psi_initial, psi_final, T);
                if (T > 1) // prihvaćamo pomak
                {
                    accepted++;
                    P[j] = psi_final;
                    printf("accepted, T: %f\n", T);
                }
                else
                {
                    random_number = ran1(&idum);
                    if (random_number <= T) // prihvaćamo pomak
                    {
                        accepted++;
                        P[j] = psi_final;
                        printf("accepted, T: %f\tRand: %f\n", T, random_number);
                    }
                    else // odbacujemo pomak
                    {
                        rejected++;
                        printf("rejected, T: %f\tRand: %f\n", T, random_number);
                        for (k = 0; k < 3; k++) // po česticama
                        {
                            x[k][j] = x_old[k];
                            y[k][j] = y_old[k];
                        }
                    }
                }
                // ovdje sad mislim računamo Lokalnu energiju - za svakog šetača, nakon što smo prihvatili (ili odbacili) korak
                E_L[j] = E_kin_L(r12, r13, r23);
                SwE = SwE + E_L[j]; // prije kraja petlje šetača
            }                       // kraj petlje šetača
            // akumulacija podataka nakon stabilizacije  - zašto je ovo ovdje? - jer je ovo nakon kraja šetača za svaki korak nakon početnih skippaninh
            if (ib > NbSkip)
            {
                SkE += SwE / Nw;
            }
        } // kraj petlje koraka
        // maksimalnu duljinu koraka podešavamo kako bi prihvaćanje bilo oko 50% - dakle tek nakon svakog bloka ovo radimo
        ratio = (double)accepted / (double)(accepted + rejected);
        printf("Prihvacenih: %d\nOdbijenih: %d\nOmjer: %f\n", accepted, rejected, ratio);

        if (ratio > 0.5)
            dxyMax = dxyMax * 1.05;
        if (ratio < 0.5)
            dxyMax = dxyMax * 0.95;

        printf("dxyMax = %f\n", dxyMax);

        if (ib >= NbSkip)
        {
            SbE += SkE / Nk;
            SbE2 += SkE * SkE / (Nk * Nk);
            // fprintf(fout, "%7d%16.8e%16.8e\n", NbEff, SkE / Nk, SbE / NbEff);
            fprintf(data, "%f\t%f\t%f\n", NbEff, SkE / Nk, SbE / NbEff);
        }
        itmp = (int)(round(ratio * 100.));
        printf("%6d. blok:  %d%% prihvacenih,  Eb = %10.2e\n", NbEff, itmp, SkE / Nk);
    } // kraj petlje blokova

    fclose(data);
    // fclose(koordinate);
    return 0;
}

// probna valna funkcija
double Psi(double r)
{
    double alpha = 4.16;
    double gamma = 2.82;
    double s = 0.0027;
    return exp(-pow(alpha / r, gamma) - s * r) / sqrt(r);
}

double E_pot_L(double r12, double r13, double r23)
{
    return U_LJ(r12) + U_LJ(r13) + U_LJ(r23);
}

double U_LJ(double r)
{
    double sigma = 4;    // angstrema - kako?
    double epsilon = 12; // dubina jame, u kelvinima preko boltzmannove konstante - kako?
    return 4 * epsilon * (pow((sigma / r), 12) - pow((sigma / r), 6));
}

double E_kin_L(double r12, double r13, double r23)
{
    double mass = 4.; // u - kako?
    double hbar = 1.; // koju vrijednost trebam ovdje?
    double D = pow(hbar, 2) / (2 * mass);
    return -D * (f_ddr(r12) + f_ddr(r13) + pow((f_dr(r12) + f_dr(r13)), 2) + f_ddr(r12) + f_ddr(r23) + pow((f_dr(r12) + f_dr(r23)), 2) + f_ddr(r13) + f_ddr(r23) + pow((f_dr(r13) + f_dr(r23)), 2));
}

double f_dr(double r)
{
    double alpha = 4.16;
    double gamma = 2.82;
    double s = 0.0027;
    return 1 / pow(r, 2) * (gamma * pow((alpha / r), gamma) - s - 1 / 2);
}

double f_ddr(double r)
{
    double alpha = 4.16;
    double gamma = 2.82;
    double s = 0.0027;
    return 1 / pow(r, 2) * (gamma * (gamma + 1) * pow((alpha / r), gamma) - s - 1 / 2);
}
