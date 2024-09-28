// Varijacijski Monte Carlo
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"
#include "global_vars.h" // import shared variables

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
#define N_r12_dist 1000                   // broj binova za distribuciju r12
#define N_r12_r13_dist 100                // broj binova za distribuciju r12 i r13 - pazi, 1000 mu bude previše
#define N_angles_dist 1000                // broj binova za distribuciju kuteva
#define L0 30. * A             // angstrema
#define mass 4. * u            // u
// // U1
// #define sigma 4 * A                    // angstrema
// #define epsilon_initial 12 * k_B *K    // dubina jame, u kelvinima preko boltzmannove konstante
// #define alpha_initial 4.55 * A         // angstrema
// #define gamma_initial 4.77             // eksponent u probnoj valnoj funkciji
// #define s_initial 0.3 / A              // eksponent u probnoj valnoj funkciji A^-1
// U2
#define sigma 8 * A                    // angstrema
#define epsilon_initial 20 * k_B *K    // dubina jame, u kelvinima preko boltzmannove konstante
#define alpha_initial 9.2 * A         // angstrema
#define gamma_initial 7.75             // eksponent u probnoj valnoj funkciji
#define s_initial 0.75 / A              // eksponent u probnoj valnoj funkciji A^-1
// // U3
// #define sigma 4 * A            // angstrema
// #define epsilon_initial 33.2 * k_B *K    // dubina jame, u kelvinima preko boltzmannove konstante
// #define alpha_initial 4.82 * A // angstrema
// #define gamma_initial 5.10     // eksponent u probnoj valnoj funkciji
// #define s_initial 0.77 / A      // eksponent u probnoj valnoj funkciji A^-1


// deklaracija funkcija
double Psi(double);                                                                     // probna valna funkcija (korelacijska funkcija)
double U_LJ(double);                                                                    // Lennard-Jonesov potencijal
double E_kin_L(double, double, double, double, double, double, double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double);                                                 // potencijalni dio lokalne energije
double f_ddr(double), f_dr(double);

// OČEKIVANA ENERGIJA BI TREBALA BITI OKO -5 KELVINA
void VMC(double *E_return, double *sigmaE_return, double *r2_return, double *sigmar2_return, int Nt, int Nw, int Nb, int NbSkip)
{
    printf("\nVMC:\n Nt=%d; Nw=%d; Nb=%d; NbSkip=%d\n", Nt, Nw, Nb, NbSkip);
#pragma region // VARIJABLE
    long idum = -1234;
    int i, ib, it, iw, k;
    double x[4][Nw + 1], y[4][Nw + 1];                               // indeks čestice, indeks šetača
    double x_old[4], y_old[4];
    double dx, dy;                                                   // promjene koordinata čestica
    double dxyMax = L0 / 100;                                        // maksimalne promjene x,y,z koordinata
    double random_number;
    double psi_initial, psi_final;                                   // ukupna inicijalna i finalna probna valna funkcija
    double P[Nw + 1];                                                // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
    double E_L[Nw + 1];                                              // lokalna energija zadnjeg koraka svakog od šetača
    double r12, r23, r13;                                            // udaljenosti između čestica
    double x1, x2, x3, y1, y2, y3;                                   // trenutne koordinate svih čestica
    double T;                                                        // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    double ratio;
    double SwE;                                                      // = suma(srednjih E) po setacima
    double StE;                                                      // = suma (srednjih E) po koracima
    double SbE;                                                      // = suma (srednjih E) po blokovima
    double SbE2;                                                     // = suma (srednjih E^2) po blokovima
    double Swr, Str, Sbr, Sbr2;                                      // za srednju kvadratnu vrijednost r
    double Sw_r2, St_r2, Sb_r2, Sb_r2_2;                                      // za srednju vrijednost r
    double AE, sigmaE;                                               // srednja vrijednost i standardna devijacija
    double Ar, sigmar;                                               // srednja vrijednost i standardna devijacija
    double Ar2, sigmar2;                                               // srednja vrijednost i standardna devijacija
    int NbEff;                                                       // efektivni indeks bloka
    int itmp;                                                        // postotak prihvaćanja
    double angle_alpha, angle_beta, angle_gamma;                                                    // kut između r12 i r13 trimera (za svakog šetača posebno)
    // double r12_dist[N_r12_dist + 1], angles_gamma_dist[N_angles_dist + 1], angles_beta_dist[N_angles_dist + 1], angles_alpha_dist[N_angles_dist + 1]; // distribucija duljina r12 i kuteva
    // double r12_r13_dist[N_r12_r13_dist + 1][N_r12_r13_dist + 1];     // distribucija duljina r12 i r13
    double r12_dist[N_r12_dist], angles_gamma_dist[N_angles_dist], angles_beta_dist[N_angles_dist], angles_alpha_dist[N_angles_dist]; // distribucija duljina r12 i kuteva
    double r12_r13_dist[N_r12_r13_dist][N_r12_r13_dist];     // distribucija duljina r12 i r13
    double max_r12 = 20, max_r13 = 20, max_angle = 180;             // max_angle je pi, koje su dobre vrijednosti max_r12 i max_r13?
    int n, m;                                                        // indeksi za distribucije
    double E_kin_calc, E_pot_calc;
    double fdr_value[3], fddr_value[3];
    double x_acos;
    double r2_MSD, r_CM, x_CM, y_CM;
#pragma endregion

    FILE *data, *data_angles, *data_r12, *data_r12_r13, *data_coordinates, *data_log;
    data = fopen("VMC_data.txt", "w");
    data_angles = fopen("VMC_data_angles.txt", "w");
    data_r12 = fopen("VMC_data_r12.txt", "w");
    data_r12_r13 = fopen("VMC_data_r12_r13.txt", "w");
    data_coordinates = fopen("VMC_data_coordinates.txt", "w");
    data_log = fopen("VMC_data_log_VMC.txt", "w");

    // // Initialize 1D arrays
    // for (int i = 0; i < N_r12_dist; i++) {
    //     r12_dist[i] = 0.0;
    // }

    // for (int i = 0; i < N_angles_dist; i++) {
    //     angles_alpha_dist[i] = 0.0;
    //     angles_beta_dist[i] = 0.0;
    //     angles_gamma_dist[i] = 0.0;
    // }

    // // Initialize 2D array
    // for (int i = 0; i < N_r12_r13_dist; i++) {
    //     for (int j = 1; j <= N_r12_r13_dist; j++) {
    //         r12_r13_dist[i][j] = 0.0;
    //     }
    // }

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
    Sbr = 0.;
    Sb_r2 = 0.;
    SbE2 = 0;
    Sbr2 = 0.;
    Sb_r2_2 = 0;
    for (ib = 1; ib <= Nb; ib++) // po blokovima
    {
        Str = 0;
        StE = 0;
        St_r2 = 0;
        NbEff = ib - NbSkip;
        accepted = 0;
        rejected = 0;
        for (it = 1; it <= Nt; it++) // po koracima
        {
            Swr = 0;
            Sw_r2 = 0;
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
                if (T > 1) // prihvaćamo pomak
                {
                    accepted++;
                    P[iw] = psi_final;
                }
                else
                {
                    random_number = ran1(&idum);
                    if (random_number <= T) // prihvaćamo pomak
                    {
                        accepted++;
                        P[iw] = psi_final;
                    }
                    else // odbacujemo pomak
                    {
                        rejected++;
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
                x_CM = (x1 + x2 + x3)/3;
                y_CM = (y1 + y2 + y3)/3;

                E_kin_calc = E_kin_L(r12, r13, r23, x1, x2, x3, y1, y2, y3);
                E_pot_calc = E_pot_L(r12, r13, r23);
                E_L[iw] = E_kin_calc + E_pot_calc; // kinetički dio + potencijalni dio
                SwE = SwE + E_L[iw];
                Swr = Swr + r12;
                r_CM = sqrt(x_CM*x_CM+y_CM*y_CM);
                r2_MSD = 0;
                r2_MSD += (x_CM-x1)*(x_CM-x1) + (y_CM-y1)*(y_CM-y1) + (x_CM-x2)*(x_CM-x2)+(y_CM-y2)*(y_CM-y2) + (x_CM-x3)*(x_CM-x3) +(y_CM-y3)*(y_CM-y3);
                r2_MSD = r2_MSD/3;
                Sw_r2 = Sw_r2 + r2_MSD;

                fdr_value[0] = f_dr(r12); 
                fdr_value[1] = f_dr(r13); 
                fdr_value[2] = f_dr(r23); 
                fddr_value[0] = f_ddr(r12); 
                fddr_value[1] = f_ddr(r13); 
                fddr_value[2] = f_ddr(r23); 

                // data log
                if(ib==150 && it==500){
                    fprintf(data_log, "iw:%d\tE=%f\tE_pot=%f\tE_kin=%f\tr12=%f(fdr=%f; fddr=%f); r13=%f(fdr=%f; fddr=%f); r23=%f(fdr=%f; fddr=%f)\n", iw, E_L[iw], E_pot_calc, E_kin_calc, r12, fdr_value[0], fddr_value[0], r13, fdr_value[1], fddr_value[1], r23, fdr_value[2], fddr_value[2]);
                }

                // // ubacujemo svakog šetača u svakom koraku u distribucije ako je simulacija stabilizirana (NbSkip blokova preskočeno)
                // if (ib > NbSkip)
                // {
                //     n = (int)(r12 / max_r12 * N_r12_dist); 
                //     // printf("nti_coord = %f < r12=%f < n+1ti_coord = %f\n", n*max_r12/N_r12_dist, r12, (n+1)*max_r12/N_r12_dist);
                //     if (n <= N_r12_dist)
                //         r12_dist[n]++;

                //     n = (int)(r12 / max_r12 * N_r12_r13_dist); 
                //     m = (int)(r13 / max_r13 * N_r12_r13_dist);
                //     if (n <= N_r12_r13_dist && m <= N_r12_r13_dist)
                //         r12_r13_dist[n][m]++;

                //     // kut gamma je kut između vektora r12 i r13
                //     x_acos = (r23 * r23 - r12 * r12 - r13 * r13) / (-2 * r12 * r13);
                //     if(x_acos > 1) x_acos = 1;
                //     if(x_acos < -1) x_acos = -1;
                //     angle_gamma = acos(x_acos); 
                //     n = (int)(angle_gamma / 3.14 * N_angles_dist);  
                //     if (n <= N_angles_dist) 
                //         angles_gamma_dist[n]++;

                //     // kut alpha je kut između vektora r13 i r23
                //     x_acos = (r12 * r12 - r23 * r23 - r13 * r13) / (-2 * r23 * r13);
                //     if(x_acos > 1) x_acos = 1;
                //     if(x_acos < -1) x_acos = -1;
                //     angle_alpha = acos(x_acos); 
                //     n = (int)(angle_alpha / 3.14 * N_angles_dist);  
                //     if (n <= N_angles_dist) 
                //         angles_alpha_dist[n]++;

                //     // kut beta je kut između vektora r12 i r23
                //     x_acos = (r13 * r13 - r12 * r12 - r23 * r23) / (-2 * r12 * r23);
                //     if(x_acos > 1) x_acos = 1;
                //     if(x_acos < -1) x_acos = -1;
                //     angle_beta = acos(x_acos); 
                //     n = (int)(angle_beta / 3.14 * N_angles_dist);  
                //     if (n <= N_angles_dist) 
                //         angles_beta_dist[n]++;                    
                // }
            } // kraj petlje šetača
            // akumulacija podataka nakon stabilizacije
            if (ib > NbSkip)
            {
                StE += SwE / Nw;
                Str += Swr / Nw;
                St_r2 += Sw_r2 / Nw;
            }
        } // kraj petlje koraka
        // nakon svakog bloka podeđavamo maksimalnu duljinu koraka kako bi prihvaćanje bilo oko 50%
        ratio = (double)accepted / (double)(accepted + rejected);

        if (ratio > 0.5)
            dxyMax = dxyMax * 1.05;
        if (ratio < 0.5)
            dxyMax = dxyMax * 0.95;

        if (ib > NbSkip)
        {
            SbE += StE / Nt;
            Sbr += Str / Nt;
            Sb_r2 += St_r2 / Nt;
            SbE2 += StE*StE / (Nt*Nt);
            Sbr2 += Str*Str / (Nt*Nt);
            Sb_r2_2 += St_r2*St_r2 / (Nt*Nt);
            fprintf(data, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff, Str / Nt, Sbr / NbEff, St_r2 / Nt, Sb_r2 / NbEff); // indeks bloka, srednji E po koracima (po jednom bloku), Srednji E po blokovima (od početka simulacije), srednja kvadratna vrijednost r12 po koracima (u jednom bloku), -||- od početka simulacije
            itmp = (int)(round(ratio * 100.));
            printf("%6d. blok:  %d%% prihvacenih,  Eb = %10.2e, rb = %f, rb2 = %f\n", NbEff, itmp, StE / Nt, Str / Nt, St_r2 / Nt);
        }
    } // kraj petlje blokova

    // double ukupni_br = Nw*NbEff*Nt;

    // for (i = 0; i < N_r12_dist; i++) // po binovima raspodjele r12
    // {
    //     fprintf(data_r12, "%f\t%f\t\n", (double)(i+0.5) * max_r12 / N_r12_dist, r12_dist[i]/ukupni_br);
    // }

    // for (i = 0; i < N_angles_dist; i++) // po binovima raspodjele kutova
    // {
    //     fprintf(data_angles, "%f\t%f\t%f\t%f\n", (double)(i+0.5) * max_angle / N_angles_dist, angles_gamma_dist[i]/ukupni_br, angles_beta_dist[i]/ukupni_br, angles_alpha_dist[i]/ukupni_br);
    // }

    // for (i = 0; i < N_r12_r13_dist; i++) // po binovima raspodjele r12 i r13
    // {
    //     for (k = 0; k < N_r12_r13_dist; k++)
    //     {
    //         fprintf(data_r12_r13, "%f\t%f\t%f\n", (double)(i+0.5) * max_r12 / 100, (double)(k+0.5) * max_r13 / 100, r12_r13_dist[i][k]/ukupni_br);
    //     }
    // }

    // zapiši konačne koordinate svih šetača u datoteku
    for (iw = 1; iw <= Nw; iw++)
    {
        fprintf(data_coordinates, "%f\t%f\t%f\t%f\t%f\t%f\n", x[1][iw], y[1][iw], x[2][iw], y[2][iw], x[3][iw], y[3][iw]);
    }

    AE = SbE / NbEff;
    Ar = Sbr / NbEff;
    Ar2 = Sb_r2 / NbEff;
    sigmaE = sqrt(fabs(SbE2 / NbEff - AE * AE) / (NbEff - 1.));
    sigmar = sqrt(fabs(Sbr2 / NbEff - Ar * Ar) / (NbEff - 1.));
    sigmar2 = sqrt(fabs(Sb_r2_2 / NbEff - Ar2 * Ar2) / (NbEff - 1.));
    printf(" konacni max. korak: %6.2e\n", dxyMax);
    printf(" alpha = %f, gamma = %f, s = %f, epsilon = %f\n", alpha, gamma_var, s, epsilon);
    printf(" E = %8.5e +- %6.2e \n", AE, sigmaE);
    printf(" r = %8.5e +- %6.2e \n", Ar, sigmar);
    printf(" r2 = %8.5e +- %6.2e \n\n", Ar2, sigmar2);
    *E_return = AE;
    *sigmaE_return = sigmaE;
    *r2_return = Ar2;
    *sigmar2_return = sigmar2;
    fclose(data);
    fclose(data_angles);
    fclose(data_r12);
    fclose(data_r12_r13);
    fclose(data_coordinates);
    fclose(data_log);
}

// probna valna funkcija
double Psi(double r)
{
    return exp(-pow(alpha / r, gamma_var) - s * r) / sqrt(r);
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
    double fdr = 1 / pow(r, 2) * (gamma_var * pow((alpha / r), gamma_var) - s * r - 0.5);
    return fdr;
}

double f_ddr(double r)
{
    double fddr = -1 / pow(r, 2) * (pow(gamma_var, 2) * pow((alpha / r), gamma_var) + s * r);
    return fddr;
}
