// Difuzijski Monte Carlo
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ran1.c"
#include "gasdev.c"
#include "global_vars.h" // import shared variables

// KONSTANTE u zadanim mjernim jedinicama i njihov preračun u SI
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

// POČETNE VRIJEDNOSTI
// #define Nt_initial 1000                                  // broj koraka
// #define Nw0_initial 100                                  // početni broj šetača
// #define Nb_initial 220                                   // broj blokova
// #define NbSkip_initial 20                                // broj prvih blokova koje preskačemo
#define percentage 0.3                                   // +/- varijacije broja šetača
// // U1
// #define sigma 4 * A                                      // angstrema
// #define epsilon 12 * k_B * K                             // dubina jame, u kelvinima preko boltzmannove konstante
// #define alpha_initial 4.55 * A                           // angstrema
// #define gamma_initial 4.77                               // eksponent u probnoj valnoj funkciji
// #define s_initial 0.3 / A                                // eksponent u probnoj valnoj funkciji A^-1
// U2
#define sigma 8 * A                                      // angstrema
#define epsilon 20 * k_B * K                             // dubina jame, u kelvinima preko boltzmannove konstante
#define alpha_initial 9.2 * A                           // angstrema
#define gamma_initial 7.75                               // eksponent u probnoj valnoj funkciji
#define s_initial 0.75 / A                                // eksponent u probnoj valnoj funkciji A^-1
// // U3 - dt=0.00001
// #define sigma 4 * A            // angstrema
// #define epsilon 33.2 * k_B *K    // dubina jame, u kelvinima preko boltzmannove konstante
// #define alpha_initial 4.82 * A // angstrema
// #define gamma_initial 5.10     // eksponent u probnoj valnoj funkciji
// #define s_initial 0.77 / A      // eksponent u probnoj valnoj funkciji A^-1

#define mass 4. * u                                      // u
#define N_r12_dist 1000         // broj binova za distribuciju r12
#define N_r12_r13_dist 100         // broj binova za distribuciju r12 i r13
#define N_angles_dist 1000      // broj binova za distribuciju kuteva
// #define dtau 1.0 * pow(10, -3) / K                       // korak vremena ∆τ (10^(-6) 1/mK)

// deklaracija funkcija
double U_LJ(double);                                                                    // Lennard-Jonesov potencijal
double E_kin_L(double, double, double, double, double, double, double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double);                                                 // potencijalni dio lokalne energije
double f_ddr(double), f_dr(double);

typedef struct
{
  double value;
  int index;
} w_pairs;

// Comparator function for sorting w_pairs structs by value
int compareByValue(const void *a, const void *b)
{
  // return ((w_pairs *)a)->value - ((w_pairs *)b)->value; // ascending - 0, 1, 2, 3, ...
  return ((w_pairs *)b)->value - ((w_pairs *)a)->value; // descending - ..., 3, 2, 1, 0
}

void DMC(double *E_return, double *sigmaE_return, int Nt, int Nw0, int Nb, int NbSkip)
{
  printf("primljeni dtau: %f\n", dtau);
#pragma region // VARIJABLE
  long idum = -1234;
  int ib, it, iw, i, j, k, l, indeks;          // indeks bloka, indeks koraka, indeks šetača, indeks čestice
  int Nw = Nw0;                                // broj šetača
  int Nw_max = (int)((1 + percentage) * Nw0);    // maksimalna duljina liste
  int Nw_min = (int)((1-percentage)*Nw0);      // minimalna duljina liste
  double x[4][Nw_max+1], y[4][Nw_max+1];       // indeks čestice, indeks šetača - bar 20% više od referentnog broja šetača
  double x_temp[4][Nw_max+1], y_temp[4][Nw_max+1]; // temporary lista
  double r12, r23, r13;                        // udaljenosti između čestica
  double dx, dy;                               // promjene koordinata čestica
  double sigma2 = hbar2 / mass * dtau;         // varijanca ???
  double stand_dev = sqrt(sigma2);             // standardna devijacija
  double Fq_x_prime[4], Fq_y_prime[4];         // x i y komponente driftne (kvantne) sile
  double Fqa_x[4], Fqa_y[4];                   // privremena driftna sila međukoraka a u x i y smjeru (Fqa) za svaku od čestica
  double Fqb_x[4], Fqb_y[4];                   // privremena driftna sila međukoraka b u x i y smjeru (Fqb) za svaku od čestica
  double x_a[4], y_a[4];                       // privremena koordinata međukoraka a svake čestice u x i y smjeru (R_a^m)
  double x_b[4], y_b[4];                       // privremena koordinata međukoraka b svake čestice u x i y smjeru (R_a^m)
  double x_pw_prime[4], y_pw_prime[4];         // srednji driftni pomak
  double E_L[Nw_max+1];                        // lokalna energija zadnjeg koraka svakog od šetača
  double E_L_temp[Nw_max+1];                   // temporary lista
  double E_L_prime;                            // lokalna energija trenutnog koraka - placeholder za novu energiju
  double E_R;                                  // energija R koja se mijenja u svkakom koraku (za svakog šetača??)
  double E_kin_calc, E_pot_calc;               // kinetički i potencijalni dio energije
  w_pairs W_Rpw[Nw_max];                       // statistička težina
  double W_Rpw_value;                          // statistička težina
  int n_w[Nw_max+1];                           // broj potomaka
  int sum_nw;                                  // suma potomaka nastalih tijekom jednog koraka
  int dodavanje, oduzimanje;                                  // suma potomaka nastalih tijekom jednog koraka
  int nw_max_remove;                           // broj šetača koje smijemo maknuti prije nego broj padne ispod (1-percentage)*Nw0
  int nw_deficit;                              // koliko šetača se treba nadodati u slučaju da je broj pao ispod (1-percentage)*Nw0
  int n_w_temp[Nw_max+1];                      // temporary lista broja potomaka
  double SwE;                                  // = suma(srednjih E) po setacima
  double StE;                                  // = suma (srednjih E) po koracima
  double SbE;                                  // = suma (srednjih E) po blokovima
  double SbE2;                                 // = suma (srednjih E^2) po blokovima
  double Swr, Str, Sbr, Sbr2;                  // za srednju kvadratnu vrijednost r
  double Sw_r2, St_r2, Sb_r2, Sb_r2_2;         // za srednju kvadratnu vrijednost r
  double Ar, sigmar;
  double Ar2, sigmar2;
  double AE, sigmaE;                           // srednja vrijednost i standardna devijacija
  int NbEff;                                   // efektivni indeks bloka
  int Nw_temp;                                 // temporary broj setaca
  double rand_num;                             // random broj
  int plus_minus = 0, deltaNw, count;          // broj dodanih/oduzetih pojedinih šetača
  float remainder;                             // ostatak dijeljenja
  double fddr_value[3], fdr_value[3];          // radi logiranja vrijedosti f_dr i f_ddr
  double max_r12 = 20, max_r13=20, max_angle = 180;      // maksimalne vrijednosti za distribucije
  double r12_dist[N_r12_dist], angles_alpha_dist[N_angles_dist], angles_beta_dist[N_angles_dist], angles_gamma_dist[N_angles_dist]; // distribucija duljina r12 i kuteva
  double r12_r13_dist[N_r12_r13_dist][N_r12_r13_dist];     // distribucija duljina r12 i r13
  int n, m;                                    // indeksi za distribucije
  double angle_alpha, angle_beta, angle_gamma; // kutovi trimera
  double x_acos;
  double x1, x2, x3, y1, y2, y3, x_CM, y_CM, r2_MSD;  // koordinate centra mase trimera
#pragma endregion

  printf("\nDMC:\n Nt=%d; Nw0=%d; Nb=%d; NbSkip=%d\n", Nt, Nw0, Nb, NbSkip);

  FILE *data, *VMC_coordinates, *data_log, *data_r12, *data_angles, *data_coordinates, *data_r12_r13;
  data = fopen("DMC_data.txt", "w");
  // VMC_coordinates = fopen("../VMC/VMC_data_coordinates.txt", "r"); // za lokalno pokretanje
  VMC_coordinates = fopen("VMC_data_coordinates.txt", "r");  // za pokretanje na udaljenom serveru
  if (VMC_coordinates == NULL) {
    printf("Error opening VMC_coordinates file.\n");
    exit(1);
  }

  data_log = fopen("data_log_DMC.txt", "w");
  data_angles = fopen("data_angles_DMC.txt", "w");
  data_r12 = fopen("data_r12_DMC.txt", "w");
  data_r12_r13 = fopen("data_r12_r13_DMC.txt", "w");

  // inicijalizacija liste radi kasnijeg sortiranja
  for (i = 0; i < Nw_max; i++)
  {
    W_Rpw[i].value = 1.;
    W_Rpw[i].index = i;
  }

  // distirbucije pucaju na udaljenom serveru
  // // Initialize 1D arrays
  // for (int i = 0; i < N_r12_dist; i++) {
  //   r12_dist[i] = 0.0;
  // }

  // for (int i = 0; i < N_angles_dist; i++) {
  //   angles_alpha_dist[i] = 0.0;
  //   angles_beta_dist[i] = 0.0;
  //   angles_gamma_dist[i] = 0.0;
  // }

  // // Initialize 2D array
  // for (int i = 0; i < N_r12_r13_dist; i++) {
  //   for (int j = 1; j <= N_r12_r13_dist; j++) {
  //     r12_r13_dist[i][j] = 0.0;
  //   }
  // }

  SwE = 0;
  // inicijalizacija koordinata čestica 
  for (iw = 1; iw <= Nw; iw++) // po šetačima
  {
    // load coordinates from VMC
    fscanf(VMC_coordinates, "%lf %lf %lf %lf %lf %lf", &x[1][iw], &y[1][iw], &x[2][iw], &y[2][iw], &x[3][iw], &y[3][iw]);
    // račun energije za pojedinog šetača(kinetički + potencijalni dio):
    r12 = sqrt(pow((x[2][iw] - x[1][iw]), 2) + pow((y[2][iw] - y[1][iw]), 2));
    r23 = sqrt(pow((x[3][iw] - x[2][iw]), 2) + pow((y[3][iw] - y[2][iw]), 2));
    r13 = sqrt(pow((x[3][iw] - x[1][iw]), 2) + pow((y[3][iw] - y[1][iw]), 2));
    E_kin_calc = E_kin_L(r12, r13, r23, x[1][iw], x[2][iw], x[3][iw], y[1][iw], y[2][iw], y[3][iw]);
    E_pot_calc = E_pot_L(r12, r13, r23);
    E_L[iw] = E_kin_calc + E_pot_calc;
    SwE += E_L[iw];
  }
  E_R = SwE / Nw; // početna prosječna energija po šetačima = -5.271490

  SbE = 0.;
  SbE2 = 0;
  Sbr = 0.;
  Sbr2 = 0;
  Sb_r2 = 0.;
  Sb_r2_2 = 0.;
  for (ib = 1; ib <= Nb; ib++) // po blokovima
  {
    StE = 0;
    Str = 0;
    St_r2 = 0;
    NbEff = ib - NbSkip;
    for (it = 1; it <= Nt; it++) // po koracima
    {
      SwE = 0;
      Swr = 0;
      Sw_r2 = 0;
      dodavanje = 0;
      deltaNw = Nw_max - Nw;
      oduzimanje = 0;
      for (iw = 1; iw <= Nw; iw++) // po šetačima
      {
        // korak 1. - gaussov pomak (Ra)
        for (k = 1; k <= 3; k++) // po česticama
        { 
          dx = gasdev(&idum) * stand_dev; // množi sa standardnom devijacijom
          dy = gasdev(&idum) * stand_dev; 
          // dx = 0; // kad se makne gausov korak
          // dy = 0; // kad se makne gausov korak
          x_a[k] = x[k][iw] + dx;      // spremamo privremeno, x koordinatu za svaku česticu ovog šetača u ovom koraku
          y_a[k] = y[k][iw] + dy;      // spremamo privremeno, y koordinatu za svaku česticu ovog šetača u ovom koraku
        }
        // korak 2. - računanje driftne sile (Fa)
        for (k = 1; k <= 3; k++) // po česticama
        {
          Fqa_x[k] = 0;
          Fqa_y[k] = 0;
          for (int l = 1; l <= 3; l++) // po česticama
          {
            if (l != k) // svi osim l=k
            {
              // 4.92 VEKTOR - privremeno spremamo x i y smjer svake čestice ovog šetača u ovom koraku
              Fqa_x[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (x_a[l] - x_a[k]); 
              Fqa_y[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (y_a[l] - y_a[k]); 
            }
          }
        }
        // korak 3. - pomoćni driftni pomak (Rb)
        for (k = 1; k <= 3; k++) // po česticama
        {
          x_b[k] = x_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * Fqa_x[k];
          y_b[k] = y_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * Fqa_y[k];
        }
        // korak 4. - ponovno računanje driftne sile (Fb)
        for (k = 1; k <= 3; k++) // po česticama
        {
          Fqb_x[k] = 0;
          Fqb_y[k] = 0;
          for (l = 1; l <= 3; l++) // po česticama
          {
            if (l != k) // svi osim l=k
            {
              Fqb_x[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (x_b[l] - x_b[k]); // u x-smjeru
              Fqb_y[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (y_b[l] - y_b[k]); // u y-smjeru
            }
          }
        }
        // korak 5. - srednji driftni pomak (R'_p(w))
        for (k = 1; k <= 3; k++) // po česticama
        {
          x_pw_prime[k] = x_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * (Fqa_x[k] + Fqb_x[k]) / 2;
          y_pw_prime[k] = y_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * (Fqa_y[k] + Fqb_y[k]) / 2;
          // pohranimo R_pw u R_b
          x_b[k] = x_pw_prime[k];
          y_b[k] = y_pw_prime[k];
        }
        // korak 6. - računanje driftne sile
        for (k = 1; k <= 3; k++) // po česticama
        {
          Fq_x_prime[k] = 0;
          Fq_y_prime[k] = 0;
          for (l = 1; l <= 3; l++) // po česticama
          {
            if (l != k) // svi osim l=k
            {
              Fq_x_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (x_pw_prime[l] - x_pw_prime[k]); // u x-smjeru
              Fq_y_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (y_pw_prime[l] - y_pw_prime[k]); // u y-smjeru
            }
          }
        }

        // korak 7. - konačni driftni pomak
        for (k = 1; k <= 3; k++) // po česticama
        {
          x[k][iw] = x_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_x_prime[k];
          y[k][iw] = y_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_y_prime[k];
        }

        // račun lokalne energije (kinetički + potencijalni dio):
        r12 = sqrt(pow((x_pw_prime[2] - x_pw_prime[1]), 2) + pow((y_pw_prime[2] - y_pw_prime[1]), 2));
        r23 = sqrt(pow((x_pw_prime[3] - x_pw_prime[2]), 2) + pow((y_pw_prime[3] - y_pw_prime[2]), 2));
        r13 = sqrt(pow((x_pw_prime[3] - x_pw_prime[1]), 2) + pow((y_pw_prime[3] - y_pw_prime[1]), 2));

        E_kin_calc = E_kin_L(r12, r13, r23, x_pw_prime[1], x_pw_prime[2], x_pw_prime[3], y_pw_prime[1], y_pw_prime[2], y_pw_prime[3]);
        E_pot_calc = E_pot_L(r12, r13, r23);
        E_L_prime = E_kin_calc + E_pot_calc;

        // korak 8. - određivanje statističke težine W(R'_p(w)) - (E_L bez crtanog je iz prošlog koraka (stara energija), E_L' je iz trenutnog koraka (nova energija))
        W_Rpw_value = exp(-(0.5 * (E_L[iw] + E_L_prime) - E_R) * dtau); // E_R ovdje je prosjek koji se mijenja kroz simulaciju pri svakom koraku

        if(W_Rpw_value>1000) {
          // printf("W_Rpw_value overflow! (extreme number)\n");
          W_Rpw_value = 1000;
        }
        if (isnan(W_Rpw_value)) {
          // Handle NaN
          // printf("Value is NaN (likely numerical overflow)\n");
          W_Rpw_value = 1000;
        }

        W_Rpw[iw-1].value = W_Rpw_value; 
        W_Rpw[iw-1].index = iw;

        // korak 9. - stohastička procjena broja potomaka
        rand_num = ran1(&idum);

        // if(n_w[iw]<1){
        //   rand_num = 0.7*rand_num;
        // }
        
        n_w[iw] = (int)(W_Rpw_value + rand_num); // trebalo bi uvest optimizaciju koja bi se riješila nekih od ovih šetača        

        if(n_w[iw]==0){
          oduzimanje++;
        }
        else{
          dodavanje += n_w[iw] - 1;
        }

        fdr_value[0] = f_dr(r12);
        fdr_value[1] = f_dr(r13);
        fdr_value[2] = f_dr(r23);
        fddr_value[0] = f_ddr(r12);
        fddr_value[1] = f_ddr(r13);
        fddr_value[2] = f_ddr(r23);
        if(ib==150 && it==500){
          fprintf(data_log, "iw:%d\tE=%f(p:%f) => w_pr=%f\trand=%f\tE_pot=%f\tE_kin=%f\tr12=%f(fdr=%f; fddr=%f); r13=%f(fdr=%f; fddr=%f); r23=%f(fdr=%f; fddr=%f)\n", iw, E_L_prime, E_L[iw], W_Rpw_value, rand_num, E_pot_calc, E_kin_calc, r12, fdr_value[0], fddr_value[0], r13, fdr_value[1], fddr_value[1], r23, fdr_value[2], fddr_value[2]);
        }

        E_L[iw] = E_L_prime;

        // printf("hello ib: %d, it:%d, iw:%d\n", ib, it, iw);


        // // distribucija radi probleme na udaljenom serveru
        // // ubacujemo svakog šetača u svakom koraku u distribucije ako je simulacija stabilizirana (NbSkip blokova preskočeno)
        // if (ib > NbSkip && n_w[iw]!=0)
        // {
        //   r12 = sqrt(pow((x[1][iw] - x[2][iw]), 2) + pow((y[1][iw] - y[2][iw]), 2));
        //   r13 = sqrt(pow((x[1][iw] - x[3][iw]), 2) + pow((y[1][iw] - y[3][iw]), 2));

        //   n = (int)(r12 / max_r12 * N_r12_dist); 
        //   // printf("nti_coord = %f < r12=%f < n+1ti_coord = %f\n", n*max_r12/N_r12_dist, r12, (n+1)*max_r12/N_r12_dist);
        //   if (n <= N_r12_dist)
        //     r12_dist[n]++;

        //   n = (int)(r12 / max_r12 * N_r12_r13_dist); 
        //   m = (int)(r13 / max_r13 * N_r12_r13_dist);
        //   if (n <= N_r12_r13_dist && m <= N_r12_r13_dist)
        //     r12_r13_dist[n][m]++;

        //   // kut gamma je kut između vektora r12 i r13
        //   x_acos = (r23 * r23 - r12 * r12 - r13 * r13) / (-2 * r12 * r13);
        //   if(x_acos > 1) x_acos = 1;
        //   if(x_acos < -1) x_acos = -1;
        //   angle_gamma = acos(x_acos); 
        //   n = (int)(angle_gamma / 3.14 * N_angles_dist);  
        //   if (n <= N_angles_dist) 
        //     angles_gamma_dist[n]++;

        //   // kut alpha je kut između vektora r13 i r23
        //   x_acos = (r12 * r12 - r23 * r23 - r13 * r13) / (-2 * r23 * r13);
        //   if(x_acos > 1) x_acos = 1;
        //   if(x_acos < -1) x_acos = -1;
        //   angle_alpha = acos(x_acos); 
        //   n = (int)(angle_alpha / 3.14 * N_angles_dist);  
        //   if (n <= N_angles_dist) 
        //     angles_alpha_dist[n]++;

        //   // kut beta je kut između vektora r12 i r23
        //   x_acos = (r13 * r13 - r12 * r12 - r23 * r23) / (-2 * r12 * r23);
        //   if(x_acos > 1) x_acos = 1;
        //   if(x_acos < -1) x_acos = -1;
        //   angle_beta = acos(x_acos); 
        //   n = (int)(angle_beta / 3.14 * N_angles_dist);  
        //   if (n <= N_angles_dist) 
        //     angles_beta_dist[n]++;                    
        // }
      } // kraj petlje šetača
      sum_nw = dodavanje - oduzimanje;
      
      // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
      if (dodavanje > deltaNw+oduzimanje) 
      {
        count = 0; 
        for (iw = 1; iw <= Nw; iw++) 
        { 
          if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
          {
            plus_minus = -1; 
          }
          if(n_w[iw] > 0){
            plus_minus = (int)((n_w[iw] - 1)*(deltaNw+oduzimanje) / dodavanje); // n_w' / slobodno = n_w / suma

            remainder = (float)(((n_w[iw] - 1)*(deltaNw+oduzimanje)) % dodavanje )/(float)dodavanje;
            if (remainder >= 0.5){
              plus_minus += 1; // round up
            }
            count += plus_minus; 

            if (count > deltaNw+oduzimanje) 
            {
              plus_minus -= count - (deltaNw+oduzimanje); // oduzmemo razliku
              count = deltaNw+oduzimanje; // spustimo na maskimalnu dozvoljenu popunjenost
            }
          }
          //dodat još jednu iteraciju dok se count sigurno ne popuni...
          n_w[iw] = 1 + plus_minus; 
        }
        // sljedeće dvije petlje su dopunjavanje do maksimalne popunjenosti liste, prvo prema prioritetu, a onda po redu
        int first_iteration = 1;
        while(count < deltaNw+oduzimanje) // dopunjavanje prema prioritetu
        {
          if(first_iteration == 1){
            for(iw = 1; iw <= Nw; iw++){
              if(n_w[iw]>1){
                n_w[iw]++;
                count++;
              }
              if (count>=deltaNw+oduzimanje){
                break;
              }
            }
            first_iteration=0;
          }
          else if(first_iteration==0 && count<deltaNw+oduzimanje) // dopunjavanje po redu
          {
            for(iw = 1; iw <= Nw; iw++){
              if(n_w[iw]>0){
                n_w[iw]++;
                count++;
              }
              if (count>=deltaNw+oduzimanje){
                break;
              }
            }
          }
        }
      }

      // ako šetači padaju na <70%*Nw0, onda +1 optimalne šetače  
      nw_max_remove = Nw - Nw_min; //razlika između trenutnog - 70*Nw0 (minimalnog broja) - koliko ih se maksimalno smije uništiti
      nw_deficit = abs(sum_nw)-nw_max_remove; // deficit šetača do 70%*Nw0
      if (sum_nw < 0 && abs(sum_nw) > nw_max_remove)
      {
        // printf("hello <70\n");

        qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // sortira W_Rpw prema values, descending
        int count = 0;
        int count_if_all_zero =0;
        while(count<nw_deficit && count_if_all_zero<nw_deficit){
          count_if_all_zero = 0;
          for (j = 0; j<Nw; j++){
            if(n_w[W_Rpw[j].index]!=0){
              n_w[W_Rpw[j].index]++;
              count++;
              if(count>=nw_deficit){
                goto exitLoops;
              }
            }
            if(n_w[W_Rpw[j].index]==0){
              count_if_all_zero++;
              if(count_if_all_zero>=Nw){
                for(int znj = 1; znj<=Nw; znj++){
                  n_w[znj] = 1;
                }
                printf("All values were zero! ib:%d; it:%d, iw%d\n", ib, it, iw);
                goto exitLoops;
              }
            }
          }
        }
        exitLoops:;
      }

      // korak 11. - kopiranje potomaka R'_p(w) u novi ansambl
      memcpy(x_temp, x, sizeof(x));
      memcpy(y_temp, y, sizeof(y));
      memcpy(E_L_temp, E_L, sizeof(E_L));
      memcpy(n_w_temp, n_w, sizeof(n_w));
      Nw_temp = Nw;
      Nw = 0;
      indeks = 0;
      for (iw = 1; iw <= Nw_temp; iw++)
      {

        if (n_w_temp[iw] != 0) // ako je n = 0 onda taj šetač biva uništen (preskačemo ga)
        {
          for (int in = 1; in <= n_w_temp[iw]; in++)
          {
            for (k = 1; k <= 3; k++)
            {
              x[k][indeks + in] = x_temp[k][iw];
              y[k][indeks + in] = y_temp[k][iw];
            }
            E_L[indeks + in] = E_L_temp[iw];
            n_w[indeks + in] = n_w_temp[iw];
          } 
        }
        indeks += n_w_temp[iw]; 

      }
      Nw = indeks; 




      // skupljanje podataka po šetačima nakon redistribucije, prije kraja koraka
      for(iw = 1; iw<=Nw; iw++){
        r12 = sqrt(pow((x[1][iw] - x[2][iw]), 2) + pow((y[1][iw] - y[2][iw]), 2));
        SwE += E_L[iw];
        Swr += r12;
        x1 = x[1][iw];
        x2 = x[2][iw];
        x3 = x[3][iw];
        y1 = y[1][iw];
        y2 = y[2][iw];
        y3 = y[3][iw];
        x_CM = (x1 + x2 + x3)/3;
        y_CM = (y1 + y2 + y3)/3;
        r2_MSD = 0;
        r2_MSD += (x_CM-x1)*(x_CM-x1) + (y_CM-y1)*(y_CM-y1) + (x_CM-x2)*(x_CM-x2)+(y_CM-y2)*(y_CM-y2) + (x_CM-x3)*(x_CM-x3) +(y_CM-y3)*(y_CM-y3);
        r2_MSD = r2_MSD/3;
        Sw_r2 += r2_MSD;
      }

      // korak 10. - akumulacija lokalnih energija nakon stabilizacije (za svaki korak)
      if (ib > NbSkip)
      {
        StE += SwE / Nw - 0.2;
        Str += Swr / Nw;
        St_r2 += Sw_r2 / Nw;
      }
      E_R = SwE / Nw; // ovo je srednja vrijednost energije po svim šetačima, računamo za svaki korak
    } // kraj petlje koraka

    if (ib > NbSkip)
    {
      SbE += StE / Nt;
      SbE2 += StE*StE / (Nt*Nt);
      Sbr += Str / Nt;
      Sbr2 += Str*Str / (Nt*Nt);
      Sb_r2 += St_r2 / Nt;
      Sb_r2_2 += St_r2 * St_r2 / (Nt*Nt);
      fprintf(data, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff, Str / Nt, Sbr / NbEff, St_r2 / Nt, Sb_r2 / NbEff); // broj efektivnih blokova, prosjek unutar bloka, prosjek od početka simulacije
      printf("%6d. blok: Nw=%d\tEb = %f\tr_b = %f\tr2_b = %f\n", NbEff, Nw, StE / Nt, Str / Nt, St_r2 / Nt);
    }
  } // kraj petlje blokova


  // // distribucije rade probleme na udaljenom serveru
  // double ukupni_br = Nw0*NbEff*Nt;

  // for (i = 0; i < N_r12_dist; i++) // po binovima raspodjele r12
  // {
  //   fprintf(data_r12, "%f\t%f\t\n", (double)(i+0.5) * max_r12 / N_r12_dist, r12_dist[i]/ukupni_br);
  // }
  // for (i = 0; i < N_angles_dist; i++) // po binovima raspodjele kutova
  // {
  //   fprintf(data_angles, "%f\t%f\t%f\t%f\n", (double)(i+0.5) * max_angle / N_angles_dist, angles_gamma_dist[i]/ukupni_br, angles_beta_dist[i]/ukupni_br, angles_alpha_dist[i]/ukupni_br);
  // }

  // for (i = 0; i < N_r12_r13_dist; i++) // po binovima raspodjele r12 i r13
  // {
  //     for (k = 0; k < N_r12_r13_dist; k++)
  //     {
  //         fprintf(data_r12_r13, "%f\t%f\t%f\n", (double)(i+0.5) * max_r12 / 100, (double)(k+0.5) * max_r13 / 100, r12_r13_dist[i][k]/ukupni_br);
  //     }
  // }

  AE = SbE / NbEff;
  Ar = Sbr / NbEff;
  Ar2 = Sb_r2 / NbEff;
  sigmaE = sqrt(fabs(SbE2 / NbEff - AE * AE) / (NbEff - 1.));
  sigmar = sqrt(fabs(Sbr2 / NbEff - Ar * Ar) / (NbEff - 1.));
  sigmar2 = sqrt(fabs(Sb_r2_2 / NbEff - Ar2 * Ar2) / (NbEff - 1.));
  printf(" alpha = %f, gamma = %f, s = %f\n", alpha, gamma_var, s);
  printf(" E = %8.5e +- %1.5e \n", AE, sigmaE);
  printf(" r = %8.5e +- %1.5e \n", Ar, sigmar);
  printf(" r2 = %8.5e +- %1.5e \n\n", Ar2, sigmar2);
  *E_return = AE;
  *sigmaE_return = sigmaE;

  fclose(data);
  fclose(VMC_coordinates);
  fclose(data_log);
  fclose(data_angles);
  fclose(data_r12);
  fclose(data_r12_r13);
}

double f_dr(double r)
{
  double fdr = 1 / pow(r, 2) * (gamma_var * pow((alpha / r), gamma_var) - s * r - 1 / 2);
  return fdr;
}

double f_ddr(double r)
{
  double fddr = -1 / pow(r, 2) * (pow(gamma_var, 2) * pow((alpha / r), gamma_var) + s * r);
  return fddr;
}

double U_LJ(double r)
{
  return 4 * epsilon * (pow((sigma / r), 12) - pow((sigma / r), 6));
}

double E_pot_L(double r12, double r13, double r23)
{
  return U_LJ(r12) + U_LJ(r13) + U_LJ(r23);
}

double E_kin_L(double r12, double r13, double r23, double x1, double x2, double x3, double y1, double y2, double y3)
{
  double D = hbar2 / (2 * mass);
  double dio1 = f_ddr(r12) + f_ddr(r13) + pow(f_dr(r12) * (x2 - x1) + f_dr(r13) * (x3 - x1), 2) + pow(f_dr(r12) * (y2 - y1) + f_dr(r13) * (y3 - y1), 2);
  double dio2 = f_ddr(r12) + f_ddr(r23) + pow(f_dr(r12) * (x1 - x2) + f_dr(r23) * (x3 - x2), 2) + pow(f_dr(r12) * (y1 - y2) + f_dr(r23) * (y3 - y2), 2);
  double dio3 = f_ddr(r13) + f_ddr(r23) + pow(f_dr(r13) * (x1 - x3) + f_dr(r23) * (x2 - x3), 2) + pow(f_dr(r13) * (y1 - y3) + f_dr(r23) * (y2 - y3), 2);
  return -D * (dio1 + dio2 + dio3);
}
