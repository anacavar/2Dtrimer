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
#define Nt 1000                                // broj koraka
#define Nw0 100                                // početni broj šetača
#define percentage 0.3                         // +/- varijacije broja šetača
#define Nw_max (int)(1 + percentage) * Nw0 + 1 // početna duljina liste
#define Nb 220                                 // broj blokova
#define NbSkip 20                              // broj prvih blokova koje preskačemo
#define sigma 4 * A                            // angstrema
#define epsilon 12 * k_B *K                    // dubina jame, u kelvinima preko boltzmannove konstante
#define L0 30. * A                             // angstrema
#define alpha_initial 4.55 * A                 // angstrema
#define gamma_initial 4.77                     // eksponent u probnoj valnoj funkciji
#define s_initial 0.3 / A                      // eksponent u probnoj valnoj funkciji A^-1
#define E_trial -5.32 * K                      // kelvina - ENERGIJA (OČEKIVANO) DłOBIVENA IZ VMC-a
#define mass 4. * u                            // u
#define dtau 1.0                               // korak vremena ∆τ (čega? u čemu su ove imaginarne jedinice? ni u čemu? veli pero 10 na minus šestu miliKelvina na minus prvu)

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
  return ((w_pairs *)a)->value - ((w_pairs *)b)->value;
}

void DMC(double *E_return, double *sigmaE_return)
{
#pragma region // VARIJABLE
  long idum = -1234;
  int ib, it, iw, i, j, k, l, indeks;          // indeks bloka, indeks koraka, indeks šetača, indeks čestice
  int Nw = Nw0;                                // broj šetača
  double x[4][Nw_max], y[4][Nw_max];           // indeks čestice, indeks šetača - bar 20% više od referentnog broja šetača
  double x_temp[4][Nw_max], y_temp[4][Nw_max]; // temporary lista
  double r12, r23, r13;                        // udaljenosti između čestica
  double dx, dy;                               // promjene koordinata čestica
  double sigma2 = hbar2 / mass * dtau;         // varijanca
  double Fq_x_prime[4], Fq_y_prime[4];         // x i y komponente driftne (kvantne) sile
  double Fqa_x[4], Fqa_y[4];                   // privremena driftna sila međukoraka a u x i y smjeru (Fqa)
  double Fqb_x[4], Fqb_y[4];                   // privremena driftna sila međukoraka b u x i y smjeru (Fqb)
  double x_a[4], y_a[4];                       // privremena koordinata međukoraka a svake čestice u x i y smjeru (R_a^m)
  double x_b[4], y_b[4];                       // privremena koordinata međukoraka b svake čestice u x i y smjeru (R_a^m)
  double x_pw_prime[4], y_pw_prime[4];         // srednji driftni pomak
  double E_L[Nw_max];                          // lokalna energija zadnjeg koraka svakog od šetača
  double E_L_temp[Nw_max];                     // temporary lista
  double E_L_prime;                            // lokalna energija trenutnog koraka - placeholder za novu energiju
  double E_R;                                  // energija R koja se mijenja u svkakom koraku (za svakog šetača??)
  w_pairs W_Rpw[Nw_max];                       // statistička težina
  int n_w[Nw_max];                             // broj potomaka
  int sum_nw;                                  // suma potomaka nastalih tijekom jednog koraka
  int nw_max_remove;                           // broj šetača koje smijemo maknuti prije nego broj padne ispod (1-percentage)*Nw0
  int n_w_temp[Nw_max];                        // temporary lista broja potomaka
  double SwE;                                  // = suma(srednjih E) po setacima
  double StE;                                  // = suma (srednjih E) po koracima
  double SbE;                                  // = suma (srednjih E) po blokovima
  double SbE2;                                 // = suma (srednjih E^2) po blokovima
  double AE, sigmaE;                           // srednja vrijednost i standardna devijacija
  int NbEff;                                   // efektivni indeks bloka
  int Nw_temp;                                 // temporary broj setaca
#pragma endregion

  FILE *data, *VMC_coordinates;
  data = fopen("data.txt", "w");
  VMC_coordinates = fopen("../VMC/data_coordinates.txt", "r");  
  if (VMC_coordinates == NULL) {
    printf("Error opening VMC_coordinates file.\n");
    exit(1);
  }

  // inicijalizacija liste radi kasnijeg sortiranja
  for (i = 0; i < Nw_max; i++)
  {
    W_Rpw[i].value = 1.;
    W_Rpw[i].index = i;
  }

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
    E_L[iw] = E_kin_L(r12, r13, r23, x[1][iw], x[2][iw], x[3][iw], y[1][iw], y[2][iw], y[3][iw]) + E_pot_L(r12, r13, r23);
    // printf("E_L[%d] = %f\n", iw, E_L[iw]);

    SwE += E_L[iw];
  }
  E_R = SwE / Nw; // prosječna energija po šetačima

  SbE = 0.;
  SbE2 = 0;
  for (ib = 1; ib <= Nb; ib++) // po blokovima
  {
    StE = 0;
    NbEff = ib - NbSkip;
    for (it = 1; it <= Nt; it++) // po koracima
    {
      SwE = 0;
      sum_nw = 0;
      // u svakom koraku ažuriram broj šetača
      // printf("Nw = %d\n", Nw);
      for (iw = 1; iw <= Nw; iw++) // po šetačima
      {
        // korak 1. - gaussov pomak (Ra)
        for (k = 1; k <= 3; k++) // po česticama
        {
          dx = gasdev(&idum) * sigma2; // množi sa varijancom (varijanca je sigma na kvadrat - double check please da nije korijen iz ovog)
          dy = gasdev(&idum) * sigma2; // množi sa varijancom
          // x_a[k] = x[k][iw] + dx;      // spremamo privremeno, x koordinatu za svaku česticu ovog šetača u ovom koraku
          // y_a[k] = y[k][iw] + dy;      // spremamo privremeno, y koordinatu za svaku česticu ovog šetača u ovom koraku
          x[k][iw] = x[k][iw] + dx; // ovo izbriši kasnije!!!!!!!
          y[k][iw] = y[k][iw] + dy; // ovo izbriši kasnije!!!!!!!!!!1
        }
        // // korak 2. - računanje driftne sile (Fa)
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   Fqa_x[k] = 0;
        //   Fqa_y[k] = 0;
        //   for (int l = 1; l <= 3; l++) // po česticama
        //   {
        //     if (l != k) // osim l!=k
        //     {
        //       // 4.92 VEKTOR
        //       // spremamo privremeno, x smjer sile za svaku česticu ovog šetača u ovom koraku
        //       Fqa_x[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (x_a[k] - x_a[l]); // u x-smjeru;
        //       // spremamo privremeno, y smjer sile za svaku česticu ovog šetača u ovom koraku
        //       Fqa_y[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (y_a[k] - y_a[l]); // u y-smjeru
        //     }
        //   }
        // }
        // // korak 3. - pomoćni driftni pomak (Rb)
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   x_b[k] = x_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * Fqa_x[k];
        //   y_b[k] = y_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * Fqa_y[k];
        // }
        // // korak 4. - ponovno računanje driftne sile (Fb)
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   Fqb_x[k] = 0;
        //   Fqb_y[k] = 0;
        //   for (l = 1; l <= 3; l++) // po česticama
        //   {
        //     if (l != k) // osim l!=k
        //     {
        //       Fqb_x[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (x_b[k] - x_b[l]); // u x-smjeru
        //       Fqb_y[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (y_b[k] - y_b[l]); // u y-smjeru
        //     }
        //   }
        // }
        // // korak 5. - srednji driftni pomak (R'p(w))
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   x_pw_prime[k] = x_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * (Fqa_x[k] + Fqb_x[k]) / 2;
        //   y_pw_prime[k] = y_a[k] + (hbar2 / (2 * mass)) * dtau / 2 * (Fqa_y[k] + Fqb_y[k]) / 2;
        //   // pohranimo R_pw u R_b
        //   x_b[k] = x_pw_prime[k];
        //   y_b[k] = y_pw_prime[k];
        // }
        // // korak 6. - računanje driftne sile
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   Fq_x_prime[k] = 0;
        //   Fq_y_prime[k] = 0;
        //   for (l = 1; l <= 3; l++) // po česticama
        //   {
        //     if (l != k) // osim l!=k
        //     {
        //       Fq_x_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (x_pw_prime[k] - x_pw_prime[l]); // u x-smjeru
        //       Fq_y_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (y_pw_prime[k] - y_pw_prime[l]); // u y-smjeru
        //     }
        //   }
        // }
        // i lokalne energije (kinetički + potencijalni dio):
        // r12 = sqrt(pow((x_pw_prime[2] - x_pw_prime[1]), 2) + pow((y_pw_prime[2] - y_pw_prime[1]), 2));
        // r23 = sqrt(pow((x_pw_prime[3] - x_pw_prime[2]), 2) + pow((y_pw_prime[3] - y_pw_prime[2]), 2));
        // r13 = sqrt(pow((x_pw_prime[3] - x_pw_prime[1]), 2) + pow((y_pw_prime[3] - y_pw_prime[1]), 2));
        // E_L_prime = E_kin_L(r12, r13, r23, x_pw_prime[1], x_pw_prime[2], x_pw_prime[3], y_pw_prime[1], y_pw_prime[2], y_pw_prime[3]) + E_pot_L(r12, r13, r23);

        // BRIŠI OVO KAD SVE RADI
        r12 = sqrt(pow((x[2][iw] - x[1][iw]), 2) + pow((y[2][iw] - y[1][iw]), 2));
        r23 = sqrt(pow((x[3][iw] - x[2][iw]), 2) + pow((y[3][iw] - y[2][iw]), 2));
        r13 = sqrt(pow((x[3][iw] - x[1][iw]), 2) + pow((y[3][iw] - y[1][iw]), 2));
        E_L[iw] = E_kin_L(r12, r13, r23, x[1][iw], x[2][iw], x[3][iw], y[1][iw], y[2][iw], y[3][iw]) + E_pot_L(r12, r13, r23);
        // printf("%f\t%f\t%f\tE=%f\n", r12, r23, r13, E_L[iw]);

        // // korak 7. - konačni driftni pomak
        // for (k = 1; k <= 3; k++) // po česticama
        // {
        //   x[k][iw] = x_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_x_prime[k];
        //   y[k][iw] = y_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_y_prime[k];
        // }
        // // korak 8. - određivanje statističke težine W(R'_p(w)) - (E_L bez crtanog je iz prošlog koraka (stara energija), E_L' je iz trenutnog koraka (nova energija))
        // W_Rpw[iw].value = exp(-(1 / 2 * (E_L[iw] + E_L_prime) - E_R) * dtau); // E_R ovdje je neki average koji se mijenja kroz simulaciju...
        // W_Rpw[iw].index = iw;
        // // korak 9. - stohastička procjena broja potomaka
        // n_w[iw] = (int)(W_Rpw[iw].value + ran1(&idum)); // trebalo bi uvest optimizaciju koja bi se riješila nekih od ovih šetača
        // sum_nw += n_w[iw] - 1;
        // printf("W_Rpw[%d] = %f => n_w[%d] = %d => sum_nw = %d, E_L[%d] = %f, E_L_prime = %f, E_R =%f \n", iw, W_Rpw[iw].value, iw, n_w[iw], sum_nw, iw, E_L[iw], E_L_prime, E_R);
        // printf("%f\n", E_L_prime);
        // E_L[iw] = E_L_prime;
        if(E_L[iw] < 5)// izbriši ovaj uvjet!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          SwE = SwE + E_L[iw];
      } // kraj petlje šetača

      // // ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi
      // if (sum_nw > Nw_max - Nw) // || sum_nw < -Nw_max??
      // {
      //   for (iw = 1; iw <= Nw; iw++)
      //   {
      //     // printf("n_w[%d] = %d,\t", iw, n_w[iw]);
      //     // 1 + koliko ih se nadodalo/oduzelo (sum_nw je isto suma nadodanih, a ne ukupnih n[w])
      //     n_w[iw] = 1 + ((int)n_w[iw] - 1) / sum_nw * (Nw_max - Nw); // n_w' / slobodno = n_w / suma
      //     // printf("n_w'[%d] = %d\n", iw, n_w[iw]);
      //   }
      // }

      // // ako ih želi pasti na niže od 30% Nw0, onda ubij samo ove do tih 30%. (pošalji u 0), 
      // // ostale koji bi otišli u 0 pošalji u 1, a ostale ostavi kakvi su bili
      // nw_max_remove = Nw - (int)((1-percentage)*Nw0);
      // if (sum_nw < 0 && abs(sum_nw) > nw_max_remove)
      // {
      //   qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // nisam ziher jel ovo sortira od najvišeg il najnižeg;
      //   // sada je W_Rpw sortiran prema values, dakle npr [(0, 56), (1, 36), (2, 56)] itd...
      //   for (j = 0; j < Nw; j++){
      //     // niš za ove do nw_max_remove
      //     if(j >= nw_max_remove){
      //       if(n_w[W_Rpw[j].index]==0){
      //         n_w[W_Rpw[j].index] = 1;
      //       }
      //     }
      //   }
      // }

      // korak 10. - akumulacija lokalnih energija nakon stabilizacije (za svaki korak)
      if (ib > NbSkip)
      {
        StE += SwE / Nw;
      }
      E_R = StE / Nt; // double check

      // // korak 11. - kopiranje potomaka R'_p(w) u novi ansambl
      // if (Nw > (1 - percentage) * Nw0) // ako broj šetača nije pao ispod 70% početnog broja šetača
      // {
      //   memcpy(x_temp, x, sizeof(x));
      //   memcpy(y_temp, y, sizeof(y));
      //   memcpy(E_L_temp, E_L, sizeof(E_L));
      //   memcpy(n_w_temp, n_w, sizeof(n_w));
      //   Nw_temp = Nw;
      //   Nw = 0;
      //   indeks = 1;
      //   for (iw = 1; iw <= Nw_temp; iw++)
      //   {
      //     if (n_w_temp[iw] != 0) // ako je n = 0 onda taj šetač biva uništen (preskačemo ga)
      //     {
      //       for (int in = 0; in < n_w_temp[iw]; in++)
      //       {
      //         for (k = 1; k <= 3; k++)
      //         {
      //           x[k][indeks + in] = x_temp[k][iw];
      //           y[k][indeks + in] = y_temp[k][iw];
      //         }
      //         E_L[indeks + in] = E_L_temp[iw];
      //         n_w[indeks + in] = n_w_temp[iw];
      //       }
      //     }
      //     if (iw != Nw_temp) // ako nije zadnji korak
      //     {
      //       indeks += n_w_temp[iw];
      //     }
      //     else if (iw == Nw_temp && n_w_temp[iw] == 0)
      //     {
      //       indeks = indeks - 1;
      //     }
      //   }
      //   Nw = indeks; // brijem ide -1 jer kao ne trebam se dalje pomaknut od zadnjeg indeksa
      // }

    } // kraj petlje koraka
    if (ib > NbSkip)
    {
      SbE += StE / Nt;
      SbE2 += StE * StE / (Nt * Nt);
      // jel se dobro računaju ove brojke???
      // broj efektivnih blokova, prosjek unutar bloka, prosjek od početka simulacije
      fprintf(data, "%d\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff); // indeks bloka, srednji E po koracima (po jednom bloku), Srednji E po blokovima (od početka simulacije)
      printf("%6d. blok: Eb = %10.2e\n", NbEff, StE / Nt);
    }
  } // kraj petlje blokova

  AE = SbE / NbEff;

  // valja li ova formula??
  sigmaE = sqrt(abs(SbE2 / NbEff - AE * AE) / (NbEff - 1.));

  printf(" alpha = %f, gamma = %f, s = %f\n", alpha, gamma_var, s);
  printf(" E = %8.5e +- %6.2e \n\n", AE, sigmaE);
  *E_return = AE;
  *sigmaE_return = sigmaE;
  fclose(data);
  fclose(VMC_coordinates);
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
