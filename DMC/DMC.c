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
#define sigma 4 * A                                      // angstrema
#define epsilon 12 * k_B * K                             // dubina jame, u kelvinima preko boltzmannove konstante
#define L0 30. * A                                       // angstrema
#define alpha_initial 4.55 * A                           // angstrema
#define gamma_initial 4.77                               // eksponent u probnoj valnoj funkciji
#define s_initial 0.3 / A                                // eksponent u probnoj valnoj funkciji A^-1
#define mass 4. * u                                      // u
#define dtau 1.0 * pow(10, -3) / K                       // korak vremena ∆τ (10^(-6) 1/mK)
#define N_r12_dist 100         // broj binova za distribuciju r12
#define N_angles_dist 100      // broj binova za distribuciju kuteva

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
  int nw_max_remove;                           // broj šetača koje smijemo maknuti prije nego broj padne ispod (1-percentage)*Nw0
  int nw_deficit;                              // koliko šetača se treba nadodati u slučaju da je broj pao ispod (1-percentage)*Nw0
  int n_w_temp[Nw_max+1];                      // temporary lista broja potomaka
  double SwE;                                  // = suma(srednjih E) po setacima
  double StE;                                  // = suma (srednjih E) po koracima
  double SbE;                                  // = suma (srednjih E) po blokovima
  double SbE2;                                 // = suma (srednjih E^2) po blokovima
  double AE, sigmaE;                           // srednja vrijednost i standardna devijacija
  int NbEff;                                   // efektivni indeks bloka
  int Nw_temp;                                 // temporary broj setaca
  double rand_num;                             // random broj
  int plus_minus = 0, deltaNw, count;          // broj dodanih/oduzetih pojedinih šetača
  float remainder;                             // ostatak dijeljenja
  double fddr_value[3], fdr_value[3];          // radi logiranja vrijedosti f_dr i f_ddr
  double max_r12 = 100, max_angle = 3.20;      // max_angle je pi, koje su dobre vrijednosti max_r12 i max_r13?
  double r12_dist[N_r12_dist + 1], angles_dist[N_angles_dist + 1]; // distribucija duljina r12 i kuteva
  int n;                                    // indeksi za distribucije
  double angle;                                                    // kut između r12 i r13 trimera (za svakog šetača posebno)
  double x_acos;
#pragma endregion

  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  system(batchScript);

  printf("\nDMC:\n Nt=%d; Nw0=%d; Nb=%d; NbSkip=%d\n", Nt, Nw0, Nb, NbSkip);

  FILE *data, *VMC_coordinates, *data_log, *data_r12, *data_angles, *data_coordinates;
  data = fopen("data.txt", "w");
  VMC_coordinates = fopen("../VMC/data_coordinates.txt", "r");  
  if (VMC_coordinates == NULL) {
    printf("Error opening VMC_coordinates file.\n");
    exit(1);
  }

  data_log = fopen("data_log_DMC.txt", "w");
  data_angles = fopen("data_angles_DMC.txt", "w");
  data_r12 = fopen("data_r12_DMC.txt", "w");
  data_coordinates = fopen("data_coordinates_DMC", "w");


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
    E_kin_calc = E_kin_L(r12, r13, r23, x[1][iw], x[2][iw], x[3][iw], y[1][iw], y[2][iw], y[3][iw]);
    E_pot_calc = E_pot_L(r12, r13, r23);
    E_L[iw] = E_kin_calc + E_pot_calc;
    SwE += E_L[iw];
  }
  E_R = SwE / Nw; // prosječna energija po šetačima = -5.271490

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
      for (iw = 1; iw <= Nw; iw++) // po šetačima
      {
        // korak 1. - gaussov pomak (Ra)
        for (k = 1; k <= 3; k++) // po česticama
        { 
          dx = gasdev(&idum) * stand_dev; // množi sa standardnom devijacijom
          dy = gasdev(&idum) * stand_dev; 
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
        W_Rpw_value = exp(-(0.5 * (E_L[iw] + E_L_prime) - E_R) * dtau);
        W_Rpw[iw-1].value = W_Rpw_value; // E_R ovdje je prosjek koji se mijenja kroz simulaciju pri svakom koraku
        W_Rpw[iw-1].index = iw;

        // korak 9. - stohastička procjena broja potomaka
        rand_num = ran1(&idum);
        n_w[iw] = (int)(W_Rpw_value + rand_num); // trebalo bi uvest optimizaciju koja bi se riješila nekih od ovih šetača        
        sum_nw += n_w[iw] - 1; // suma nadodanih

        fdr_value[0] = f_dr(r12);
        fdr_value[1] = f_dr(r13);
        fdr_value[2] = f_dr(r23);
        fddr_value[0] = f_ddr(r12);
        fddr_value[1] = f_ddr(r13);
        fddr_value[2] = f_ddr(r23);
        if(ib==150 && it==500){
          fprintf(data_log, "iw:%d\tE=%f(p:%f) => w_pr=%f\tE_pot=%f\tE_kin=%f\tr12=%f(fdr=%f; fddr=%f); r13=%f(fdr=%f; fddr=%f); r23=%f(fdr=%f; fddr=%f)\n", iw, E_L_prime, E_L[iw], W_Rpw_value, E_pot_calc, E_kin_calc, r12, fdr_value[0], fddr_value[0], r13, fdr_value[1], fddr_value[1], r23, fdr_value[2], fddr_value[2]);
        }

        E_L[iw] = E_L_prime;
      } // kraj petlje šetača
      
      // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
      deltaNw = Nw_max - Nw;
      if (sum_nw > deltaNw) // što za sum_nw < -Nw_max ?
      {
        // VJEROJATNO TU TREBA NEKA SOFISTICIRANIJA METODA KOJA ĆE UZETI U OBZIR KOJE WALKERE TREBA PRIORITIZIRATI!!
        count = 0; 
        for (iw = 1; iw <= Nw; iw++) 
        { 
          if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
          {
            plus_minus = -1; 
            count += plus_minus;
          }
          if(n_w[iw] > 0){
            plus_minus = (int)((n_w[iw] - 1)*deltaNw / sum_nw); // n_w' / slobodno = n_w / suma
            remainder = (float)(((n_w[iw] - 1)*deltaNw) % sum_nw )/(float)sum_nw;
            if (remainder >= 0.5){
              plus_minus += 1; // round up
            }
              count += plus_minus;
            if (count > deltaNw) 
            {
              plus_minus -= count - deltaNw; // oduzmemo razliku
              count = deltaNw; // spustimo na maskimalnu dozvoljenu popunjenost
            }
          }
          n_w[iw] = 1 + plus_minus; 
        }
      }

      // ako šetači padaju na <70%*Nw0, onda +1 optimalne šetače  
      nw_max_remove = Nw - Nw_min; //razlika između trenutnog - 70*Nw0 (minimalnog broja) - koliko ih se maksimalno smije uništiti
      nw_deficit = abs(sum_nw)-nw_max_remove; // deficit šetača do 70%*Nw0
      if (sum_nw < 0 && abs(sum_nw) > nw_max_remove)
      {
        qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // sortira W_Rpw prema values, descending
        int count = 0;
        while(count<nw_deficit){
          for (j = 0; j<Nw; j++){
            if(n_w[W_Rpw[j].index]!=0){
              n_w[W_Rpw[j].index]++;
              count++;
              if(count>=nw_deficit){
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
      indeks = 1;
      for (iw = 1; iw <= Nw_temp; iw++)
      {
        if (n_w_temp[iw] != 0) // ako je n = 0 onda taj šetač biva uništen (preskačemo ga)
        {
          for (int in = 0; in < n_w_temp[iw]; in++)
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
        if (iw != Nw_temp) // ako nije zadnji korak
        {
          indeks += n_w_temp[iw];
        }
        else if (iw == Nw_temp && n_w_temp[iw] == 0)
        {
          indeks = indeks - 1; // brijem ide -1 jer kao ne trebam se dalje pomaknut od zadnjeg indeksa
        }

        // // ubacujemo svakog šetača u svakom koraku u distribucije ako je simulacija stabilizirana (NbSkip blokova preskočeno)
        // if (ib > NbSkip)
        // {
        //   n = (int)(r12 / max_r12 * 100); // puca jer ih nema tolko unutra u polju - podijelit s nekim reasonable brojem, mislim, i ako je veće od toga valjda se sam zanemari..
        //   if (n <= N_r12_dist)
        //     r12_dist[n]++;          
        //   angle = acos((r23 * r23 - r12 * r12 - r13 * r13) / (-2 * r12 * r13)); // double checkaj ovu formulu
        //   n = (int)(angle / max_angle * 100);                                   // podijelit tipa s pi mislim
        //   if (n <= N_angles_dist)
        //     angles_dist[n]++;
        // }
      }
      Nw = indeks; 

      for(iw = 1; iw<=Nw; iw++){
        SwE += E_L[iw];
      }

      // korak 10. - akumulacija lokalnih energija nakon stabilizacije (za svaki korak)
      if (ib > NbSkip)
      {
        StE += SwE / Nw;
      }
      E_R = SwE / Nw; // ovo je srednja vrijednost energije po svim šetačima, računamo za svaki korak
    } // kraj petlje koraka

    if (ib > NbSkip)
    {
      SbE += StE / Nt;
      SbE2 += StE*StE / (Nt*Nt);
      fprintf(data, "%d\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff); // broj efektivnih blokova, prosjek unutar bloka, prosjek od početka simulacije
      printf("%6d. blok: Nw=%d\tEb = %f\n", NbEff, Nw, StE / Nt);
    }
  } // kraj petlje blokova

  // ubacujemo svakog šetača u svakom koraku u distribucije ako je simulacija stabilizirana (NbSkip blokova preskočeno)
  for (iw = 1; iw < Nw; iw++)
  {
    r12 = sqrt(pow((x[2][iw] - x[1][iw]), 2) + pow((y[2][iw] - y[1][iw]), 2));
    r23 = sqrt(pow((x[3][iw] - x[2][iw]), 2) + pow((y[3][iw] - y[2][iw]), 2));
    r13 = sqrt(pow((x[3][iw] - x[1][iw]), 2) + pow((y[3][iw] - y[1][iw]), 2));
    n = (int)(r12 / max_r12 * 100); // puca jer ih nema tolko unutra u polju - podijelit s nekim reasonable brojem, mislim, i ako je veće od toga valjda se sam zanemari..
    if (n <= N_r12_dist)
      r12_dist[n]++;          

    x_acos = (r23 * r23 - r12 * r12 - r13 * r13) / (-2 * r12 * r13);
    if(x_acos > 1) x_acos = 1;
    if(x_acos < -1) x_acos = -1;
    angle = acos(x_acos); // double checkaj ovu formulu
    n = (int)(angle / max_angle * 100);                                   // podijelit tipa s pi mislim
    if (n <= N_angles_dist)
      angles_dist[n]++;
  }

  // zapiši konačne koordinate svih šetača u datoteku
  for (iw = 1; iw <= Nw; iw++)
  {
    fprintf(data_coordinates, "%f\t%f\t%f\t%f\t%f\t%f\n", x[1][iw], y[1][iw], x[2][iw], y[2][iw], x[3][iw], y[3][iw]);
  }
  for (i = 1; i <= N_r12_dist; i++) // po binovima raspodjele r12
  {
    fprintf(data_r12, "%f\t%f\t\n", (double)i * max_r12 / 100.0, r12_dist[i]);
  }
  for (i = 1; i <= N_angles_dist; i++) // po binovima raspodjele kutova
  {
    fprintf(data_angles, "%f\t%f\t\n", (double)i * max_angle / 100, angles_dist[i]);
  }

  AE = SbE / NbEff;
  sigmaE = sqrt(fabs(SbE2 / NbEff - AE * AE) / (NbEff - 1.));
  printf(" alpha = %f, gamma = %f, s = %f\n", alpha, gamma_var, s);
  printf(" E = %8.5e +- %1.5e \n\n", AE, sigmaE);
  *E_return = AE;
  *sigmaE_return = sigmaE;

  fclose(data);
  fclose(VMC_coordinates);
  fclose(data_log);
  fclose(data_angles);
  fclose(data_r12);
  fclose(data_coordinates);
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

