// Difuzijski Monte Carlo
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"
#include "gasdev.c"

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
#define hbar2_si 1.112121717 * pow(10, −68)            // (Js)^2 = (m2kg/s)^2

// POČETNE VRIJEDNOSTI
#define Nt 100              // broj koraka
#define Nw 100              // broj šetača - al brijem da ovo mora bit fluidno, je li Pero reko oko +/- 30%?
#define Nb 250              // broj blokova
#define NbSkip 0            // broj prvih blokova koje preskačemo
#define sigma 4 * A         // angstrema
#define epsilon 12 * k_B *K // dubina jame, u kelvinima preko boltzmannove konstante
#define L0 30. * A          // angstrema
#define alpha 4.16 * A      // angstrema
#define gamma 2.82          // eksponent u probnoj valnoj funkciji
#define s 0.0027            // eksponent u probnoj valnoj funkciji
#define E_trial 5. * K      // kelvina - ENERGIJA (OČEKIVANO) DOBIVENA IZ VMC ANSAMBLA (MENI DOĐE OKO 2 THO TRENUTNO)
#define mass 4. * u         // u
#define dtau 1.0            // korak vremena ∆τ (čega? u čemu su ove imaginarne jedinice? ni u čemu? veli pero 10 na minus šestu miliKelvina na minus prvu)

// deklaracija funkcija
double Psi(double);                                                                     // probna valna funkcija (korelacijska funkcija)
double U_LJ(double);                                                                    // Lennard-Jonesov potencijal
double E_kin_L(double, double, double, double, double, double, double, double, double); // kinetički dio lokalne energije
double E_pot_L(double, double, double);                                                 // potencijalni dio lokalne energije
double f_ddr(double), f_dr(double);

int main(void)
{
#pragma region // VARIJABLE
  long idum = -1234;
  int ib, it, iw, k, l;                      // indeks bloka, indeks koraka, indeks šetača, indeks čestice
  double x[4][Nw + 1], y[4][Nw + 1];         // indeks čestice, indeks šetača - bar 20% više od referentnog broja šetača
  double x_pww[4][Nw + 1], y_pww[4][Nw + 1]; // indeks čestice, indeks šetača
  double r12, r23, r13;                      // udaljenosti između čestica
  double P[Nw + 1];                          // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
  double dx, dy;                             // promjene koordinata čestica
  double dxyMax = L0 / 100;                  // maksimalne promjene x,y,z koordinata
  double gauss_total[Nw + 1];                // ukupna Gaussova distribucija (za svakog šetača posebno?)
  double psi_total;                          // ukupna probna valna funkcija (za svakog šetača posebno?)
  double sigma2 = hbar2 / mass * dtau;       // varijanca
  // double Fq_x[4][Nw + 1], Fq_y[4][Nw + 1]; // x i y komponente driftne (kvantne) sile
  double Fq_x_prime[4], Fq_y_prime[4]; // x i y komponente driftne (kvantne) sile
  double Fqa_x[4], Fqa_y[4];           // privremena driftna sila međukoraka a u x i y smjeru (Fqa)
  double Fqb_x[4], Fqb_y[4];           // privremena driftna sila međukoraka b u x i y smjeru (Fqb)
  double x_a[4], y_a[4];               // privremena koordinata međukoraka a svake čestice u x i y smjeru (R_a^m)
  double x_b[4], y_b[4];               // privremena koordinata međukoraka b svake čestice u x i y smjeru (R_a^m)
  double x_pw_prime[4], y_pw_prime[4]; // srednji driftni pomak
  double E_L[Nw + 1];                  // lokalna energija zadnjeg koraka svakog od šetača
  double E_L_prime;                    // lokalna energija trenutnog koraka - placeholder za novu energiju
  double E_R;                          // energija R koja se mijenja u svkakom koraku (za svakog šetača??)
  double W_Rpw;                        // statistička težina
  int n_w[Nw];                         // broj potomaka - jel ovo treba bit polje baš? Možda da da znam kasnije za kolko umnožit walkere...
  double SwE;                          // = suma(srednjih E) po setacima
  double StE;                          // = suma (srednjih E) po koracima
  double SbE;                          // = suma (srednjih E) po blokovima
  double SbE2;                         // = suma (srednjih E^2) po blokovima
  int NbEff;                           // efektivni indeks bloka
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
    for (it = 1; it <= Nt; it++) // po koracima
    {
      SwE = 0;
      // u svakom koraku ažuriram broj šetača - vjv prema n_w[iw]
      for (iw = 1; iw <= Nw; iw++) // po šetačima
      {
        // korak 1. - gaussov pomak (Ra)
        for (k = 1; k <= 3; k++) // po česticama
        {
          dx = gasdev(&idum) * sigma2; // množi sa varijancom (varijanca je sigma na kvadrat - double check please da nije korijen iz ovog)
          dy = gasdev(&idum) * sigma2; // množi sa varijancom
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
            if (l != k) // osim l!=k
            {
              // 4.92 VEKTOR
              // spremamo privremeno, x smjer sile za svaku česticu ovog šetača u ovom koraku
              Fqa_x[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (x_a[k] - x_a[l]); // u x-smjeru;
              // spremamo privremeno, y smjer sile za svaku česticu ovog šetača u ovom koraku
              Fqa_y[k] += -2 * f_dr(sqrt(pow((x_a[k] - x_a[l]), 2) + pow((y_a[k] - y_a[l]), 2))) * (y_a[k] - y_a[l]); // u y-smjeru
            }
          }
        }
        // korak 3. pomoćni driftni pomak (Rb)
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
            if (l != k) // osim l!=k
            {
              Fqb_x[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (x_b[k] - x_b[l]); // u x-smjeru
              Fqb_y[k] += -2 * f_dr(sqrt(pow((x_b[k] - x_b[l]), 2) + pow((y_b[k] - y_b[l]), 2))) * (y_b[k] - y_b[l]); // u y-smjeru
            }
          }
        }
        // korak 5. srednji driftni pomak (R'p(w))
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
            if (l != k) // osim l!=k
            {
              Fq_x_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (x_pw_prime[k] - x_pw_prime[l]); // u x-smjeru
              Fq_y_prime[k] += -2 * f_dr(sqrt(pow((x_pw_prime[k] - x_pw_prime[l]), 2) + pow((y_pw_prime[k] - y_pw_prime[l]), 2))) * (y_pw_prime[k] - y_pw_prime[l]); // u y-smjeru
            }
          }
        }
        // i lokalne energije:
        r12 = sqrt(pow((x_pw_prime[2] - x_pw_prime[1]), 2) + pow((y_pw_prime[2] - y_pw_prime[1]), 2));
        r23 = sqrt(pow((x_pw_prime[3] - x_pw_prime[2]), 2) + pow((y_pw_prime[3] - y_pw_prime[2]), 2));
        r13 = sqrt(pow((x_pw_prime[3] - x_pw_prime[1]), 2) + pow((y_pw_prime[3] - y_pw_prime[1]), 2));
        // lokalna energija = kinetički + potencijalni dio
        // PRIME JE SAMO KAO, OZNAKA ZA NOVU VRIJEDNOST U OVOM KORAKU, A KASNIJE SE PREPIŠE U E_L
        E_L_prime = E_kin_L(r12, r13, r23, x_pw_prime[1], x_pw_prime[2], x_pw_prime[3], y_pw_prime[1], y_pw_prime[2], y_pw_prime[3]) + E_pot_L(r12, r13, r23);

        // korak 7. - konačni driftni pomak
        for (k = 1; k <= 3; k++) // po česticama
        {
          x[k][iw] = x_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_x_prime[k];
          y[k][iw] = y_a[k] + (hbar2 / (2 * mass)) * dtau * Fq_y_prime[k];
        }
        // korak 8. - određivanje statističke težine W(R'_p(w))
        // E_L bez crtanog je iz prošlog koraka (stara energija), E_L' je iz trenutnog koraka (nova energija)
        W_Rpw = exp(-(1 / 2 * (E_L[iw] + E_L_prime) - E_R) * dtau); // E_R ovdje još nisam pridala vrijednost!!! mislim da je to neki average...
        // korak 9. - stohastička procjena broja potomaka
        n_w[iw] = (int)(W_Rpw + ran1(&idum));
        // kako se točno kopiraju ovi šetači u nove? Valjda nakon petlje šetača se sve umnoži according to nw[iw]...

        // i brijem da sad spremim energiju... jel se ovo ikad odbacuje, il?
        E_L[iw] = E_L_prime;
        SwE = SwE + E_L[iw];
        // a što je s E_R? to je neki prosjek brijem? prosjek energije čega... svega? svih šetača? ili samo ovog šetača? stavila da je StE/Nt
      } // kraj petlje šetača
      // korak 10. - akumulacija lokalnih energija nakon stabilizacije (za svaki korak)
      if (ib > NbSkip)
      {
        StE += SwE / Nw;
      }
      // mislim da je E_R u biti StE/Nt ?
      E_R = StE / Nt;
      // korak 11. - kopiranje potomaka R'_p(w) u novi ansambl
      // gdje ovo ide?? valjda na kraj ili na početak petlje koraka...
      // TODO za smislit :)
    } // kraj petlje koraka
    if (ib > NbSkip)
    {
      SbE += StE / Nt;
      SbE2 += StE * StE / (Nt * Nt);
      fprintf(data, "%d\t%f\t%f\n", NbEff, StE / Nt, SbE / NbEff); // indeks bloka, srednji E po koracima (po jednom bloku), Srednji E po blokovima (od početka simulacije)
      printf("%6d. blok: Eb = %10.2e\n", NbEff, StE / Nt);
    }
  } // kraj petlje blokova
  fclose(data);
  return 0;
}

// probna valna funkcija
double Psi(double r)
{
  return exp(-pow(alpha / r, gamma) - s * r) / sqrt(r);
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
