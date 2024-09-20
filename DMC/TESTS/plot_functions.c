#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_POINTS 1000  // Broj točaka
#define X_MIN 0.1      // Početak intervala
#define X_MAX 101      // Kraj intervala

// konstante
#define k_B 1.                 // boltzmannova konstanta
#define A 1.                   // angstrem
#define K 1.                   // kelvin

// za Lennard-Jones:
#define sigma 4 * A            // angstrema
#define epsilon 12 * k_B *K    // dubina jame, u kelvinima preko boltzmannove konstante
// za korelacijsku f-ju:
#define alpha 4.55 * A // angstrema
#define gamma 4.77     // eksponent u probnoj valnoj funkciji
#define s 0.3 / A      // 1/A - eksponent u probnoj valnoj funkciji

double Psi(double r) // korelacijska funkcija
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

double U_LJ(double r) // lennard-jones potencijal
{
  return 4 * epsilon * (pow((sigma / r), 12) - pow((sigma / r), 6));
}

int main(){
  FILE *output_file;
  output_file = fopen("output.txt", "w"); 

  double step = (X_MAX - X_MIN) / (double)(N_POINTS - 1);  // Izračunaj korak između točaka

  for (int i = 0; i < N_POINTS; i++) {
    double x = X_MIN + i * step;    // Izračunaj trenutnu vrijednost x
    fprintf(output_file, "%lf\t%lf\t%lf\t%lf\t%lf\n", x, U_LJ(x), 1200*Psi(x)*Psi(x), f_dr(x), f_ddr(x));  // Zapiši x i y u datoteku
  }

  fclose(output_file);
  return 0;  
}
