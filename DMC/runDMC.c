#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

// double alpha, gamma_var, s;
// #define dtau 1.0 * pow(10, -3) / K                       // korak vremena ∆τ (10^(-6) 1/mK)

double gamma_var, alpha, s, dtau;

int main()
{
  double d_dtau;
  double dtau_min = 0.1 * pow(10, -3), dtau_max = 1.0 * pow(10, -3); // su kolki??
  double E, sigmaE;

  gamma_var = 4.77;
  alpha = 4.55;
  s = 0.3;
  // int N = 10;
  int N = 5;
  d_dtau = (dtau_max - dtau_min) / N;

  FILE *dataE_dtau;
  dataE_dtau = fopen("dataE_dtau.txt", "w");
  for (int idt = 0; idt < N; idt++)
  {
    dtau = dtau_min + idt * d_dtau;
    printf("dtau=%f\n", dtau);
    // DMC(&E, &sigmaE, 1000, 100, 220, 20);
    DMC(&E, &sigmaE, 1000, 300, 350, 50);
    fprintf(dataE_dtau, "%f\t%f\t%f\t%f\n", E, sigmaE, dtau);
  }
  fclose(dataE_dtau);
}