#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

// double alpha, gamma_var, s;
// #define dtau 1.0 * pow(10, -3) / K                       // korak vremena ∆τ (10^(-6) 1/mK)

double gamma_var, alpha, s, dtau;

int main()
{
  // int Nt=1000, Nw=300, Nb=250, NbSkip=50;
  int Nt=500, Nw=100, Nb=220, NbSkip=20;

  // vjerojatno ovo treba prebaciti u singlerun i run datoteke.
  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  // system(batchScript);

  double d_dtau;
  // double dtau_min = 0.1 * pow(10, -5), dtau_max = 1.0 * pow(10, -2); // su kolki??
  // double dtau_min = 0.1 * pow(10, -3), dtau_max = 1.0 * pow(10, -3); // su kolki??
  double E, sigmaE;
  FILE *dataE_dtau;
  dataE_dtau = fopen("dataE_dtau.txt", "w");

  gamma_var = 4.77;
  alpha = 4.55;
  s = 0.3;
  int N = 10;
  double dtau_mins[5]={0.000001, 0.00001, 0.0001, 0.001, 0.01};
  double dtau_maxes[5]={0.000011, 0.00011, 0.0011, 0.011, 0.11};

  for (int i = 0; i < 5; i++){
    d_dtau = (dtau_maxes[i] - dtau_mins[i]) / N;
    printf("%d. d_dtau=%f\n", i, d_dtau);
    for (int idt = 0; idt < N; idt++)
    {
      dtau = dtau_mins[i] + idt * d_dtau;
      printf("dtau=%f\n", dtau);
      DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
      fprintf(dataE_dtau, "%f\t%f\t%f\n", dtau, E, sigmaE);
    }
  }

  fclose(dataE_dtau);
}