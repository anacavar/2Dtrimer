#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double gamma_var, alpha, s, dtau;

int main()
{
  // int Nt=500, Nw=150, Nb=150, NbSkip=50;
  int Nt=500, Nw=100, Nb=200, NbSkip=100;

  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  // system(batchScript); // za prethodno pokretanje VMC, kako bi odgovarali početni parametri šetača

  double d_dtau;
  double E, sigmaE;
  FILE *dataE_dtau;
  dataE_dtau = fopen("dataE_dtau.txt", "w");

  // // U1
  // gamma_var = 4.77;
  // alpha = 4.55;
  // s = 0.3;
  // U2
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;

  int N = 10;
  int nn = 5;
  double dtau_mins[5]={0.000001, 0.00001, 0.0001, 0.001, 0.01};
  double dtau_maxes[5]={0.000011, 0.00011, 0.0011, 0.011, 0.11};
  // double dtau_mins[2]={0.000001, 0.00001};
  // double dtau_maxes[2]={0.000011, 0.00011};

  int count = 0;
  for (int i = 0; i < nn; i++){
    d_dtau = (dtau_maxes[i] - dtau_mins[i]) / N;    for (int idt = 0; idt < N; idt++)
    {
      count++;
      dtau = dtau_mins[i] + idt * d_dtau;
      printf("%d. dtau=%f\n", count, dtau);
      DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
      fprintf(dataE_dtau, "%f\t%f\t%f\n", dtau, E, sigmaE);
    }
  }
  fclose(dataE_dtau);
}