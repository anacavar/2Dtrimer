#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double gamma_var, alpha, s, dtau;

int main()
{
  int Nt=500, Nw=100, Nb=220, NbSkip=20;

  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  // system(batchScript); // za prethodno pokretanje VMC, kako bi odgovarali početni parametri šetača

  double d_dtau;
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