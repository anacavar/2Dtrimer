
#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;
double dtau;

int main(int argc, char *argv[])
{
  // int Nt=500, Nw=150, Nb=300, NbSkip=100;
  int Nt=500, Nw=150, Nb=300, NbSkip=100;
  // int Nt=500, Nw=1000, Nb=550, NbSkip=50;

  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  // system(batchScript); // za prethodno pokretanje VMC, kako bi odgovarali početni parametri šetača

  double E, sigmaE;
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;

  // dtau = 1.0 * pow(10, -6);
  // dtau = 1.0 * pow(10, -5);
  // dtau = 1.0 * pow(10, -4);
  dtau = 1.0 * pow(10, -3);
  // dtau = 0.000005;
  // dtau = 1.0 * pow(10, -2);
  // dtau = 1.0 * pow(10, -1);

  if(argc == 5){
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else if(argc == 1){
    printf("Started with default parameters:\n");
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    printf("sigma=%f, epsilon=%f, alpha=%f, gamma=%f, s=%f\n", sigma, epsilon, alpha, gamma_var, s);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }
}