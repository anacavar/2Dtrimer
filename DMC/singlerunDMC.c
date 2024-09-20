
#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;
double dtau;

int main(int argc, char *argv[])
{
  // int Nt = 1000, Nw = 500, Nb = 100, NbSkip = 0;
  // int Nt = 1000, Nw = 500, Nb = 1050, NbSkip = 50;
  // int Nt = 1000, Nw = 100, Nb = 220, NbSkip = 20;
  // int Nt = 1000, Nw = 300, Nb = 300, NbSkip = 50;
  int Nt = 1000, Nw = 300, Nb = 350, NbSkip = 50;
  // int Nt = 1000, Nw = 100, Nb = 100, NbSkip = 0;
  char batchScript[256];
  snprintf(batchScript, sizeof(batchScript), "C:\\repos\\2Dtrimer\\DMC\\prerunVMC.bat %d %d %d %d", Nt, Nw, Nb, NbSkip);
  system(batchScript);

  double E, sigmaE;
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;
  // double gamma_var = 4.8;
  // double alpha = 4.50;
  // double s = 0.28;
  // dtau = 1.0 * pow(10, -3);
  dtau = 1.0 * pow(10, -5);
  // dtau = 1.0 * pow(10, -5);
  if(argc == 5){
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else if(argc == 1){
    printf("Started with default parameters:\n");
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }
}