#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  double E, sigmaE;
  double r2, sigmar2;
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;
  epsilon = epsilon_initial;

  // int Nt=500, Nw=50, Nb=200, NbSkip=100;
  // int Nt=500, Nw=80, Nb=200, NbSkip=100;
  // int Nt=500, Nw=100, Nb=200, NbSkip=100;
  // int Nt=500, Nw=120, Nb=200, NbSkip=100;
  int Nt=500, Nw=150, Nb=200, NbSkip=100;

  if(argc == 5){
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
  }
  else if(argc == 1){
    printf("Started with default parameters:\n");
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    printf("sigma=%f, epsilon=%f, alpha=%f, gamma=%f, s=%f\n", sigma, epsilon, alpha, gamma_var, s);
    VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }
}