
#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  double E, sigmaE;
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;
  if(argc == 5){
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else if(argc == 1){
    printf("Started with default parameters:\n");
    // int Nt = 1000, Nw = 100, Nb = 220, NbSkip = 20;
    // int Nt = 1000, Nw = 300, Nb = 300, NbSkip = 50;
    int Nt = 1000, Nw = 300, Nb = 300, NbSkip = 0;
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    DMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }
}