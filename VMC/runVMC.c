#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  double dg;
  double ds;
  double da;
  double gamma_min = 4.5, gamma_max = 4.6;
  double alpha_min = 4.54, alpha_max = 4.75;
  double s_min = 0.2, s_max = 0.4;
  double E, sigmaE;
  int N = 5;
  dg = (gamma_max - gamma_min) / N;
  ds = (s_max - s_min) / N;
  da = (alpha_max - alpha_min) / N;
  FILE *dataEs;
  dataEs = fopen("dataEs.txt", "w");

  if(argc == 5){
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    // ako variramo uz konstantni s
    s = s_initial;
    for (int ig = 0; ig < N; ig++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
      gamma_var = gamma_min + ig * dg;
      for (int ia = 0; ia < N; ia++)
      {
        alpha = alpha_min + ia * da;
        // izvrti program..
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else if(argc == 1){
    printf("Started with default parameters\n");
    int Nt = 1000, Nw = 100, Nb = 220, NbSkip = 20;
    for (int ig = 0; ig < N; ig++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
      gamma_var = gamma_min + ig * dg;
      for (int ia = 0; ia < N; ia++)
      {
        alpha = alpha_min + ia * da;
        // izvrti program..
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }

  fclose(dataEs);
}