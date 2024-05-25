#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main()
{
  double dg;
  double ds;
  double da;
  double gamma_min = 4.5, gamma_max = 4.6;
  double alpha_min = 4.54, alpha_max = 4.75;
  // double alpha_min = 4.7, alpha_max = 4.8;
  double s_min = 0.2, s_max = 0.4;
  double E, sigmaE;
  int N = 5;
  // int N = 10;
  dg = (gamma_max - gamma_min) / N;
  ds = (s_max - s_min) / N;
  da = (alpha_max - alpha_min) / N;

  // ako variramo uz konstantni s
  FILE *dataEs;
  dataEs = fopen("dataEs.txt", "w");
  s = s_initial;
  for (int ig = 0; ig < N; ig++)
  {
    // za zadani s, za svaki gamma, ili gamma? variraš alpha...
    gamma_var = gamma_min + ig * dg;
    for (int ia = 0; ia < N; ia++)
    {
      alpha = alpha_min + ia * da;
      // izvrti program..
      VMC(&E, &sigmaE);
      fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
    }
  }
  fclose(dataEs);
}