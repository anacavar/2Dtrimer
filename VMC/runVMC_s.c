#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  int Nt = 500, Nw = 100, Nb = 150, NbSkip = 50;

  epsilon = epsilon_initial;
  s = s_initial;

  double dg;
  double da;

  double alpha_min = 4.1, alpha_max = 5.1;
  double gamma_min = 3.5, gamma_max = 6.0;
  double E, sigmaE;
  double r2, sigmar2;
  int N_alpha = 5;
  int N_gamma = 10;
  int count;
  da = (alpha_max - alpha_min) / N_alpha;
  dg = (gamma_max - gamma_min) / N_gamma;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEs.txt", "w");
  parameters_log = fopen("parameters_log_s.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    for (int ia = 0; ia < N_alpha; ia++)
    {
      alpha = alpha_min + ia * da;
      for (int ig = 0; ig < N_gamma; ig++)
      {
        count++;
        gamma_var = gamma_min + ig * dg;
        fprintf(parameters_log, "%d\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
    for (int ia = 0; ia < N_alpha; ia++)
    {
      alpha = alpha_min + ia * da;
      for (int ig = 0; ig < N_gamma; ig++)
      {
        gamma_var = gamma_min + ig * dg;
        count++;
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        fprintf(parameters_log, "%d.\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        // izvrti program..
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }

  fclose(dataEs);
  fclose(parameters_log);
}