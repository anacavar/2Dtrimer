#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  int Nt = 500, Nw = 150, Nb = 150, NbSkip = 30;

  // epsilon = epsilon_initial;
  epsilon = 6;
  s = 0.082;

  double dg;
  double da;

  double alpha_min = 4.32, alpha_max = 8.9;
  double gamma_min = 4.6, gamma_max = 4.7;
  double E, sigmaE;
  double r2, sigmar2;
  int N_alpha = 1;
  int N_gamma = 6;
  int count;
  dg = (gamma_max - gamma_min) / N_gamma;
  da = (alpha_max - alpha_min) / N_alpha;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEs.txt", "w");
  parameters_log = fopen("parameters_log_s.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    for (int ia = 0; ia < N_alpha; ia++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
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
    for (int ig = 0; ig < N_gamma; ig++)
    {
      gamma_var = gamma_min + ig * dg;
      for (int ia = 0; ia < N_alpha; ia++)
      {
        alpha = alpha_min + ia * da;
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