#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  int Nt = 1000, Nw = 100, Nb = 250, NbSkip = 50;

  double dg;
  double da;
  // double ds;
  // double gamma_min = 4.67, gamma_max = 5.00;
  // double alpha_min = 4.48, alpha_max = 4.63;
  // double gamma_min = 8, gamma_max = 10;
  double gamma_min = 9, gamma_max = 15;
  double alpha_min = 7.8, alpha_max = 8.9;
  // double s_min = 0.2, s_max = 0.4;
  double E, sigmaE;
  // int N_gamma = 5;
  int N_gamma = 5;
  int N_alpha = 10;
  // int N_alpha = 10;
  // int N_s = 5;
  int count;
  dg = (gamma_max - gamma_min) / N_gamma;
  da = (alpha_max - alpha_min) / N_alpha;
  // ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEs.txt", "w");
  parameters_log = fopen("parameters_log_s.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    // ako variramo uz konstantni s
    s = s_initial;
    for (int ig = 0; ig < N_gamma; ig++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
      gamma_var = gamma_min + ig * dg;
      for (int ia = 0; ia < N_alpha; ia++)
      {
        count++;
        alpha = alpha_min + ia * da;
        fprintf(parameters_log, "%d\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
    s = s_initial; // za konstantni s
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
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
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