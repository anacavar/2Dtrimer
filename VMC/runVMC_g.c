#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  int Nt = 1000, Nw = 100, Nb = 250, NbSkip = 50;

  // double dg;
  double ds;
  double da;
  // double alpha_min = 4.40, alpha_max = 4.75;
  // double s_min = 0.25, s_max = 0.45;
  // double alpha_min = 8.5, alpha_max = 10;
  // double s_min = 0.10, s_max = 0.30;
  double s_min = 0.35, s_max = 0.55;
  double alpha_min = 8.2, alpha_max = 9;
  double E, sigmaE;
  // int N_s = 5;
  int N_s = 3;
  int N_alpha = 5;
  int count;
  da = (alpha_max - alpha_min) / N_alpha;
  ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEg.txt", "w");
  parameters_log = fopen("parameters_log_g.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    // ako variramo uz konstantni s
    gamma_var = gamma_initial;
    for (int is = 0; is < N_s; is++)
    {
      // za zadani s, variraš alpha...
      // gamma_var = gamma_min + ig * dg;
      s = s_min + is * ds;
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
    gamma_var = gamma_initial; // za konstantni gamma
    for (int is = 0; is < N_s; is++)
    {
      s = s_min + is * ds;
      for (int ia = 0; ia < N_alpha; ia++)
      {
        alpha = alpha_min + ia * da;
        count++;
        fprintf(parameters_log, "%d.\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var, s);
      }
    }
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }

  fclose(dataEs);
  fclose(parameters_log);
}