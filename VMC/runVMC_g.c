#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  int Nt = 500, Nw = 100, Nb = 150, NbSkip = 50;

  epsilon = epsilon_initial;
  // gamma_var = gamma_initial;
  gamma_var = 7.4;

  double ds;
  double da;
  double s_min = 0.6, s_max = 0.75;
  double alpha_min = 9, alpha_max = 9.7;
  double E, sigmaE;
  double r2, sigmar2;
  int N_s = 5;
  int N_alpha = 10;
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
    for (int is = 0; is < N_s; is++)
    {
      // za zadani s, variraš alpha...
      s = s_min + is * ds;
      for (int ia = 0; ia < N_alpha; ia++)
      {
        count++;
        alpha = alpha_min + ia * da;
        fprintf(parameters_log, "%d\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var, s, r2);
      }
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
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
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var, s, r2);
      }
    }
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }

  fclose(dataEs);
  fclose(parameters_log);
}