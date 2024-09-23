#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  int Nt = 500, Nw = 80, Nb = 80, NbSkip = 30;


  // epsilon = epsilon_initial;
  // epsilon = 8;
  // epsilon = 16.4;
  epsilon = 6;
  // alpha = alpha_initial;
  alpha = 4.3;

  double dg;
  double ds;
  double gamma_min = 4.65, gamma_max = 15.5;
  double s_min = 0.05, s_max = 0.11;
  double E, sigmaE;
  double r2, sigmar2;
  int N_gamma = 1;
  int N_s = 5;
  int count;
  dg = (gamma_max - gamma_min) / N_gamma;
  ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEa.txt", "w");
  parameters_log = fopen("parameters_log_a.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    for (int ig = 0; ig < N_gamma; ig++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
      gamma_var = gamma_min + ig * dg;
      for (int is = 0; is < N_s; is++)
      {
        count++;
        s = s_min + is * ds;
        fprintf(parameters_log, "%d\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var, s);
      }
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
    for (int ig = 0; ig < N_gamma; ig++)
    {
      gamma_var = gamma_min + ig * dg;
      for (int is = 0; is < N_s; is++)
      {
        count++;
        s = s_min + is * ds;
        fprintf(parameters_log, "%d.\tg=%f\ta=%f\ts=%f\n", count, gamma_var, alpha, s);
        printf("%d. alpha=%f; gamma=%f; s=%f\n", count, alpha, gamma_var, s);
        // izvrti program..
        VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
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