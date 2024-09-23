#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  int Nt = 1000, Nw = 100, Nb = 250, NbSkip = 50;

  // double da;
  double dg;
  double ds;
  // run 1
  // double gamma_min = 8, gamma_max = 10;
  // double s_min = 0.1, s_max = 0.3;
  double gamma_min = 12, gamma_max = 15.5;
  // double s_min = 0.2, s_max = 0.4;
  // double s_min = 0.4, s_max = 0.7;
  // double s_min = 0, s_max = 0.35;
  double s_min = 0.2, s_max = 0.5;
  double E, sigmaE;
  // int N_gamma = 5;
  int N_gamma = 3;
  int N_s = 5;
  int count;
  dg = (gamma_max - gamma_min) / N_gamma;
  // da = (alpha_max - alpha_min) / N_alpha;
  ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEa.txt", "w");
  parameters_log = fopen("parameters_log_a.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    // ako variramo uz konstantni s
    alpha = alpha_initial;
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
        VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
        fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
      }
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
    alpha = alpha_initial; // za konstantni gamma
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