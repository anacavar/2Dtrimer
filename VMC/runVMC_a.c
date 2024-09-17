#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  double dg;
  // double da;
  double ds;
  double gamma_min = 4.68, gamma_max = 4.92;
  // double alpha_min = 4.40, alpha_max = 4.70;
  double s_min = 0.2, s_max = 0.4;
  double E, sigmaE;
  int N_gamma = 5;
  int N_alpha = 5;
  int N_s = 5;
  int count;
  dg = (gamma_max - gamma_min) / N_gamma;
  // da = (alpha_max - alpha_min) / N_alpha;
  ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEs.txt", "w");
  parameters_log = fopen("parameters_log.txt", "w");

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
    int Nt = 1000, Nw = 100, Nb = 220, NbSkip = 20;
    alpha = alpha_initial; // za konstantni gamma
    for (int ig = 0; ig < N_s; ig++)
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