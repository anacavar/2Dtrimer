#include <stdio.h>
// #include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main(int argc, char *argv[])
{
  double dg;
  double ds;
  double da;
  double gamma_min = 4.79, gamma_max = 4.89;
  double alpha_min = 4.50, alpha_max = 4.60;
  double s_min = 0.2, s_max = 0.4;
  double E, sigmaE;
  int N_gamma = 5;
  int N_alpha = 10;
  int N_s = 3;
  dg = (gamma_max - gamma_min) / N_gamma;
  da = (alpha_max - alpha_min) / N_alpha;
  ds = (s_max - s_min) / N_s;
  FILE *dataEs, *parameters_log;
  dataEs = fopen("dataEs.txt", "w");
  parameters_log = fopen("parameters_log.txt", "w");

  // if(argc == 5){
  //   int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
  //   printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
  //   // ako variramo uz konstantni s
  //   s = 0.3;
  //   for (int ig = 0; ig < N_gamma; ig++)
  //   {
  //     // za zadani s, za svaki gamma, ili gamma? variraš alpha...
  //     gamma_var = gamma_min + ig * dg;
  //     for (int ia = 0; ia < N_alpha; ia++)
  //     {
  //       alpha = alpha_min + ia * da;
  //       fprintf(parameters_log, "%d\tg=%f\ta=%f\n", ig*N_alpha+ia+1, gamma_var, alpha);
  //       printf("%d. alpha=%f; gamma=%f; s=%f\n", ig*N_alpha+ia+1, alpha, gamma_var, s);
  //       // izvrti program..
  //       // VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
  //       fprintf(dataEs, "%f\t%f\t%f\t%f\n", E, sigmaE, alpha, gamma_var);
  //     }
  //   }
  // }
  if(argc == 1){
    printf("Started with default parameters\n");
    int Nt = 1000, Nw = 100, Nb = 220, NbSkip = 20;
    s = 0.3;
    int count = 0;
    for (int ig = 0; ig < N_gamma; ig++)
    {
      // za zadani s, za svaki gamma, ili gamma? variraš alpha...
      gamma_var = gamma_min + ig * dg;
      for (int ia = 0; ia < N_alpha; ia++)
      {
        count++;

        alpha = alpha_min + ia * da;
        printf("%d*%d+%d+1 = %d != %d. alpha=%f; gamma=%f; s=%f\n", ig, N_gamma, ia, ig*N_gamma+ia+1, count, alpha, gamma_var, s);
        fprintf(parameters_log, "%d.\tg=%f\ta=%f\n", ig*N_gamma+ia+1, gamma_var, alpha);
        // izvrti program..
        // VMC(&E, &sigmaE, Nt, Nw, Nb, NbSkip);
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