#include <stdio.h>
#include "VMC.c"
#include "global_vars.h"

double alpha, gamma_var, s, epsilon;

int main(int argc, char *argv[])
{
  int Nt = 1000, Nw = 100, Nb = 250, NbSkip = 50;
  int N_epsilon = 5;
  int count;
  double epsilon_min = 8, epsilon_max = 50;
  double de = (epsilon_max - epsilon_min)/N_epsilon;
  double E, sigmaE, r2, sigmar2;
  alpha = alpha_initial;
  gamma_var = gamma_initial;
  s = s_initial;

  FILE *dataE_epsilon, *parameters_log;
  dataE_epsilon = fopen("dataE_epsilon.txt", "w");
  parameters_log = fopen("parameters_log_epsilon.txt", "w");

  if(argc == 5){
    count=0;
    int Nt = atoi(argv[1]), Nw = atoi(argv[2]), Nb = atoi(argv[3]), NbSkip = atoi(argv[4]);
    printf("%d %d %d %d\n", Nt, Nw, Nb, NbSkip);
    // ako variramo uz konstantni s
    for (int ie = 0; ie < N_epsilon; ie++)
    {
      
      count++;
      epsilon = epsilon_min + ie * de;
      fprintf(parameters_log, "%d\tepsilon=%f\n", count, epsilon);
      printf("%d. epsilon=%f\n", count, epsilon);
      // izvrti program..
      VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
      fprintf(dataE_epsilon, "%f\t%f\t%f\n", E, sigmaE, r2, sigmar2, epsilon);
      
    }
  }
  else if(argc == 1){
    count=0;
    printf("Started with default parameters\n");
    alpha = alpha_initial; // za konstantni gamma
    gamma_var = gamma_initial;
    s = s_initial;
    for (int ie = 0; ie < N_epsilon; ie++)
    {
      count++;
      epsilon = epsilon_min + ie * de;
      fprintf(parameters_log, "%d\tepsilon=%f\n", count, epsilon);
      printf("%d. epsilon=%f\n", count, epsilon);
      // izvrti program..
      VMC(&E, &sigmaE, &r2, &sigmar2, Nt, Nw, Nb, NbSkip);
      fprintf(dataE_epsilon, "%f\t%f\t%f\t%f\t%f\n", E, sigmaE, r2, sigmar2, epsilon);
    }
  }
  else{
    printf("Nt=?; Nw=?; Nb=?; NbSkip=?\n");
  }

  fclose(dataE_epsilon);
  fclose(parameters_log);
}