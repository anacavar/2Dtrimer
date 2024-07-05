
#include <stdio.h>
#include "DMC.c"
#include "global_vars.h"

double alpha, gamma_var, s;

int main()
{
  double E, sigmaE;
  gamma_var = gamma_initial;
  alpha = alpha_initial;
  s = s_initial;
  DMC(&E, &sigmaE);
}