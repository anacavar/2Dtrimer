// Difuzijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

// početne vrijednosti
#define Nt 100 // broj koraka
#define Nw 100 // broj šetača
#define Nb 250 // broj blokova
// #define NbSkip 0            // broj prvih blokova koje preskačemo

int main(void)
{
#pragma region // VARIJABLE
  long idum = -1234;
  int ib, it, iw, k;                 // indeks bloka, indeks koraka, indeks šetača, indeks čestice
  double x[4][Nw + 1], y[4][Nw + 1]; // indeks čestice, indeks šetača
#pragma endregion

  FILE *data;
  data = fopen("data.txt", "w");

  // inicijalizacija koordinata čestica gdje je gustoća Psi*Psi znacajna
  for (iw = 1; iw <= Nw; iw++) // po šetačima
  {
    for (k = 1; k <= 3; k++) // po česticama
    {
      // x[k][iw] = (2 * ran1(&idum) - 1) * L0;
      // y[k][iw] = (2 * ran1(&idum) - 1) * L0;
    }
  }

  for (ib = 1; ib <= Nb; ib++) // po blokovima
  {
    for (it = 1; it <= Nt; it++) // po koracima
    {
      for (iw = 1; iw <= Nw; iw++) // po šetačima
      {
        for (k = 1; k <= 3; k++) // po česticama
        {
          // dx = (ran1(&idum) * 2 - 1) * dxyMax;
          // dy = (ran1(&idum) * 2 - 1) * dxyMax;
          // x[k][iw] += dx;
          // y[k][iw] += dy;
        } // kraj petlje čestica
      }   // kraj petlje šetača
    }     // kraj petlje koraka
  }       // kraj petlje blokova
  fclose(data);
  return 0;
}