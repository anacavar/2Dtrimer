// Varijacijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nw 1000 // broj šetača
#define Nk 50   // broj koraka
// #define Nb 100   // broj blokova
// #define Nbskip 0 // broj blokova koje preskačemo radi stabilizacije

float lennardJones(float x1, float x2, float y1, float y2)
{
    float sigma = 1;   // 3.4 * 10^(-10) [m]
    float epsilon = 1; // 1.65 * 10^(-21) [J]
    float Ulj;
    float x = sigma / sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2)); // x = sigma/r
    x = pow(x, 6);                                                 // x = (sigma/r)**6
    if (x <= 2)
        Ulj = (pow(x, 2) - x);
    else
        Ulj = 0;
    Ulj = 4 * epsilon * Ulj;
    return Ulj;
}

float potencijal(float *x, float *y)
{
    float U = 0;
    for (int i = 1; i <= 3; i++)
    {
        for (int j = i + 1; j <= 3; j++)
        {
            U += lennardJones(x[i], x[j], y[i], y[j]);
        }
    }
    return U;
}

int main(void)
{
#pragma region // KONSTANTE & VARIJABLE
    long idum = -1234;
    int i, j, k, ib;
    // konstante:
    // float sigma = 1;   // 3.4 * 10^(-10) [m]
    // float epsilon = 1; // 1.65 * 10^(-21) [J]
    // veličine:

    // svaki šetač je jedna distribucija čestica
    float x[3][Nw + 1]; // indeks čestice, indeks šetača - zašto Nw+1?
    float y[3][Nw + 1]; // indeks čestice, indeks šetača
    float L0 = 9;

    // // ili da inicijaliziram početne pozicije čestica pa onda pustim da se simulacija sama stabilizira u trokut prema potencijalu?
    // float x_trimer[3] = {L0 / 3, 2 * L0 / 3, L0 / 2};                 // trimer
    // float y_trimer[3] = {L0 / 3, L0 / 3, L0 / 3 * (1 + sqrt(3) / 2)}; // trimer L0/3 + korijen(3)/2 * L0/3
    // // možda su šetači pozicije od po tri čestice???
    // // ŠETAČI SU POZCIJE OD PO TRI ČESTICE!!!!

#pragma endregion

    FILE *data;
    FILE *koordinate;
    data = fopen("data.txt", "w");
    koordinate = fopen("koordinate.txt", "w");

    // inicijalizacija čestica
    for (i = 1; i <= Nw; i++) // po šetačima
    {
        for (j = 0; j < 3; j++) // po česticama
        {
            x[j][i] = ran1(&idum) * L0;
            y[j][i] = ran1(&idum) * L0;
        }
        fprintf(koordinate, "%f\t%f\t%f\t%f\t%f\t%f\n", x[0][i], y[0][i], x[1][i], y[1][i], x[2][i], y[2][i]);
    }

    fclose(data);
    fclose(koordinate);
    return 0;
}