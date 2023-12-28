// Varijacijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nw 1000 // broj šetača
#define Nk 100  // broj koraka
// #define Nb 100   // broj blokova
// #define Nbskip 0 // broj blokova koje preskačemo radi stabilizacije

float lennardJones(float x1, float x2, float y1, float y2)
{
    float sigma = 1;   // 3.4 * 10^(-10) [m]
    float epsilon = 1; // 1.65 * 10^(-21) [J]
    float Ulj;
    float x = sigma / sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2)); // x = sigma/r
    x = pow(x, 6);                                                 // x = (sigma/r)**6
    if (x <= 1)                                                    // koji broj tu ide i zašto?
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
    float x[3][Nw]; // indeks čestice, indeks šetača - zašto Nw+1?
    float y[3][Nw]; // indeks čestice, indeks šetača
    float L0 = 9;
    float dx, dy;            // promjene koordinata čestica
    float dxyMax = L0 / 100; // maksimalne promjene x,y,z koordinata

#pragma endregion

    FILE *data;
    FILE *koordinate;
    data = fopen("data.txt", "w");
    koordinate = fopen("koordinate.txt", "w");

    // inicijalizacija čestica
    for (i = 0; i < Nw; i++) // po šetačima
    {
        for (j = 0; j < 3; j++) // po česticama
        {
            x[j][i] = ran1(&idum) * L0;
            y[j][i] = ran1(&idum) * L0;
        }
        fprintf(koordinate, "%f\t%f\t%f\t%f\t%f\t%f\n", x[0][i], y[0][i], x[1][i], y[1][i], x[2][i], y[2][i]);
    }

    for (i = 0; i < Nk; i++) // po koracima
    {
        for (j = 0; j < Nw; j++) // po šetačima
        {
            for (k = 0; k < 3; k++) // po česticama
            {
                dx = (ran1(&idum) * 2 - 1) * dxyMax;
                dy = (ran1(&idum) * 2 - 1) * dxyMax;
                x[k][j] += dx;
                y[k][j] += dy;
            }
        }
        fprintf(data, "%f\t%f\t%f\t%f\t%f\t%f\n", x[0][0], y[0][0], x[1][0], y[1][0], x[2][0], y[2][0]);
    }

    fclose(data);
    fclose(koordinate);
    return 0;
}