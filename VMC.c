// Varijacijski Monte Carlo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

#define Nw 10 // broj šetača
#define Nk 10 // broj koraka

int main(void)
{
#pragma region // KONSTANTE & VARIJABLE
    long idum = -1234;
    int i, j, k, ib;
    float x[3][Nw]; // indeks čestice, indeks šetača
    float y[3][Nw]; // indeks čestice, indeks šetača
    float L0 = 9;
    float dx, dy;            // promjene koordinata čestica
    float dxyMax = L0 / 100; // maksimalne promjene x,y,z koordinata
    float V[Nw];             // potencijal između 3 čestice šetača u poslijednjem koraku ?
    float delta_V, random_number;
    float x_old[3], y_old[3];
    float Psi(float);             // probna valna funkcija
    float psi_initial, psi_final; // ukupna inicijalna i finalna probna valna funkcija
    float r12, r23, r13;          // udaljenosti između čestica
    float T;                      // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    float ratio;
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
                x_old[k] = x[k][j];
                y_old[k] = y[k][j];
                x[k][j] += dx;
                y[k][j] += dy;
            }
            // METROPOLIS ALGORITAM
            // T = |P(Rf)|^2/|P(Ri)|^2; za svakog šetača posebno
            r12 = sqrt(pow((x_old[0] - x_old[1]), 2) + pow((y_old[0] - y_old[1]), 2));
            r23 = sqrt(pow((x_old[1] - x_old[2]), 2) + pow((y_old[1] - y_old[2]), 2));
            r13 = sqrt(pow((x_old[0] - x_old[2]), 2) + pow((y_old[0] - y_old[2]), 2));
            psi_initial = Psi(r12) * Psi(r23) * Psi(r13);
            r12 = sqrt(pow((x[0][j] - x[1][j]), 2) + pow((y[0][j] - y[1][j]), 2));
            r23 = sqrt(pow((x[1][j] - x[2][j]), 2) + pow((y[1][j] - y[2][j]), 2));
            r13 = sqrt(pow((x[0][j] - x[2][j]), 2) + pow((y[0][j] - y[2][j]), 2));
            psi_final = Psi(r12) * Psi(r23) * Psi(r13);
            printf("r12: %f\tr23: %f\tr13: %f\n", r12, r23, r13);
            printf("Psi(R12): %f\n", Psi(r12));
            printf("Psi_initial: %f\tPsi_final: %f\n", psi_initial, psi_final);
            T = psi_final * psi_final / psi_initial * psi_initial;
            if (T > 1)
            {
                // prihvaćamo pomak
                accepted++;
            }
            else
            {
                random_number = ran1(&idum);
                if (random_number < T) // je li sigurno ovako ili je možda obrnuto?
                {
                    // prihvaćamo pomak
                    accepted++;
                }
                else
                {
                    // odbacujemo pomak
                    rejected++;
                    printf("T: %f\tRand: %f\n", T, random_number);
                    for (k = 0; k < 3; k++) // po česticama
                    {
                        x[k][j] = x_old[k];
                        y[k][j] = y_old[k];
                    }
                }
            }
        }
    }

    ratio = (float)accepted / (float)(accepted + rejected);
    printf("Prihvacenih: %d\nOdbijenih: %d\nOmjer: %f\n", accepted, rejected, ratio);

    fclose(data);
    fclose(koordinate);
    return 0;
}

// probna valna funkcija
float Psi(float r)
{
    // šta je s ovima? otkud te varijable? brijem da ću ih varirat u VMCu
    float alpha = 1;
    float gamma = 1;
    float s = 1;
    // jel koristim ovu probnu ili neku drugu? PITAJ PERU
    return exp(-pow(alpha / r, gamma) - s * r) / r;
}