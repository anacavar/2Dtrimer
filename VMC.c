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
    double x[3][Nw];          // indeks čestice, indeks šetača
    double y[3][Nw];          // indeks čestice, indeks šetača
    double L0 = 10;           // angstromi (10^(−10) m))
    double dx, dy;            // promjene koordinata čestica
    double dxyMax = L0 / 100; // maksimalne promjene x,y,z koordinata
    double x_old[3], y_old[3];
    double random_number;
    double Psi(double);            // probna valna funkcija
    double P[Nw];                  // vjerojatnost prijelaza Ri -> Rf za svakog šetača posebno
    double psi_initial, psi_final; // ukupna inicijalna i finalna probna valna funkcija
    double r12, r23, r13;          // udaljenosti između čestica
    double T;                      // vjerojatnost prijelaza Ri -> Rf
    int accepted = 0, rejected = 0;
    double ratio;
    // double sigma; // "najzgodnije je sigmu mjerit u angstremima"
    // double energija; // "energiju mjerimo u K preko boltzmannove konstante"
#pragma endregion

    FILE *data;
    FILE *koordinate;
    data = fopen("data.txt", "w");
    koordinate = fopen("koordinate.txt", "w");

    // inicijalizacija čestica gdje je gustoća Psi*Psi znacajna
    for (i = 0; i < Nw; i++) // po šetačima
    {
        for (j = 0; j < 3; j++) // po česticama
        {
            x[j][i] = (2 * ran1(&idum) - 1) * L0;
            y[j][i] = (2 * ran1(&idum) - 1) * L0;
        }
        fprintf(koordinate, "%f\t%f\t%f\t%f\t%f\t%f\n", x[0][i], y[0][i], x[1][i], y[1][i], x[2][i], y[2][i]);
        r12 = sqrt(pow((x[0][i] - x[1][i]), 2) + pow((y[0][i] - y[1][i]), 2));
        r23 = sqrt(pow((x[1][i] - x[2][i]), 2) + pow((y[1][i] - y[2][i]), 2));
        r13 = sqrt(pow((x[0][i] - x[2][i]), 2) + pow((y[0][i] - y[2][i]), 2));
        P[i] = Psi(r12) * Psi(r23) * Psi(r13);
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
            psi_initial = P[j];
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
                P[j] = psi_final;
                printf("T: %f\n", T);
            }
            else
            {
                random_number = ran1(&idum);
                if (random_number <= T)
                {
                    // prihvaćamo pomak
                    accepted++;
                    P[j] = psi_final;
                    // printf("T: %f\tRand: %f\n", T, random_number);
                }
                else
                {
                    // odbacujemo pomak
                    rejected++;
                    // printf("T: %f\tRand: %f\n", T, random_number);
                    for (k = 0; k < 3; k++) // po česticama
                    {
                        x[k][j] = x_old[k];
                        y[k][j] = y_old[k];
                    }
                }
            }
        }
    }

    ratio = (double)accepted / (double)(accepted + rejected);
    printf("Prihvacenih: %d\nOdbijenih: %d\nOmjer: %f\n", accepted, rejected, ratio);
    fclose(data);
    fclose(koordinate);
    return 0;
}

// probna valna funkcija
double Psi(double r)
{
    // double alpha = 4.16; // 1;
    double alpha = 4.16; // 1;
    // double alpha = 1; // 1;
    double gamma = 2.82; // 1;
    // double gamma = 1; // 1;
    double s = 0.0027; // 1;
    // double s = 1; // 1;
    return exp(-pow(alpha / r, gamma) - s * r) / r;
}