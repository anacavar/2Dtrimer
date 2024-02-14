// TRIMER HELIJA Josip
// Primjer koda za HFDB potencijal i VMC metoda.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran1.c"

const long int Ns = 1000;  // broj koraka
const long int Nw = 5000;  // broj setaca
const long int Nb = 100;   // broj blokova
const long int NbSkip = 1; // broj prvih blokova koje preskacemo

const double Dk = 6059.6499694573; // D je konstanta h^2/(2*m*amu*k_b)
// https://www.wolframalpha.com/input?i=+%286.582119569e-
16 % 29 % 5E2 % 2F % 282 * 4.002602 * % 28931.49410242e6 % 2F % 28299792458 * 1e10 % 29 % 5E2 % 29 * 8.61 7333262e-8 % 29
    // kad se uvrste sve mjerne jedinice i pokrate, energiju dobijemo u mK.
    // HFDB konstante
    const double A = 184431.01;
const double epsilon = 10948; // mK
const double rm = 2.963;      // Å
const double alfa = 10.43329537;
const double beta = -2.27965105;
const double C6 = 1.36745214;
const double C8 = 0.42123807;
const double C10 = 0.17473318;
const double D = 1.4826;
// varijacijski parametri
double g;
double s;
double a;
void fprovjera(FILE *dat)
{
  if (dat == NULL)
  {
    printf("Error while opening the file.\n");
    exit(EXIT_FAILURE);
  }
}
// pomoćna fje za potencijal
double F(double x)
{
  if (x <= D)
  {
    return exp(-pow(((D / x) - 1.), 2));
  }
  else
  {
    return 1.0;
  }
}
// hfdb potencijal
double V(double r)
{
  return epsilon * (A * exp(-alfa * (r / rm) + beta * pow(r / rm, 2)) - F(r / rm) * (C6 / pow(r / rm, 6) +
                                                                                     C8 / pow(r / rm, 8) + C10 / pow(r / rm, 10)));
}
// valna f2(r)
20

    double
    Psi(double r)
{

  return (exp(-pow((a / r), g) - s * r)) / r;
}

double fdr(double r)
{

  return (-s * r + g * pow((a / r), g) - 1.0) / pow(r, 2);
}

double fddr(double r)
{

  return (-2.0 * s * r - 1.000 + (-pow(g, 2) + g) * pow((a / r), g)) / pow(r, 2);
}
// lokalna energija (4.96)
double energija(double r[4][4], double x[4][4][4])
{
  double Vuk = (V(r[1][2]) + V(r[1][3]) + V(r[2][3])); // potencijal po svim parovima

  return Vuk -
         Dk * (fddr(r[1][2]) + fddr(r[1][3]) + (pow((fdr(r[1][2]) * x[1][2][1] + fdr(r[1][3]) * x[1][3][1]), 2) + pow((fdr(r[1][2]) * x[1][2][2] + fdr(r[1][3]) * x[1][3][2]), 2) + pow((fdr(r[1][2]) * x[1][2][3] + fdr(r[1][3]) * x[1][3][3]), 2)) + fddr(r[2][1]) + fddr(r[2][3]) + (pow((fdr(r[2][1]) * x[2][1][1] + fdr(r[2][3]) * x[2][3][1]), 2) + pow((fdr(r[2][1]) * x[2][1][2] + fdr(r[2][3]) * x[2][3][2]), 2) + pow((fdr(r[2][1]) * x[2][1][3] + fdr(r[2][3]) * x[2][3][3]), 2)) + fddr(r[3][1]) + fddr(r[3][2]) + (pow((fdr(r[3][1]) * x[3][1][1] + fdr(r[3][2]) * x[3][2][1]), 2) + pow((fdr(r[3][1]) * x[3][1][2] + fdr(r[3][2]) * x[3][2][2]), 2) + pow((fdr(r[3][1]) * x[3][1][3] + fdr(r[3][2]) * x[3][2][3]), 2)));
}

#include "VMChfdb3.c"
int main()
{
  double grange[] = {4.14, 4.18}, dg;  // g
  double arange[] = {2.82, 2.83}, da;  // a
  double srange[] = {0.025, 0.03}, ds; // s

  double E[2]; // E i sigma E (varijacija)
  int N = 1;   // N^3 parametara

  FILE *param;
  param = fopen("parametri2.txt", "w");
  fprovjera(param);

  // koraci u parametarskom prostoru
  dg = (grange[1] - grange[0]) / N;
  da = (arange[1] - arange[0]) / N;
  ds = (srange[1] - srange[0]) / N;

  int ig, ia, i, j, is;

  for (is = 0; is < N; is++)
  {
    s = srange[0] + is * ds;

    for (ig = 0; ig < N; ig++)
    {

      g = grange[0] + ig * dg;

      for (ia = 0; ia < N; ia++)
      {
        a = arange[0] + ia * da;

        VMC(E);
        fprintf(param, "%e\t%e\t%e\t%e\t%e\n", g, a, s, E[0], E[1]);
      }
    }
  }
  fclose(param);
}
KOD ZA VODIK PRILAGOĐEN ZA TRIMER HELIJA
    /* ====================================================================
       _______________________ ______  __   ______________________  ___
       __  ___/___  __/__  __ \___  / / /   __  ___/____  _/___   |/  /
       _____ \ __  /   _  / / /__  /_/ /    _____ \  __  /  __  /|_/ /
       ____/ / _  /    / /_/ / _  __  /     ____/ / __/ /   _  /  / /
       /____/  /_/     \____/  /_/ /_/      /____/  /___/   /_/  /_/
          Prirodoslovno-matematicki fakultet u Splitu
          Stohasticke simulacije u klasicnoj i kvantnoj fizici
          Varijacijski Monte Carlo :: Osnovno stanje vodikova atoma
          2018/2019, 2019/2020
    Koristene oznake:
          Nw    = broj setaca (number of walkers)
          iw    = indeks setaca (index of walker)
          Nb    = broj blokova (number of blocks)
          ib    = indeks bloka (index of block)
          Nt    = broj koraka (number of time-steps)
          it    = indeks vremenskog koraka (time-step)
          k     = koordinata (coordinate)
          NbSkip= broj preskocenih blokova
          E     = lokalna energija (local energy)
          SwE   = suma (srednjih E) po setacima
          StE   = suma (srednjih E) po koracima
          SbE   = suma (srednjih E) po blokovima
          accept= brojac prihvacenih koraka
          acc_ib= udio prihvacenih koraka
          x     = koordinate posa setaca
          dx    = promjena koordinate nasumicno od -dk do dk
          dk    = maksimalna duljina koraka
          xp    = koordinate probnog posa setaca
          r1    = modul radijvektora probnog posa
          r2    = r1*r1
          P     = vjerojatnost nalazenja na posu x
          Pp    = vjerojatnost nalazenja u probnom posu
          T     = vjerojatnost prijelaza x -> xp
          Psi   = valna funkcija
    ======================================================================= */

    void
    VMC(double *E_ret)
{
  long idum = (-455);
  int is, ib, iw, k, itmp, at, d, rdis, i, kut10, n = 40000;
  double accept, acc_is, AE, sigmaE, kut;
  double x[4][4][4], xp[4][4][4], dk[4], E[Nw + 1], P[Nw + 1],
      pos[4][4][Nw], newpos[4][4];
  double dx, r[4][4], rnew[4][4], SwE, SsE, SbE, SbE2, Pp, T;
  float *distribucija;
  float *kutdis;
  distribucija = (float *)malloc(n * sizeof(float));
  kutdis = (float *)malloc(n * sizeof(float)); // distribucija kuta
  accept = 0.;                                 // prihvacanje
  dk[1] = dk[2] = dk[3] = 2.1;                 // maksimalne promjene koordinata
  FILE *fout;
  FILE *dis;
  FILE *tocke;
  dis = fopen("distribucija-kuta-duljina.dat", "w");
  fprovjera(dis);
  tocke = fopen("distribucijatocaka.dat", "w");
  fprovjera(tocke);

  char ime[70];

  sprintf(ime, "E_hfdb_g_%.3lf_a_%.3lf_s_%.3lf_.dat", g, a, s);
  fout = fopen(ime, "w"); // datoteka za pohranu srednjih vrijednosti

  for (i = 0; i < 40000; i++)
  {
    distribucija[i] = 0; // početne vrij.
    kutdis[i] = 0;
  }

  // inicijalizacija polozaja gdje je gustoca Psi*Psi znacajna
  for (iw = 1; iw <= Nw; iw++)
  {
    r[1][2] = 0.;
    r[1][3] = 0.;
    r[2][3] = 0.;
    r[2][1] = 0.;
    r[3][1] = 0.;
    r[3][2] = 0.;
    pos[1][1][iw] = 0.;
    pos[2][1][iw] = 0.;
    pos[1][2][iw] = 0.;
    pos[2][2][iw] = 0.;
    pos[1][3][iw] = 0.;
    pos[2][3][iw] = 0.;

    for (at = 1; at <= 3; at++)
    { // petlje po atomima

      for (k = 1; k <= 3; k++)
      { // petlja po koordinatama

        pos[at][k][iw] = 150.0 * ((ran1(&idum)) - 0.5); // poč.koordinate
      }
    }

    for (k = 1; k <= 3; k++)
    {

      // dvočestične komponenete. x[a][d][k] gdje je [a] = prvi atom [d]= drugi atom,
      [k] - komponenta x, y, z.

                             x[1][2][k] = pos[2][k][iw] - pos[1][k][iw]; // x2-x1, y2-y1 i z2-z1
      x[1][3][k] = pos[3][k][iw] - pos[1][k][iw];
      x[2][3][k] = pos[3][k][iw] - pos[2][k][iw];
      x[2][1][k] = pos[1][k][iw] - pos[2][k][iw];
      x[3][1][k] = pos[1][k][iw] - pos[3][k][iw];
      x[3][2][k] = pos[2][k][iw] - pos[3][k][iw];
      // r12=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) i analogno za r23 i r13
      r[1][2] += pow(x[1][2][k], 2);
      r[1][3] += pow(x[1][3][k], 2);
      r[2][3] += pow(x[2][3][k], 2);
      r[2][1] += pow(x[2][1][k], 2);
      r[3][1] += pow(x[3][1][k], 2);
      r[3][2] += pow(x[3][2][k], 2);
    }
    // magnituda duljina vektora r12, r23 i r13.
    r[1][2] = sqrt(r[1][2]);
    r[1][3] = sqrt(r[1][3]);
    r[2][3] = sqrt(r[2][3]);
    r[2][1] = sqrt(r[2][1]);
    r[3][1] = sqrt(r[3][1]);
    r[3][2] = sqrt(r[3][2]);

    P[iw] = pow(Psi(r[1][2]) * Psi(r[1][3]) * Psi(r[2][3]), 2); // vjerojatnost
    E[iw] = energija(r, x);                                     // lokalna energija
  }
  SbE = 0.;
  SbE2 = 0.;
  for (ib = 1; ib <= Nb; ib++) // blokovi
  {
    SsE = 0.;

    for (is = 1; is <= Ns; is++) // koraci
    {
      SwE = 0;

      for (iw = 1; iw <= Nw; iw++) // setaci
      {
        rnew[1][2] = 0.;
        rnew[1][3] = 0.;
        rnew[2][3] = 0.;
        rnew[2][1] = 0.;
        rnew[3][1] = 0.;
        rnew[3][2] = 0.;
        newpos[1][1] = 0.;
        newpos[2][1] = 0.;
        newpos[1][2] = 0.;
        newpos[2][2] = 0.;
        newpos[1][3] = 0.;
        newpos[2][3] = 0.;

        for (at = 1; at <= 3; at++)
        { // petlje po parovima atomima ukupno ih
          je 3. r12, r13 i r23

              for (k = 1; k <= 3; k++)
          {

            newpos[at][k] = pos[at][k][iw] + dk[k] * (2.0 * (ran1(&idum) - 0.5));
          }
        }

        for (k = 1; k <= 3; k++)
        {

          // međučestične komponenete.INDEX x[a][d][k] gdje je [a] = prvi atom [d]= drugi
          atom, [k] - komponenta x, y, z.

                                       xp[1][2][k] = newpos[2][k] - newpos[1][k]; // x2-x1, y2-y1 i z2-z1
          xp[1][3][k] = newpos[3][k] - newpos[1][k];
          xp[2][3][k] = newpos[3][k] - newpos[2][k];
          xp[2][1][k] = newpos[1][k] - newpos[2][k];
          xp[3][1][k] = newpos[1][k] - newpos[3][k];
          xp[3][2][k] = newpos[2][k] - newpos[3][k];
          // r12=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) i analogno za r23 i r13
          rnew[1][2] += pow(xp[1][2][k], 2);
          rnew[1][3] += pow(xp[1][3][k], 2);
          rnew[2][3] += pow(xp[2][3][k], 2);
          rnew[2][1] += pow(xp[2][1][k], 2);
          rnew[3][1] += pow(xp[3][1][k], 2);
          rnew[3][2] += pow(xp[3][2][k], 2);
        }
        // magnituda duljina vektora r12, r23 i r13.
        rnew[1][2] = sqrt(rnew[1][2]);
        rnew[1][3] = sqrt(rnew[1][3]);
        rnew[2][3] = sqrt(rnew[2][3]);
        rnew[2][1] = sqrt(rnew[2][1]);
        rnew[3][1] = sqrt(rnew[3][1]);
        rnew[3][2] = sqrt(rnew[3][2]);

        Pp = pow(Psi(rnew[1][2]) * Psi(rnew[1][3]) * Psi(rnew[2][3]), 2); // vjerojatnost u
        probnom pol
            T = Pp / P[iw]; // prijelaz

        // Metropolis algoritam
        if (T >= 1)
        {
          for (at = 1; at <= 3; at++)
          {
            for (k = 1; k <= 3; k++)
            {
              pos[at][k][iw] = newpos[at][k];
            }
          }

          if (ib > 0)
          {
            rdis = (int)(rnew[1][2] * 100.); // pretvorba duljine*100 u cijeli broj
            distribucija[rdis] += 1.0;       // brojač dodaj +1 za taj index duljine
            // račun kuta iz skalarnog produkta

            kut = acos((xp[1][2][1] * xp[1][3][1] + xp[1][2][2] * xp[1][3][2] + xp[1][2][3] * xp[1][3][3]) / (rnew[1][2] * rnew[1][3]));
            kut10 = (int)(kut * 100.);
            kutdis[kut10] += 1.0; // brojač kuteva
            fprintf(tocke, "%16.8e %16.8e %16.8e\n", rnew[1][2],
                    rnew[1][3], Pp); // duljine parova atoma i vjerojatnost.
          }

          accept += 1.;

          P[iw] = Pp;
          E[iw] = energija(rnew, xp);
        }
        else if (ran1(&idum) <= T)
        {
          for (at = 1; at <= 3; at++)
          {
            for (k = 1; k <= 3; k++)
            {
              pos[at][k][iw] = newpos[at][k];
            }
          }

          if (ib > 0)
          {
            rdis = (int)(rnew[1][2] * 100.);
            distribucija[rdis] += 1.0;

            kut = acos((xp[1][2][1] * xp[1][3][1] + xp[1][2][2] * xp[1][3][2] + xp[1][2][3] * xp[1][3][3]) / (rnew[1][2] * rnew[1][3]));
            kut10 = (int)(kut * 100.);
            kutdis[kut10] += 1.0;
            fprintf(tocke, "%16.8e %16.8e %16.8e\n", rnew[1][2],
                    rnew[1][3], Pp);
          }
          accept += 1.;
          P[iw] = Pp;
          E[iw] = energija(rnew, xp);
        }

        SwE = SwE + E[iw];
      } // setaci
      if (is % 100 == 0)
      {
        acc_is = accept / ((ib - 1) * Nw * Ns + is * Nw);
        if (acc_is > 0.5)
        {
          for (k = 1; k <= 3; k++)
          {
            dk[k] = dk[k] * 1.05;
          }
        }
        if (acc_is < 0.5)
        {
          for (k = 1; k <= 3; k++)
          {
            dk[k] = dk[k] * 0.95;
          }
        }
      }
      // akumulacija podataka nakon stabilizacije

      if (ib > NbSkip)
      {
        SsE += SwE / Nw;
      }
    } // koraci

    if (ib > NbSkip) // akumulacija podataka nakon stabilizacije
    {
      SbE += SsE / Ns;
      SbE2 += SsE * SsE / (Ns * Ns);

      fprintf(fout, "%7d %16.8e %16.8e\n", ib - NbSkip, SsE / Ns, SbE / (ib - NbSkip));
    }
    itmp = (int)(round(acc_is * 100.));

    printf("%6d. blok:  %d%% prihvacenih,  Eb = %10.2e\n", ib - NbSkip,
           itmp, SsE / Ns);
  } // blokovi

  for (i = 0; i < 40000; i++)
  {

    fprintf(dis, "%f\t%f\t%lf\t%lf\n", (float)(i) * 180.0 / (100. * M_PI), (float)(i) / 100.0, distribucija[i], kutdis[i]);
  }

  AE = SbE / (Nb - NbSkip);
  sigmaE = sqrt(fabs((SbE2 / (Nb - NbSkip) - AE * AE) / ((Nb - NbSkip) - 1.)));

  E_ret[0] = AE;
  E_ret[1] = sigmaE;
  accept = accept / (Nw * Ns * Nb);

  printf("postotak prihvacenih koraka: %4.1f\n", accept * 100.);
  printf("\n konacni max. koraci: %6.2e %6.2e %6.2e\n", dk[1], dk[2],
         dk[3]);
  printf("\n g: %6.2e  a %6.2e s %6.2e\n", g, a, s);
  printf("\n E = %8.5e +- %6.2e \n\n", AE, sigmaE);
  fclose(fout);
  fclose(dis);
}