// TESTIRAMO LISTE KOJE PRELAZE BOUNDS

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
  double value;
  int index;
} w_pairs;

int compareByValue(const void *a, const void *b)
{
  return ((w_pairs *)b)->value - ((w_pairs *)a)->value; // descending - ..., 3, 2, 1, 0
}

int main(){

    int Nw = 20, Nw_min = 15, Nw_max = 25, deltaNw, sum_nw = 0, plus_minus, count;
    float remainder;

    w_pairs W_Rpw[Nw_max]; // statistička težina


    deltaNw = Nw_max - Nw;

    int n_w[Nw_max+1]; // broj potomaka
    int arr[20] = {0, 5, 1, 4, 3, 1, 1, 2, 0, 0, 4, 0, 3, 0, 3, 2, 0, 1, 2, 4};

    // printf("prije\n");
    for(int i=1; i<=Nw; i++){
        n_w[i] = arr[i-1];
        // printf("n_w[%d]=%d\n", i, n_w[i]);
        sum_nw += n_w[i] - 1;
    }

    printf("suma nadodanih = %d / Nw_max-Nw=%d\n", sum_nw, Nw_max-Nw);

    // inicijalizacija liste radi kasnijeg sortiranja
    int cunt = 0;
    for (int znj = 0; znj < Nw_max; znj++)
    {
        W_Rpw[znj].value = 1.;
        W_Rpw[znj].index = znj;
        cunt ++;
    }

    //setiranje vrijednosti
    printf("postavljene vrijednosti:\n");
    for (int znj = 0; znj < Nw; znj++) // PAZI TU INAČE ODE UKURAC POČNE ASSIGNAT OVE KOJE NE POSTOJE...
    {
        W_Rpw[znj].value = (double)n_w[znj+1]; // + 1
        W_Rpw[znj].index = znj;
        // printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
    }

    qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // sortira W_Rpw prema values, descending

    // printf("nakon sortiranja:\n");
    // for (int znj = 0; znj < Nw; znj++) // PAZI TU INAČE ODE UKURAC POČNE ASSIGNAT OVE KOJE NE POSTOJE...
    // {
    //     printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
    // }


    // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
    if (sum_nw > deltaNw) // što za sum_nw < -Nw_max ?
    {
        // VJEROJATNO TU TREBA NEKA SOFISTICIRANIJA METODA KOJA ĆE UZETI U OBZIR KOJE WALKERE TREBA PRIORITIZIRATI!!
        count = 0; // <= deltaNw
        // for (int iw = 1; iw <= Nw; iw++) // TU STAVI DA SE ŠEĆE PO SORTIRANOJ LISTI, A NE PO OVAK I ONAK
        for (int iw = 0; iw < Nw; iw++) // TU STAVI DA SE ŠEĆE PO SORTIRANOJ LISTI, A NE PO OVAK I ONAK
        { 
            // if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
            if(W_Rpw[iw].value == 0) // svi koji idu u nulu moraju ići u nulu
            {
                plus_minus = -1; 
                count += plus_minus; // ?
            }

            // if(n_w[iw] > 0){
            if(W_Rpw[iw].value > 0){
                // (5-1)*5/16 = 20/16 = 5/4= 1
                plus_minus = (int)(((int)W_Rpw[iw].value - 1)*deltaNw / sum_nw); // n_w' / slobodno = n_w / suma
                printf("%d\t%d\t%d\t%d\n", (int)W_Rpw[iw].value, plus_minus, deltaNw, sum_nw);

                // remainder = ((n_w[iw] - 1)*deltaNw) % sum_nw;
                remainder = (((int)W_Rpw[iw].value - 1)*deltaNw) % sum_nw;
                if (remainder >= 0.5){
                    plus_minus += 1; // round up
                }
                count += plus_minus;
                if (count > deltaNw) 
                {
                    plus_minus -= count - deltaNw; // oduzmemo razliku
                    count = deltaNw; // spustimo na maskimalnu dozvoljenu popunjenost
                }
            }
            // n_w[iw] = 1 + plus_minus; 
            W_Rpw[iw].value = 1 + plus_minus; 
        }
    }

    sum_nw = 0;

    // printf("poslije\n");
    // for(int i=1; i<=Nw; i++){
    //     printf("n_w[%d]=%d\n", i, n_w[i]);
    //     sum_nw += n_w[i] - 1;
    // }

    printf("poslije\n");
    for(int i=0; i<Nw; i++){
        printf("W_rp[%d].value=%d\n", i, W_Rpw[i].value);
        sum_nw += W_Rpw[i].value - 1;
    }

    printf("sum_nw poslije: %d / Nw_max-Nw:%d\n", sum_nw, Nw_max-Nw);
  
  return 0;  
}

// što bude ako je onaj bug kojeg sam ispravila?
// što bude ako na listama nije +1?