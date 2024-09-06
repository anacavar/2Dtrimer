// TESTIRAMO LISTE KOJE PRELAZE BOUNDS

#include <stdio.h>
#include <stdlib.h>

int main(){

    int Nw = 20, Nw_min = 15, Nw_max = 25, deltaNw, sum_nw = 0, plus_minus, count;
    float remainder;

    deltaNw = Nw_max - Nw;

    int n_w[Nw_max+1]; // broj potomaka
    int arr[20] = {0, 5, 1, 4, 3, 1, 1, 2, 0, 0, 4, 0, 3, 0, 3, 2, 0, 1, 2, 4};

    printf("prije\n");
    for(int i=1; i<=Nw; i++){
        n_w[i] = arr[i-1];
        printf("n_w[%d]=%d\n", i, n_w[i]);
        sum_nw += n_w[i] - 1;
    }

    printf("%d\n", sum_nw);


    // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
    if (sum_nw > deltaNw) // što za sum_nw < -Nw_max ?
    {
        count = 0; // <= deltaNw
        for (int iw = 1; iw <= Nw; iw++) // nije dobro testirano
        { 
            if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
            {
                plus_minus = -1; 
                count += plus_minus; // ?
            }

            if(n_w[iw] > 0){
                plus_minus = (int)((n_w[iw] - 1)*deltaNw / sum_nw); // n_w' / slobodno = n_w / suma
                remainder = ((n_w[iw] - 1)*deltaNw) % sum_nw;
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
            printf("prije: n_w[%d]=%d\n", iw, n_w[iw]);
            n_w[iw] = 1 + plus_minus; 
            printf("poslije: n_w[%d]=%d\n", iw, n_w[iw]);
        }
    }

    sum_nw = 0;

    printf("poslije\n");
    for(int i=1; i<=Nw; i++){

        printf("n_w[%d]=%d\n", i, n_w[i]);
        sum_nw += n_w[i] - 1;
    }

    printf("%d\n", sum_nw);
  
  return 0;  
}

// što bude ako je onaj bug kojeg sam ispravila?
// što bude ako na listama nije +1?