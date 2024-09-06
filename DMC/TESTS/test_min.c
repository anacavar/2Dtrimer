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

    int Nw = 20, Nw_min = 15, Nw_max = 25, sum_nw = 0;
    int nw_max_remove, nw_deficit;

    w_pairs W_Rpw[Nw_max]; // statistička težina

    int n_w[Nw_max+1]; // broj potomaka
    int arr[20] = {0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0, 1, 2, 1};

    printf("prije\n");
    for(int i=1; i<=Nw; i++){
        n_w[i] = arr[i-1];
        printf("n_w[%d]=%d\n", i, n_w[i]);
        sum_nw += n_w[i] - 1;
    }

    // inicijalizacija liste radi kasnijeg sortiranja
    int cunt = 0;
    for (int znj = 0; znj < Nw_max; znj++)
    {
        W_Rpw[znj].value = 1.;
        W_Rpw[znj].index = znj;
        printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
        cunt ++;
    }
    printf("total pairs: %d\n", cunt);

    //setiranje vrijednosti
    printf("postavljene vrijednosti:\n");
    for (int znj = 0; znj < Nw; znj++) // PAZI TU INAČE ODE UKURAC POČNE ASSIGNAT OVE KOJE NE POSTOJE...
    {
        W_Rpw[znj].value = (double)n_w[znj+1]; // + 1
        W_Rpw[znj].index = znj;
        printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
    }

    qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // sortira W_Rpw prema values, descending

    printf("sortirano descending (prvo najveći onda najmanji):\n");
    for (int znj = 0; znj < Nw_max; znj++)
    {
        printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
    }



    
    
    // ako šetači padaju na <70%*Nw0, onda +1 optimalne šetače  
    nw_max_remove = Nw - Nw_min; //razlika između trenutnog - 70*Nw0 (minimalnog broja) - koliko ih se maksimalno smije uništiti
    nw_deficit = abs(sum_nw)-nw_max_remove; // deficit šetača do 70%*Nw0

    printf("suma: %d / %d (max remove)\n", sum_nw, nw_max_remove);
    
    if (sum_nw < 0 && abs(sum_nw) > nw_max_remove)
    {
        qsort(W_Rpw, Nw_max, sizeof(w_pairs), compareByValue); // sortira W_Rpw prema values, descending
        for (int znj = 0; znj < Nw_max; znj++)
        {
            printf("W_pairs[%d]: value = %f; index = %d\n", znj, W_Rpw[znj].value, W_Rpw[znj].index);
        }
    
        int count = 0;
        while (count < nw_deficit){
            for (int j = 0; j<Nw; j++){
                if(n_w[W_Rpw[j].index+1]>0){
                    n_w[W_Rpw[j].index+1]++; // ovdje ide +1! tu da, ali ne i u kodu... briem. jer je indeks tamo već pomjeren
                    count++;
                    if(count>=nw_deficit){
                        goto exitLoops;
                    }
                }
            }
        }
        exitLoops:;
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