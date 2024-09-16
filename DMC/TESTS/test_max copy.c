// TESTIRAMO LISTE KOJE PRELAZE BOUNDS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){

    // int Nw = 20, Nw_min = 15, Nw_max = 25, deltaNw, sum_nw = 0, plus_minus, count;
    int Nw = 390, Nw_min = 210, Nw_max = 390, deltaNw, sum_nw = 0, plus_minus, count;
    float remainder;
    int deltaNw_initial = Nw_max - Nw;
    deltaNw = deltaNw_initial;

    int n_w[Nw_max+1], n_w_temp[Nw_max+1]; // broj potomaka
    // int arr[20] = {6, 9, 5, 4, 3, 5, 1, 2, 6, 6, 4, 0, 3, 0, 3, 2, 0, 1, 2, 0};
    int arr[390] = {[0 ... 389] = 1};

    arr[68] = 0;
    arr[191] = 2;
    arr[193] = 2;

    int oduzimanje, dodavanje;
    oduzimanje = 0;
    dodavanje = 0;

    // printf("prije\n");
    for(int i=1; i<=Nw; i++){
        n_w[i] = arr[i-1];
        // printf("n_w[%d]=%d\n", i, n_w[i]);
        if(n_w[i]==0){
            // deltaNw++;
            oduzimanje++;
        }
        else{
            // sum_nw += n_w[i] - 1;
            dodavanje += n_w[i] - 1;
        }
    }
    
    printf("oduzimanje=%d\tdodavanje=%d\tdeltaNw=%d\n", oduzimanje, dodavanje, deltaNw);

    // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
    if (dodavanje > deltaNw+oduzimanje) // što za sum_nw < -Nw_max ?
    {
        // VJEROJATNO TU TREBA NEKA SOFISTICIRANIJA METODA KOJA ĆE UZETI U OBZIR KOJE WALKERE TREBA PRIORITIZIRATI!!
        count = 0; 
        for (int iw = 1; iw <= Nw; iw++) 
        { 
            if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
            {
                plus_minus = -1; 
                count += plus_minus;
            }
            if(n_w[iw] > 0){
                plus_minus = (int)((n_w[iw] - 1)*(deltaNw+oduzimanje) / dodavanje); // n_w' / slobodno = n_w / suma
                remainder = (float)(((n_w[iw] - 1)*(deltaNw+oduzimanje)) % dodavanje )/(float)dodavanje;
                if (remainder >= 0.5){
                    plus_minus += 1; // round up
                }
                count += plus_minus; // jer ovo može bit negativno ako je previše nula a svi ovi ostali se vrate u 1
                if (count > deltaNw+oduzimanje) 
                {
                    printf("count: %d\n", count);
                    plus_minus -= count - (deltaNw+oduzimanje); // oduzmemo razliku
                    count = deltaNw+oduzimanje; // spustimo na maskimalnu dozvoljenu popunjenost
                }
            }
            //dodat još jednu iteraciju dok se count sigurno ne popuni...
            n_w[iw] = 1 + plus_minus; 
        }
        int first_iteration = 1;
        // tbh mislim da jako rijetko ako ikad ulazi u ove petlje...
        while(count < deltaNw+oduzimanje){
            if(first_iteration == 1){
                for(int iw = 1; iw <= Nw; iw++){
                    if(n_w[iw]>1){
                        n_w[iw]++;
                        count++;
                    }
                    if (count>=deltaNw+oduzimanje){
                        break;
                    }
                }
                first_iteration=0;
            }
            else if(first_iteration=0 && count<deltaNw+oduzimanje){
                for(int iw = 1; iw <= Nw; iw++){
                    if(n_w[iw]>0){
                        n_w[iw]++;
                        count++;
                    }
                    if (count>=deltaNw+oduzimanje){
                        break;
                    }
                }
            }
        }
    }

    // // (>130%) ako je suma dodanih šetača veća od preostalog slobodnog prostora u listi (>130%)
    if (sum_nw > deltaNw) // što za sum_nw < -Nw_max ?
    {
        count = 0; // <= deltaNw
        for (int iw = 1; iw <= Nw; iw++) // nije dobro testirano
        { 
            if(n_w[iw] == 0) // svi koji idu u nulu moraju ići u nulu
            {
                plus_minus = -1; 
            }

            if(n_w[iw] > 0){
                plus_minus = (int)((n_w[iw] - 1)*deltaNw / sum_nw); // n_w' / slobodno = n_w / suma
                remainder = (float)(((n_w[iw] - 1)*deltaNw) % sum_nw )/(float)sum_nw;
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
            n_w[iw] = 1 + plus_minus; 
        }
        while(count < deltaNw){
          for(int i = 1; i <= Nw; i++){
            if(n_w[i]>1){
              n_w[i]++;
              count++;
            }
            if (count>=deltaNw){
              break;
            }
          }
        }
    }

    printf("smanjena suma\n");
    sum_nw = 0;
    for(int i=1; i<=Nw; i++){
        printf("n_w[%d]=%d\n", i, n_w[i]);
        if(n_w[i]!=0){
            sum_nw += n_w[i] - 1;
        }
    }

    printf("sum_nw poslije: %d\n", sum_nw);

    // // korak 11. - kopiranje potomaka R'_p(w) u novi ansambl
    // // memcpy(x_temp, x, sizeof(x));
    // // memcpy(y_temp, y, sizeof(y));
    // // memcpy(E_L_temp, E_L, sizeof(E_L));
    // // memcpy(n_w_temp, n_w, sizeof(n_w));
    // memcpy(n_w_temp, n_w, sizeof(n_w));
    // int Nw_temp = Nw;
    // Nw = 0;
    // int indeks = 1;
    // // printf("Nw_temp=%d\n", Nw_temp);
    // for (int iw = 1; iw <= Nw_temp; iw++)
    // {
    //     // printf("iw=%d\n", iw);
    //     if (n_w_temp[iw] != 0) // ako je n = 0 onda taj šetač biva uništen (preskačemo ga)
    //     {
    //         for (int in = 0; in < n_w_temp[iw]; in++)
    //         {
    //             n_w[indeks + in] = n_w_temp[iw];
    //         } 
    //     }
    //     // if (iw != Nw_temp) // ako nije zadnji korak
    //     // {
    //     // printf("Nw_temp=%d\tiw=%d\n", Nw_temp, iw);
    //     printf("indeks=%d + nw_temp[%d]=%d\n", indeks, iw, n_w_temp[iw]);
    //     indeks += n_w_temp[iw];
    //     // }
    //     // else if (iw == Nw_temp && n_w_temp[iw] == 0) // ako je zadnji korak
    //     // if (iw == Nw_temp) // ako je zadnji korak
    //     // {
    //     //     // printf("Nw_temp=%d\tiw=%d\n", Nw_temp, iw);
    //     //     // indeks += n_w_temp[iw];
    //     //     indeks = indeks - 1; // brijem ide -1 jer kao ne trebam se dalje pomaknut od zadnjeg indeksa
    //     // }
    // }
    // Nw = indeks-1; 

    printf("Nw=%d\n", Nw);

    printf("nakon preraspodijele\n");
    for(int i=1; i<=Nw; i++){
        printf("n_w[%d]=%d\n", i, n_w[i]);
        if(n_w[i]!=0){
            sum_nw += n_w[i] - 1;
        }
    }
  
  return 0;  
}

// što bude ako je onaj bug kojeg sam ispravila?
// što bude ako na listama nije +1?