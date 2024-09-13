#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

    double dg;
    double ds;
    double da;
    double alpha, gamma, s;
    double gamma_min = 4.52, gamma_max = 4.58;
    double alpha_min = 4.57, alpha_max = 4.75;
    double s_min = 0.2, s_max = 0.4;
    double E, sigmaE;
    int N = 5;
    dg = (gamma_max - gamma_min) / N;
    ds = (s_max - s_min) / N;
    da = (alpha_max - alpha_min) / N;

    for (int i=0; i<N; i++){
        gamma = gamma_min + i*dg;
        for(int j=0; j<N; j++){
            alpha = alpha_min + j*da;
            printf("%d. g=%f; a=%f\n",  i*N+j+1, gamma, alpha);
        }  
    }

    return 0;  
}
