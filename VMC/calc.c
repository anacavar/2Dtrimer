#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// konstante u zadanim mjernim jedinicama i njihov preračun u SI
#define k_B 1.                                         // boltzmannova konstanta
#define k_B_si 1.3806 * pow(10, -23)                   // J/K = m2 kg s-2 K-1
#define A 1.                                           // angstrem
#define A_si 1.0 * pow(10, -10)                        // m
#define K 1.                                           // kelvin
#define K_si 1.                                        // kelvin
#define u 1.                                           // jedinica mase
#define u_si 1.6605 * pow(10, -27)                     // kg
#define hbar2 4.851159708 * 10 * k_B *K *pow(A, 2) * u // (reducirana planckova konstanta)^2
#define hbar2_si 1.112121717 * pow(10, -68)            // (Js)^2 = (m2kg/s)^2
// početne vrijednosti
#define sigma 4 * A                                    // angstrema
#define epsilon_initial 1 * k_B *K                     // dubina jame, u kelvinima preko boltzmannove konstante
#define mass 4. * u                                    // u

int main() {
    double x, y, epsilon, rescale_E, rescale_r2, C_epsilon, r_min, r_min3;
    double s, alpha, gamma; // nebitno

    FILE *inputFile, *outputFile;
    inputFile = fopen("final_data.txt", "r");
    outputFile = fopen("data_rescaled.txt", "w");
    
    // Read from the input file, rescale, and write to the output file
    // while (fscanf(inputFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &epsilon, &s, &alpha, &gamma) == 6) {
    for (int i=0; i<9; i++) {

        fscanf(inputFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &epsilon, &s, &alpha, &gamma);
        printf("%f, %f, %f\n", x, y, epsilon);

        C_epsilon = sqrt(2*mass*epsilon/hbar2);
        r_min = pow(2, 1.0/6)*sigma;
        r_min3 = pow(r_min, 3);

        rescale_E = hbar2/(mass*C_epsilon*r_min3);  // Example: rescale first column by A
        rescale_r2 = C_epsilon*r_min3;  // Example: rescale second column by B

        fprintf(outputFile, "%lf %lf\n", x/rescale_E, y/rescale_r2);
    }

    // Close the files
    fclose(inputFile);
    fclose(outputFile);

    printf("Rescaling completed and written to output.txt\n");

    return 0;
}
