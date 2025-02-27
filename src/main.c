#include <stdio.h>
#include "jacob.h"

#define N_SPECIES 10

// static const int n_species = 10;
// static const int n_reactions = 40;

int main() {
    const double t = 0.0;
    const double pres = 101325.0;
    
    // mass fraction H2:O2 phi 1.0
    double y[N_SPECIES] = {
        1000, 0., 0.04030952, 0., 0.95969048, 0., 0., 0., 0., 0.};
    double jac[N_SPECIES*N_SPECIES] = {0.0};

    printf("Calling eval_jacob...\n");
    
    eval_jacob(t, pres, y, jac);

    printf("eval_jacob finished.\n");

    // 結果を表示
    printf("Output values: ");
    for (int i = 0; i < N_SPECIES*N_SPECIES; i++) {
        printf("%e ", jac[i]);
    }
    printf("\n");

    return 0;
}
