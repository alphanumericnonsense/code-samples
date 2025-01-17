#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

#include "mceliece-small.h"

int main () {

    srand(time(NULL));
    int N = 10;

    bool key_gen_success = false; // for key_gen loop
    FILE *pkptr;
    FILE *skptr;
    char pkstr[64];
    char skstr[64];
    //registers
    GF g[GOPPA_T + 1];
    GF support[CODE_N];
    uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T];

    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("* Generating %i key pairs and writing to file.\n", N);
    printf("* Parameters m = %i, t = %i, n = %i (small).\n", GF_M, GOPPA_T, CODE_N);
    printf("* \n");

    for (int i = 0; i < N; i++){
        key_gen_success = false;
        printf("*\tGenerating keys (expect 0 - 5 key rejections):\n");
        while (!key_gen_success) {
            // create secret key
            sk_gen(g, support);

            // create public key from secret key
            key_gen_success = pk_gen(g, support, K);
            if (!key_gen_success) {
                printf("*\t\tKeys rejected...\n");
            }
        }
        printf("*\n");
        printf("*\tKey pair successfully created.\n");

        printf("*\tWriting to pk%i.txt, sk%i.txt\n", i, i);
        printf("*\n");

        sprintf(pkstr, "small-keys/pk%i.txt", i);
        sprintf(skstr, "small-keys/sk%i.txt", i);

        pkptr = fopen(pkstr, "w");
        skptr = fopen(skstr, "w");

        for (int j = 0; j < GOPPA_T*GF_M; j++){
            for (int k = 0; k < CODE_N - GF_M*GOPPA_T; k++){
                //printf("*\tpk loop %i %i\n", j, k);
                fprintf(pkptr, "%i", K[j][k]);
            }
            fprintf(pkptr, "\n");
        }

        fprintf(skptr, "goppa polynomial degree %i\n", GOPPA_T);
        for (int j = 0; j <= GOPPA_T; j++){
            fprintf(skptr, "%i\n", g[j]);
        }
        fprintf(skptr, "support length %i\n", CODE_N);
        for (int j = 0; j < CODE_N; j++){
            fprintf(skptr, "%i\n", support[j]);
        }

        fclose(pkptr);
        fclose(skptr);
    }
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return 0;
}
