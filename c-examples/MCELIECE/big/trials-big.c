/*********************************************************************************************
 * For larger parameter sets, arrays for public key need more than 8 MB and must be
 * explicitly allocated.  Various functions are rewritten to accommadate this.
 * *******************************************************************************************/


#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

#include "mceliece-big.h"

int main () {

    srand(time(NULL));
    int N = 50;
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("* Generate secret/public key and encrypt/decrypt %i random messages.\n", N);
    printf("* Parameters m = 13, t = 128, n = 6688 (big).\n");
    printf("* Expect 0 - 5 key rejections.\n");

    bool key_gen_success = false; // for key_gen loop
    //registers
    GF g[GOPPA_T + 1];
    GF support[CODE_N];

    // memory allocation larger than 8 MB
    // array of pointers?
    uint8_t (*K)[CODE_N - GF_M*GOPPA_T] = malloc(sizeof(uint8_t[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T]));

    while (!key_gen_success) {
        // create secret key
        sk_gen(g, support);
        //printf("* any key to continue\n");
        //getchar();
        // create public key from secret key
        key_gen_success = pk_gen(g, support, K);
        if (!key_gen_success) {
            printf("*\tKeys rejected...\n");
        }
    }
    printf("* Key pair successfully created.\n");

    uint8_t message[CODE_N];
    for (int i = 0; i < N; i++) {
        random_message(message);
        uint8_t c[GF_M*GOPPA_T];
        encrypt(message, K, c);
        uint8_t plain[CODE_N];
        decrypt(c, g, support, plain);

        int counter = 0;
        for (int j = 0; j < CODE_N; j++) {
            if (plain[j] != message[j]) {
                if (counter == 0) {
                    printf("*\tTrial %i: Failure :'(\n", i);
                }
                printf("*\t\tFailed bit (j, alpha_j) = (%i, %i)\n", j, support[j]);
                counter++;
            }
        }
        if (counter == 0) {
            printf("*\tTrial %i: Success!\n", i);
        }
    }
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");

    free(K);
    return 0;
}
//-----------------------------------------//
// debugging bullshit
//-----------------------------------------//
    /*
    // testing single encryption/decryption test vector
    printf("\nEncrypt/decrypt 10 random messages with Goppa poly and support G_TEST, L_TEST (in constants.c) and report success/failure:\n");
    GF g[GOPPA_T + 1];
    GF support[CODE_N];
    for (int i = 0; i < GOPPA_T + 1; i++) {
        g[i] = G_TEST[i];
    }
    for (int i = 0; i < CODE_N; i++) {
        support[i] = L_TEST[i];
    }
    uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T];
    pk_gen(g, support, K);

    // random plaintext
    uint8_t message[CODE_N];
    for (int i = 0; i < 10; i++) {
        random_message(message);

        uint8_t c[GF_M*GOPPA_T];
        encrypt(message, K, c);

        uint8_t plain[CODE_N];
        decrypt(c, g, support, plain);

        int counter = 0;
        for (int j = 0; j < CODE_N; j++) {
            if (plain[j] != message[j]) {
                //printf("decryption failure at %i!\n", i);
                counter++;
            }
        }
        if (counter == 0) {
            printf("\t%i Success!\n", i);
        } else {
            printf("\t%i Failures: %i\n", i, counter);
        }
    }
    printf("End program.\n\n");
    */

    /* // fixed test plaintext
    u i*nt8_t c[GF_M*GOPPA_T];
    encrypt(PLAIN_TEST, K, c);

    uint8_t plain[CODE_N];
    decrypt(c, g, support, plain);

    for (int j = 0; j < CODE_N; j++) {
        printf("%i %i\n", plain[j], PLAIN_TEST[j]);
        if (plain[j] != PLAIN_TEST[j]) {
            printf("decryption failure at %i!\n", j);
    }
    }
    */
    /*
    u i*nt16_t x = 1234;
    uint8_t xbits[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (int i = 0; i<16; i++) {
        xbits[i] = x&1;
        x>>=1;
    }
    int s = 0;
    for (int i = 0; i<16; i++) {
        s += (1<<i) * xbits[i];
    }
    printf("%i\n", s);
    */

    /*
    uint8_t c[GF_M*GOPPA_T];
    encrypt(PLAIN_TEST, K, c);

    // testing cipher, works!
    for (int j = 0; j < GF_M*GOPPA_T; j++) {
        printf("%i %i\n", c[j], CIPHER_TEST[j]);
        if (c[j] != CIPHER_TEST[j]) {
            printf("cipher failure at %i!\n", j);
        }
    }

    uint8_t plain[CODE_N];
    decrypt(c, g, support, plain);

    for (int j = 0; j < CODE_N; j++) {
        if (plain[j] != PLAIN_TEST[j]) {
            printf("decryption failure at %i!\n", j);
        }
    }
    */

    /* // testing gf_rref()
    GF M[3][3] = {{124,765,343},{565,1234,568},{1,2,3}};
    bool OK;
    OK = gf_rref(3,3,M);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%i ", M[i][j]);
        }
        printf("\n");
    }
    */
    /* // testing permute_field
    printf("RAND_MAX = %i\n", RAND_MAX);
    GF field[1<<GF_M];
    for (int i =0; i < 1 << GF_M; i++) {
        field[i] = i;
    }
    permute_field(field);
    for (int i =0; i < 1 << GF_M; i++) {
        printf("%i\n", field[i]);
    }
    */

    /* // testing pk_gen
    // the right answer
    uint8_t answer[10][10] =
    {{1,0,0,0,0,1,0,0,1,0},
    {0,1,1,0,0,1,0,1,1,1},
    {0,1,0,1,0,1,1,1,0,0},
    {1,1,0,1,1,0,0,0,0,1},
    {1,0,0,0,1,1,0,0,1,1},
    {1,1,0,0,1,1,0,1,0,0},
    {0,0,1,1,1,0,0,0,0,0},
    {0,1,0,1,1,1,0,1,0,1},
    {1,1,1,1,0,0,1,0,0,0},
    {0,0,0,1,1,0,1,1,0,1}};
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            printf("%i ", answer[i][j]);
        }
        printf("\n");
    }

    //rref says...

    GF g[GOPPA_T + 1];
    GF support[CODE_N];

    for (int i = 0; i < GOPPA_T + 1; i++) {
        g[i] = G_TEST[i];
    }

    for (int i = 0; i < CODE_N; i++) {
        support[i] = L_TEST[i];
    }

    uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T];

    pk_gen(g, support, K);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j<10; j++) {
            printf("%i ", K[i][j]);
        }
        printf("\n");
    }
    */

    /* // testing binary_rref
    uint8_t M[3][6] = {{1,0,1,1,0,0},{0,1,0,0,1,0},{1,1,0,0,0,1}};
    binary_rref(3,6,M);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j<6; j++) {
            printf("%i ", M[i][j]);
        }
        printf("\n");
    }
    */
