/***************************************************************************************************
 * Code for encrypting weight t plaintext.
 * Includes random weight t plaintext generator and random permutation.
 * Byte-packing done elsewhere, bits are uint8_t.
 ***************************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "mceliece-big.h"

// encrypts plain to cipher with public key K, cipher = (I|K)*plain
void encrypt(uint8_t plain[CODE_N], uint8_t (*K)[CODE_N - GF_M*GOPPA_T], uint8_t cipher[GF_M*GOPPA_T]) {
    // initialize cipher to zero
    for (int i = 0; i < GF_M*GOPPA_T; i++) {
        cipher[i] = 0;
    }
    // (I|K)*plain = cipher
    for (int i = 0; i < GF_M*GOPPA_T; i++) {
        cipher[i] ^= plain[i];
        for (int j = 0; j < CODE_N - GF_M*GOPPA_T; j++) {
            //printf("%i %i in encrypt\n", i, j);
            cipher[i] ^= K[i][j] & plain[GF_M*GOPPA_T + j];
        }
    }
}

// creates a random {0,1} message of weight t and length n
void random_message(uint8_t plain[CODE_N]) {
    // initialize plain to zero
    for (int i = 0; i < CODE_N; i++) {
        plain[i] = 0;
    }
    // list [0, CODE_N), then permute
    uint16_t L[CODE_N];
    for (int i = 0; i < CODE_N; i++) {
        L[i] = i;
    }
    random_permutation(CODE_N, L);
    // take first segment for weighted positions
    for (int i = 0; i < GOPPA_T; i++) {
        plain[L[i]] = 1;
    }
}

// randomly permute elements of L
void random_permutation(int N, uint16_t L[N]) {
    int j;
    int temp;
    for (int i = N; i > 0; i--) {
        j = rand() % i;
        temp = L[j];
        L[j] = L[i];
        L[i] = temp;
    }
}
