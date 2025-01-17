#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

#include "mceliece-small.h"

// go through rows of K, packing bits to bytes and store in 1-D array.
// tested
void pack_pk(uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T], uint8_t packed_pk[]) {
    int nrows = GF_M*GOPPA_T;
    int ncols = CODE_N - nrows;
    int p_ncols = ncols/8;
    int s; // for holding partial byte
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < p_ncols; j++) {
            s = 0;
            for (int k = 0; k < 8; k++) {
                s += (1 << k) * K[i][8*j + k];
            }
            packed_pk[p_ncols*i + j] = s;
        }
    }
}

// inverse of pack_pk
// tested
void unpack_pk(uint8_t packed_pk[], uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T]) {
    int nrows = GF_M*GOPPA_T;
    int ncols = CODE_N - nrows;
    int p_ncols = ncols/8;
    uint8_t temp;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < p_ncols; j++) {
            temp = packed_pk[p_ncols*i + j];
            for (int k = 0; k < 8; k++) {
                K[i][8*j + k] = temp & 1;
                temp >>= 1;
            }
        }
    }
}

void pack_sk(GF g[], GF support[], uint8_t packed_sk[]) {

}

void unpack_sk(uint8_t packed_sk, GF g[], GF support[]) {

}

void pk_to_file() {

}

void sk_to_file() {

}

void pk_from_file() {

}

void sk_from_file() {

}
