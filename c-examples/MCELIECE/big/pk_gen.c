/***************************************************************************************************
 * Code for generateing public key from secret key.
 * Includes rref over GF(2) and polynomial evaluation (Horner, not AFFT yet).
 *
 * Byte-packing done elsewhere, bits are uint8_t.
 ***************************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include "mceliece-big.h"

// Generate public key as {0,1}-valued uint8_t array
bool pk_gen(GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t (*K)[CODE_N - GF_M*GOPPA_T]) {
    bool good2go = false;

    // larger than 8 MB
    // array of pointers?
    uint8_t (*Hprime)[CODE_N] = malloc(sizeof(uint8_t[GF_M*GOPPA_T][CODE_N]));
    GF (*H)[CODE_N] = malloc(sizeof(GF[GOPPA_T][CODE_N]));

    canonical_pcheck(g, support, H); // generate parity-check matrix
    restrict_pcheck(H, Hprime); // restrict from GF(2^m) to GF(2)
    good2go = binary_rref(GF_M*GOPPA_T, CODE_N, Hprime); // try to reduce to systematic form
    if (good2go) {
        for (int i = 0; i < GF_M*GOPPA_T; i++) {
            for (int j = 0; j < CODE_N - GF_M*GOPPA_T; j++) {
                // public K is right side of rref(Hprime) = (I|K)
                K[i][j] = Hprime[i][j + GF_M*GOPPA_T];
            }
        }
        // free allocated memory
        free(H);
        free(Hprime);
        return true;
    } else {
        // free allocated memory
        free(H);
        free(Hprime);
        return false;
    }
}

// Evaluates degree deg polynomial f at x
GF horner_eval(GF x, GF f[], int deg) {
    GF s = 0;
    for (int i = deg; i > -1; i--) {
        s = gf_mult(s, x);
        s = gf_add(s, f[i]);
    }
    return s;
}

// generate the canonical parity-check matrix for the Goppa code
void canonical_pcheck(GF g[GOPPA_T + 1], GF support[CODE_N], GF (*H)[CODE_N]) {
    GF alpha_j;
    GF alpha;
    GF g_at_alpha_j;
    for (int j = 0; j < CODE_N; j++) {
        //printf("canonical_pcheck loop column %i \n", j);
        alpha_j = support[j];
        alpha = 1;
        g_at_alpha_j = horner_eval(alpha_j, g, GOPPA_T);
        //printf("%i ", g_at_alpha_j);
        for (int i = 0; i < GOPPA_T; i++){
            H[i][j] = gf_div(alpha, g_at_alpha_j); // alpha_j^i/g(alpha_j)
            alpha = gf_mult(alpha, alpha_j); // alpha_j^i
        }
    }
    //printf("exiting canonical_pcheck...\n");
}

// Puts binary M into rref if "systematic" and returns true,
// aborts early if not systematic, i.e. of the form (I|K)
bool binary_rref(int nrows, int ncols, uint8_t (*M)[ncols]) {
    // to flag no pivot in column
    bool zero_pivot;
    // for temporary rows
    uint8_t tempi[ncols];
    uint8_t tempj[ncols];
    // first do upper triangular to speed up non-systematic early abort
    for (int i = 0; i < nrows; i++){
        if (M[i][i] == 0) { // need a non-zero pivot
            zero_pivot = true;
            for (int j = i + 1; j < nrows; j++) { // look below for pivot
                if (M[j][i] == 1) { // switch rows
                    for (int k = 0; k < ncols; k++) {
                        tempi[k] = M[i][k];
                        tempj[k] = M[j][k];
                    }
                    for (int k = 0; k < ncols; k++) {
                        M[i][k] = tempj[k];
                        M[j][k] = tempi[k];
                    }
                    zero_pivot = false;
                    break;
                }
            }
            if (zero_pivot) { // no non-zero pivot, not full-rank
                return false;
            }
        }
        //printf("clearing out below %i\n", i);
        for (int j = i + 1; j < nrows; j++) { // clear out below
            if (M[j][i] == 1) { // add correct multiple of row i to row j
                for (int k = 0; k < ncols; k++) {
                    tempj[k] = M[j][k];
                }
                for (int k = 0; k < ncols; k++) {
                    M[j][k] = (tempj[k] ^ M[i][k]) & 1;
                }
            }
        }
        /*
        // for testing
        if (i > nrows - 10) {
            printf("row %i\n", i);
            for (int k = nrows - 10; k < nrows; k++) {
                for (int l = nrows - 10; l < nrows; l++) {
                    printf("%i ", M[k][l]);
                }
                printf("\n");
            }
        printf("\n");
        }
        */
    }
    // now in systematic upper triangular
    // clear out above
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < i; j++) {
            if (M[j][i] == 1) {
                for (int k = 0; k < ncols; k++) {
                    tempj[k] = M[j][k];
                }
                for (int k = 0; k < ncols; k++) {
                    M[j][k] = (tempj[k] ^ M[i][k]) & 1;
                }
            }
        }
    }
    return true;
}

// Restricts coefficients of parity-check over GF(2^m) to a binary matrix,
// replacing entries with columns according to minimal polynomial
void restrict_pcheck(GF (*H)[CODE_N], uint8_t (*Hprime)[CODE_N]) {
    GF temp;
    for (int i = 0; i < GOPPA_T; i++) {
        for (int j = 0; j < CODE_N; j++) {
            temp = H[i][j];
            for (int k = 0; k < GF_M; k++){
                Hprime[GF_M * i + k][j] = temp & 1;
                temp >>= 1;
            }
        }
    }
}

