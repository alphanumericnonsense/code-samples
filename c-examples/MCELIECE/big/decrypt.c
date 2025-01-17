/***************************************************************************************************
 * Code for decryption.
 * Has some bad/redundant code; could be refactored.
 * Byte-packing done elsewhere, bits are uint8_t.
 ***************************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "mceliece-big.h"

void decrypt(uint8_t cipher[GF_M*GOPPA_T], GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t plain[CODE_N]) {

    // generate and restrict the canonical parity-check matrix for g^2 code
    GF g2[2*GOPPA_T + 1];
    g_to_g2(g, g2); // fill g2

    // probably bigger than 8 MB
    // array of pointers?
    GF (*H2)[CODE_N] = malloc(sizeof(GF[2*GOPPA_T][CODE_N]));
    canonical_pcheck2(g2, support, H2);// fill H2 - can improve; only need m*t columns

    // probably bigger than 8 MB
    // array of pointers?
    uint8_t (*H2prime)[CODE_N] = malloc(sizeof(uint8_t[2*GF_M * GOPPA_T][CODE_N]));
    restrict_pcheck2(H2, H2prime); // fill H2prime

    // generate "double syndrome"
    uint8_t c2[2*GF_M*GOPPA_T];
    for (int i = 0; i < 2*GF_M*GOPPA_T; i++) {
        c2[i] = 0;
    }
    // only need first m*t columns
    for (int i = 0; i < 2*GF_M*GOPPA_T; i++) {
        for (int j = 0; j < GF_M*GOPPA_T; j++) {
            c2[i] ^= H2prime[i][j] & cipher[j];
        }
    }

    // unrestrict double syndrome to GF(2^m), S is syndrome polynomial
    GF S[2*GOPPA_T];
    for (int i = 0; i < 2*GOPPA_T; i++) { //initialize to zero
        S[i] = 0;
    }
    for (int i = 0; i < 2*GOPPA_T; i++) {
        for (int j = 0; j < GF_M ; j++) {
            S[i] += (1 << j) * c2[GF_M*i + j]; // bits to GF
        }
    }

    // syndrome decoding
    GF sigma[GOPPA_T + 1];
    decode(S, sigma);
    // switching from inverse roots to roots in error locator
    // easier, avoids conditional for root at infinity
    GF sigma_polar[GOPPA_T + 1];
    for (int i = 0; i < GOPPA_T + 1; i++) {
        sigma_polar[i] = sigma[GOPPA_T - i];
    }
    // evaluate sigma_polar along support
    GF value;
    for (int i = 0; i < CODE_N; i++) {
        value = horner_eval(support[i], sigma_polar, GOPPA_T);
        if (value == 0) {
            plain[i] = 1;
        } else {
            plain[i] = 0;
        }
    }
    // free allocated memory
    free(H2);
    free(H2prime);
}

// fill the canonical parity-check matrix for g^2 code
// redundant code from canonical_pcheck() in pk_gen.c
void canonical_pcheck2(GF g2[2*GOPPA_T + 1], GF support[CODE_N], GF (*H2)[CODE_N]) {
    GF alpha_j;
    GF alpha;
    GF g2_at_alpha_j;
    for (int j = 0; j < CODE_N; j++) {
        alpha_j = support[j];
        alpha = 1;
        g2_at_alpha_j = horner_eval(alpha_j, g2, 2*GOPPA_T);
        for (int i = 0; i < 2*GOPPA_T; i++){
            H2[i][j] = gf_div(alpha, g2_at_alpha_j); // alpha_j^i/g^2(alpha_j)
            alpha = gf_mult(alpha, alpha_j); // alpha_j^i
        }
    }
}

// square the goppa polynomial (sum ai*x^i)^2 = sum ai^2 * x^(2i)
void g_to_g2(GF g[GOPPA_T + 1], GF g2[2*GOPPA_T + 1]) {
    // initialize g2 to 0
    for (int i = 0; i < 2*GOPPA_T + 1; i++) {
        g2[i] = 0;
    }
    // fill g2
    for (int i = 0; i < GOPPA_T + 1; i++) {
        g2[2*i] = gf_mult(g[i], g[i]);
    }
}

// Restricts coefficients of parity-check over GF(2^m) to a binary matrix,
// replacing entries with columns according to minimal polynomial.
// redundant code from canonical_pcheck() in pk_gen.c
void restrict_pcheck2(GF (*H2)[CODE_N], uint8_t (*H2prime)[CODE_N]) {
    GF temp;
    for (int i = 0; i < 2*GOPPA_T; i++) {
        for (int j = 0; j < CODE_N; j++) {
            temp = H2[i][j];
            for (int k = 0; k < GF_M; k++){
                H2prime[GF_M * i + k][j] = temp & 1;
                temp >>= 1;
            }
        }
    }
}

// Syndrome decoding for code from (g^2, support).
// sigma is the error locator polynomial (roots are inverses of elements alpha_j
// in the support whose positions j are where the errors are located).
void decode(GF S[2*GOPPA_T], GF sigma[GOPPA_T + 1]) {
    // initialize variables sigma, beta, l, delta, d
    for (int i = 0; i < GOPPA_T + 1; i++) {
        sigma[i] = 0;
    }
    sigma[0] = 1;
    GF beta[GOPPA_T + 1];
    for (int i = 0; i < GOPPA_T + 1; i++) {
        beta[i] = 0;
    }
    beta[1] = 1;
    int l = 0;
    GF delta = 1;
    GF d = 0;
    GF scalar;
    GF temp[GOPPA_T + 1];
    // main loop
    for (int k = 0; k < 2*GOPPA_T; k++) {
        d = 0;
        for (int i = 0; i < GOPPA_T + 1 && i <= k; i++) {
                d = gf_add(d, gf_mult(sigma[i], S[k - i]));
        }
        scalar = gf_mult(d, gf_inv(delta));
        if (d == 0 || k < 2*l) {
            // sigma
            for (int i = 0; i < GOPPA_T + 1; i++) {
                sigma[i] = gf_add(sigma[i], gf_mult(scalar, beta[i]));
            }
            // beta = x*beta, shift right
            for (int i = GOPPA_T; i > 0; i--) {
                beta[i] = beta[i-1];
            }
            beta[0] = 0;
            l = l;
            delta = delta;
        } else {
            // temp sigma
            for (int i = 0; i < GOPPA_T + 1; i++) {
                temp[i] = sigma[i];
            }
            // sigma
            for (int i = 0; i < GOPPA_T + 1; i++) {
                sigma[i] = gf_add(sigma[i], gf_mult(scalar, beta[i]));
            }
            // beta = x*sigma, shift right
            for (int i = GOPPA_T; i > 0; i--) {
                beta[i] = temp[i-1];
            }
            beta[0] = 0;
            l = k - l + 1;
            delta = d;
        }
    }
}
