#include <stdio.h>
#include <stdint.h>
#include "mceliece-big.h"

// Field addition
GF gf_add(GF a, GF b) {
    return a ^ b;
}

// Field multiplication, most significant exponent
GF gf_mult(GF a, GF b) {
    GF s = 0;
    int eps = 0;
    for (int i = GF_M - 1; i > -1; i--) {
        eps = (s >> (GF_M - 1)) & 1;
        s = ((s << 1) & MASK) ^ (eps * REDUCER) ^ (((b >> i) & 1) * a);
    }
    return s & MASK;
}

// Field division, assumes denominator is non-zero
GF gf_div(GF a, GF b) {
    return gf_mult(a, gf_inv(b));
}

/******************************************
* Use Fermat inverse to create lookup table
* for field inverses.
*******************************************/
GF gf_fermat_inv(GF a) {
    return gf_pow(a, (1 << GF_M) - 2);
}

void create_inverse_table(GF* inverse_table) {
    for (int i = 0; i < (1 << GF_M); i++) {
        inverse_table[i] = gf_fermat_inv(i);
    }
}
// code used to generate table
//GF INVERSE_TABLE[1 << GF_M];
//create_inverse_table(INVERSE_TABLE);
/*****************************************/

// Field inverses using LUT
GF gf_inv(GF a) {
    if (a == 0) {
        printf("inverting zero!\n");
        return 0;
    } else {
        return INVERSE_TABLE[a];
    }
}

// Field exponentiation
GF gf_pow(GF a, int e) {
    if (e == 0){
        return 1;
    } else if (e == 1) {
        return a;
    } else if (e < 0) {
        return gf_pow(gf_fermat_inv(a), -e);//return gf_pow(gf_inv(a), -e);
    } else {
        if (e % 2 == 0) {
            return gf_pow(gf_mult(a, a), e/2);
        } else {
            return gf_mult(a, gf_pow(gf_mult(a, a), e/2));
        }
    }
}
