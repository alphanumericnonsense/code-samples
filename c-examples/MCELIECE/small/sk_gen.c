/***************************************************************************************************
 * Code for generateing secret key.
 * Includes rref over GF(2^m), permutation of GF(2^m), some polynomial arithmetic in GF(2^m)[x]
 *
 *
 ***************************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "mceliece-small.h"


void sk_gen(GF g[GOPPA_T + 1], GF support[CODE_N]) {
    //printf("in sk_gen\n");
    // randomly permute the field and
    // take first segment for support
    GF field[1<<GF_M];
    for (int i = 0; i < (1 << GF_M); i++) {
        field[i] = i;
    }
    permute_field(field);
    for (int i = 0; i < CODE_N; i++) {
        support[i] = field[i];
    }

    // generate irreducible Goppa polynomial through
    // rejection sampling of random element in degree t extension
    // which is usually a generator
    bool irred = false;
    GF R[GOPPA_T][GOPPA_T + 1];
    GF r[GOPPA_T];
    GF rpow[GOPPA_T]; // holds powers of r modulo EXTRA_MIN_POLY
    GF temp[GOPPA_T]; // register for product r*rpow
    while (!irred) {
        //printf("Attempting sk_gen...\n");
        // fill r with random field elements
        for (int i = 0; i < GOPPA_T; i++) {
            r[i] = rand() % (1 << GF_M);
        }
        //set rpow = r
        for (int i = 0; i < GOPPA_T; i++) {
            rpow[i] = r[i];
        }
        // create R, the rpower matrix
        //first column
        for (int i = 0; i < GOPPA_T; i++) {
            R[i][0] = 0;
        }
        R[GOPPA_T - 1][0] = 1;
        // rest of columns
        for (int i = 1; i < GOPPA_T + 1; i++) { // columns
            for (int j = 0; j < GOPPA_T; j++) { // rows
                R[j][i] = rpow[GOPPA_T - 1 - j];
            }
            // temp = rpow * r
            poly_mod_prod(r, GOPPA_T - 1, rpow, GOPPA_T - 1, EXTRA_MIN_POLY, GOPPA_T, temp);
            // rpow = temp
            for (int j = 0; j < GOPPA_T; j++) {
                rpow[j] = temp[j];
            }
        }
        /*
        // print left part of R
        for (int k = 0; k < GOPPA_T; k++) {
            for (int l = 0; l < 10; l++) {
                printf("%i ", R[k][l]);
            }
            printf("\n");
        }
        */
        irred = gf_rref(GOPPA_T, GOPPA_T + 1, R);
        /*
        // print reduced R
        for (int k = 0; k < GOPPA_T; k++) {
            for (int l = 0; l < GOPPA_T + 1; l++) {
                printf("%i ", R[k][l]);
            }
            printf("\n");
        }
        */
    }
    // Goppa poly is last column of R
    for (int i = 0; i < GOPPA_T; i++) {
        g[i] = R[i][GOPPA_T];
    }
    g[GOPPA_T] = 1; // g is monic
}

// random permutation of [0, 2^m)
// not quite uniform due to RNG
void permute_field(GF L[1 << GF_M]) {
    int j;
    GF temp;
    for (int i = (1 << GF_M) - 1; i > 0; i--) {
        j = rand() % (i + 1);
        temp = L[j];
        L[j] = L[i];
        L[i] = temp;
    }
}

// multiplies a, b in GF(2^m)[x] modulo f.
// most significant exponent with reduction after each multiplication by x.
// degree of output is deg(f) - 1.
// assumes f is monic!
void poly_mod_prod(GF a[], int d_a, GF b[], int d_b, GF f[], int d_f, GF result[]) {
    //initialize temporary result to zero
    GF temp[d_f + 1];
    for (int i = 0; i < d_f + 1; i++) {
        temp[i] = 0;
    }
    /*
    for (int l = 0; l < d_f+1; l++) {
        printf("%i ", temp[l]);
    }
    printf("\n");
    */
    // x^d_f = reducer
    GF reducer[d_f];
    for (int i = 0; i < d_f; i++) {
        reducer[i] = f[i];
    }
    // MSE with reduction
    for (int i = d_b; i > -1; i--) {
        // multiply by x, i.e. right shift
        for (int j = d_f; j > 0; j--) {
            temp[j] = temp[j-1];
        }
        temp[0] = 0;
        /*
        for (int l = 0; l < d_f+1; l++) {
            printf("%i ", temp[l]);
        }
        printf("\n");
        */
        // reduce x^d_f term if it is non-zero
        if (temp[d_f] != 0) {
            for (int j = 0; j < d_f; j++) {
                temp[j] = gf_add(temp[j], gf_mult(temp[d_f], reducer[j]));
            }
            temp[d_f] = 0;
            /*
            for (int l = 0; l < d_f+1; l++) {
                printf("%i ", temp[l]);
            }
            printf("\n");
            */
        }
        // add b[i]*a to temp
        for (int j = 0; j < d_a + 1; j++) {
            temp[j] = gf_add(temp[j], gf_mult(a[j], b[i]));
        }
        /*
        for (int l = 0; l < d_f+1; l++) {
            printf("%i ", temp[l]);
        }
        printf("\n");
        */
    }
    // fill result from temp
    for (int i = 0; i < d_f; i++) {
        result[i] = temp[i];
    }
}



// multiplication in GF(2^m)[x]
// most significant exponent
void poly_mult(GF a[], int d_a, GF b[], int d_b, GF c[]) {
    int d_c = d_a + d_b; // degree of the product
    for (int i = 0; i < d_c + 1; i++) {
        c[i] = 0;
    }
    for (int i = d_b + 1; i > -1; i--) {
        // multiply by x, i.e. shift
        for (int j = d_c; j > 0; j--) {
            c[j] = c[j-1];
        }
        c[0] = 0;
        // add b[i]*a to c
        for (int j = 0; j < d_a + 1; j++) {
            c[j] = gf_add(c[j], gf_mult(b[i], a[j]));
        }
    }
}

// Puts GF(2^m) M into rref if "systematic" and returns true,
// aborts early if not systematic, i.e. not of the form (I|M0)
bool gf_rref(int nrows, int ncols, GF M[nrows][ncols]) {
    //printf("in gf_rref\n");
    // to flag no pivot in column
    bool zero_pivot;
    // for some scaling
    GF scalar;
    // for temporary rows
    GF tempi[ncols];
    GF tempj[ncols];
    // first do upper triangular to speed up non-systematic early abort
    for (int i = 0; i < nrows; i++){
        if (M[i][i] == 0) { // need a non-zero pivot
            zero_pivot = true;
            for (int j = i + 1; i < nrows; j++) { // look below for pivot
                if (M[j][i] != 0) { // switch rows
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
        // scale row
        scalar = M[i][i];
        for (int j = i; j < ncols; j++) {
            M[i][j] = gf_mult(gf_inv(scalar), M[i][j]);
        }
        // clearing out below
        for (int j = i + 1; j < nrows; j++) {
            // add correct multiple of row i to row j
            if (M[j][i] != 0) {
                for (int k = 0; k < ncols; k++) {
                    tempj[k] = M[j][k];
                }
                for (int k = i; k < ncols; k++) {
                    M[j][k] = gf_add(tempj[k], gf_mult(M[i][k], tempj[i]));
                }
            }
        }
    }
    // now in systematic upper triangular
    // clear out above
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < i; j++) {
            if (M[j][i] != 0) {
                for (int k = 0; k < ncols; k++) {
                    tempj[k] = M[j][k];
                }
                for (int k = i; k < ncols; k++) {
                    M[j][k] = gf_add(tempj[k], gf_mult(M[i][k], tempj[i]));
                }
            }
        }
    }
    return true;
}
