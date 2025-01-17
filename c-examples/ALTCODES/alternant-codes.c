// General purpose alternant decoder.
//
// Alternant codes are defined by:
// integers q, m, n, r < n
// a list of n non-zero elements (h_0, ..., h_{n-1}) in GF(q^m)^n
// a list of n distinct elements (alpha_0, ..., alpha_{n-1}) in GF(q^m)^n
//
// Over GF(q^m), we have a canonical parity-check matrix H = (h_j*alpha_j^i)_{i=0 to r-1, j=0 to n}
// which has rank r.  We then take the GF(q) subfield subcode, getting a parity check matrix H'
// over GF(q) by replacing entries of H with columns using a GF(q)-linear isomorphism GF(q^m) \equiv GF(q)^m,
// deleting redundant rows.
//
// We have inequalities n - m*r < k <= n - r and d >= r + 1
// for the dimension and minimum distance of the code.
//
// In the random codes constructed, we force k = n - m*r
// by imposing the restricted parity-check matrix reduces to (I|K)
// through rejection sampling (for convenience to control code parameters).
//
// These codes include the classes of Reed-Solomon, BCH, and Goppa codes.
//
// The decoders below formulate the error locator as (1-alpha_{i_0}*X)*...*(1-alpha_{i_t}*X)
// where there errors are at indices i_k, 0 <= k < t.  However, this makes detecting an error at alpha = 0
// difficult unless you know the value of t (so you can see the degree is too small by one).
// This is fine for BCH codes where alpha != 0, but the alternant construction allows alpha = 0.
// Couldn't get around this yet, so excluded alpha = 0 (thought I could just reverse the syndrome).

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

/*******************************************************
 * Code parameters need to be cajoled so that
 * the random code construction actually works,
 * i.e. the FLD_DIM*CODE_R - by - CODE_N binary matrix
 * needs to rref to (I|K) with some probability (asymptotically 1/3ish)
 *
 * Need 2^FLD_DIM > CODE_N and
 * CODE_N > CODE_R*FLD_DIM.
 *
 * Corrects less than r/n < 1/(2m) fraction of errors
 * *****************************************************/

const int FLD_DIM = 10;
const int CODE_N = 500;
const int CODE_R = 40;
const int CODE_T = CODE_R/2; // guaranteed to correct up to CODE_T errors
const int MASK = (1 << FLD_DIM) - 1;
const int REDUCER = (1<<3) + 1; // from minimal polynomial
typedef uint16_t GF; // make larger than 16-bits if needed

/*********************************************************************************************
 * Some minimal polynomials
 *
 * first number is degree
 * following numbers are nonzero degrees
 * except for lead (monic) and
 * constant term 1
 *
 * 2,1 : x^2 + x + 1
 * 3,1 : x^3 + x + 1
 * 4,1
 * 5,2
 * 6,1
 * 7,1
 * 8,4,3,1 : x^8 + x^4 + x^3 + x + 1
 * 9,1
 * 10,3
 * 11,2
 * 12,3
 * 13,4,3,1
 * 14,5
 * 15,1
 *
 * upgrade GF typedef to uint32_t for bigger fields
 *
 * 16,5,3,1
 * 17,3
 * 18,3
 * 19,5,2,1
 * 20,3
 * 21,2
 * 22,1
 * 23,5
 * 24,4,3,1
 * 25,3
 * 26,4,3,1
 * 27,5,2,1
 * 28,1
 * 29,2
 * 30,1
 * *******************************************************************************************/

/*********************************************************************************************
 *
 * function prototypes
 *
 * *******************************************************************************************/
// starting at line 125
GF gf_add(GF a, GF b);
GF gf_mult(GF a, GF b);
GF gf_div(GF a, GF b, GF lut[]);
GF gf_fermat_inv(GF a);
void create_inverse_table(GF* inverse_table);
GF gf_inv(GF a, GF lut[]);
GF gf_pow(GF a, int e);
// starting at line 200
void print_poly(GF f[], int d);
void poly_add(GF a[], GF b[], GF c[], int d);
void poly_scale(GF f[], GF s, GF f_scaled[], int d);
void poly_mult(GF a[], GF b[], GF c[], int d);
GF horner_eval(GF x, GF f[], int deg);
void poly_div_rem(GF a[], GF b[], GF quot[], GF rem[], int n, GF inverse_table[]);
void EEA(GF f[], GF g[], GF u[], GF v[], GF r[], int n, int thresh, GF inverse_table[]);
// starting at line 375
bool binary_rref(int nrows, int ncols, uint8_t M[nrows][ncols]);
void restrict_pcheck(GF H[CODE_R][CODE_N], uint8_t Hprime[FLD_DIM * CODE_R][CODE_N]);
void binary_to_gf(uint8_t v[], int n, GF vprime[]);
// starting at line 460
void random_permutation(int n, int L[]);
void permute_units(GF units[]);
GF rand_gf();
void random_binary_error(int t, int n, uint8_t e[]);
bool random_alternant_code(uint8_t Hprime[CODE_R*FLD_DIM][CODE_N], GF h_list[], GF alpha_list[]);
// starting at line 555
void syndrome(uint8_t e[CODE_N], uint8_t Hprime[CODE_R*FLD_DIM][CODE_N], uint8_t S[CODE_R*FLD_DIM]);
void BM_syndrome_decode(uint8_t S[], uint8_t e[], GF alpha_list[], GF inverse_table[]);
void euclid_syndrome_decode(uint8_t S[], uint8_t e[], GF alpha_list[], GF inverse_table[]);
void actual_sigma(uint8_t e[], int t, GF sigma[], GF alpha_list[]);
// main() at line 730

/*********************************************************************************************
 *
 * Field arithmetic
 *
 * field elements GF are uint16_t
 *
 * *******************************************************************************************/
// Field addition
GF gf_add(GF a, GF b) {
    return a ^ b;
}

// Field multiplication, most significant exponent
GF gf_mult(GF a, GF b) {
    GF s = 0;
    int eps = 0;
    for (int i = FLD_DIM - 1; i > -1; i--) {
        eps = (s >> (FLD_DIM - 1)) & 1;
        s = ((s << 1) & MASK) ^ (eps * REDUCER) ^ (((b >> i) & 1) * a);
    }
    return s & MASK;
}

// Field division, assumes denominator is non-zero
GF gf_div(GF a, GF b, GF lut[]) {
    return gf_mult(a, gf_inv(b, lut));
}

/******************************************
 * Use Fermat inverse to create lookup table
 * for field inverses.
 *******************************************/
GF gf_fermat_inv(GF a) {
    return gf_pow(a, (1 << FLD_DIM) - 2);
}

void create_inverse_table(GF* inverse_table) {
    for (int i = 0; i < (1 << FLD_DIM); i++) {
        inverse_table[i] = gf_fermat_inv(i);
    }
}

// Field inverses using LUT
// from create_inverse_table()
GF gf_inv(GF a, GF lut[]) {
    if (a == 0) {
        //printf("inverting zero!\n");
        return 0;
    } else {
        return lut[a];
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

/*********************************************************************************************
 *
 * polynomial arithmetic (GF coefficients)
 *
 * *******************************************************************************************/
void print_poly(GF f[], int d){
    for (int i = 0; i <= d; i++){
        printf("%i ", f[i]);
    }
    printf("\n");
}

// assumes a and b are of degree d
void poly_add(GF a[], GF b[], GF c[], int d){
    for (int i = 0; i <= d; i++){
        c[i] = gf_add(a[i], b[i]);
    }
}

void poly_scale(GF f[], GF s, GF f_scaled[], int d){
    for (int i = 0; i <= d; i++){
        f_scaled[i] = gf_mult(f[i], s);
    }
}

// multiplication in GF(2^FLD_DIM)[x]
// most significant exponent
// a, b, c degree d (c is truncated from c_large)
void poly_mult(GF a[], GF b[], GF c[], int d) {
    GF c_large[2*d+1];
    for (int i = 0; i < 2*d + 1; i++) {
        c_large[i] = 0;
    }
    for (int i = d; i > -1; i--) {
        // multiply by x, i.e. shift
        for (int j = 2*d; j > 0; j--) {
            c_large[j] = c_large[j-1];
        }
        c_large[0] = 0;
        // add b[i]*a to c
        for (int j = 0; j < d + 1; j++) {
            c_large[j] = gf_add(c_large[j], gf_mult(b[i], a[j]));
        }
    }
    for (int i = 0; i < d + 1; i++) {
        c[i] = c_large[i];
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

int get_degree(GF f[], int d){
    for (int i = d; i >= 0; i--){
        if (f[i] != 0){
            return i;
        }
    }
    return -1;
}

// n is length of arrays, max_degree - 1
// seems to work
void poly_div_rem(GF a[], GF b[], GF quot[], GF rem[], int n, GF inverse_table[]){

    GF monomial[n];
    GF b_scaled[n];
    GF denom;
    int shift_deg;
    for (int i = 0; i < n; i++){
        quot[i] = 0;
        rem[i] = 0;
        rem[i] = a[i];
    }
    GF scalar;
    int d_a = get_degree(a, n-1);
    int d_b = get_degree(b, n-1);

    if (d_b > -1){
        denom = gf_inv(b[d_b], inverse_table);
    }
    for (int i = d_a; i >= d_b; i--){
        scalar = gf_mult(rem[i], denom);
        if (scalar != 0){
            poly_scale(b, scalar, b_scaled, n);
            shift_deg = i - d_b;
            for (int j = 0; j < n; j++){
                monomial[j] = 0;
            }
            monomial[shift_deg] = scalar;
            poly_add(quot, monomial, quot, n - 1);
            // shift b_scaled
            for (int j = n; j >= shift_deg ; j--){
                b_scaled[j] = b_scaled[j - shift_deg];
            }
            for (int j = 0; j < shift_deg; j++){
                b_scaled[j] = 0;
            }
            poly_add(rem, b_scaled, rem, n - 1);
        }
    }
}

// garbage, but maybe works
void EEA(GF f[], GF g[], GF u[], GF v[], GF r[], int n, int thresh, GF inverse_table[]){
    GF r1[n];
    GF r0[n];
    GF u1[n];
    GF u0[n];
    GF v1[n];
    GF v0[n];

    GF quot[n];
    GF rem[n];

    GF temp[n];
    int check_deg_0, check_deg_1;

    for (int i = 0; i < n; i++){
        r0[i] = f[i];
        r1[i] = g[i];
        u0[i] = 0;
        v0[i] = 0;
        u1[i] = 0;
        v1[i] = 0;
    }
    u0[0] = 1;
    v1[0] = 1;

    check_deg_0 = get_degree(r0, n - 1);
    check_deg_1 = get_degree(r1, n - 1);

    while (!(check_deg_0 >= thresh && check_deg_1 < thresh)){
        poly_div_rem(r0, r1, quot, rem, n, inverse_table);
        // updates
        poly_mult(quot, r1, temp, n-1);
        poly_add(temp, r0, temp, n-1);
        for (int i = 0; i < n; i++){
            r0[i] = r1[i];
            r1[i] = temp[i];
        }
        poly_mult(quot, u1, temp, n-1);
        poly_add(temp, u0, temp, n-1);
        for (int i = 0; i < n; i++){
            u0[i] = u1[i];
            u1[i] = temp[i];
        }
        poly_mult(quot, v1, temp, n-1);
        poly_add(temp, v0, temp, n-1);
        for (int i = 0; i < n; i++){
            v0[i] = v1[i];
            v1[i] = temp[i];
        }
        check_deg_0 = get_degree(r0, n - 1);
        check_deg_1 = get_degree(r1, n - 1);
    }
    for (int i = 0; i < n; i++){
        v[i] = v1[i];
        u[i] = u1[i];
        r[i] = r1[i];
    }
}


/*********************************************************************************************
 *
 * linear algebra
 *
 * *******************************************************************************************/

// Checks if Hprime has full rank on its left side
// i.e. could be reduced to (I|K)
bool binary_rref(int nrows, int ncols, uint8_t Hprime[nrows][ncols]) {
    // to flag no pivot in column
    bool zero_pivot;
    // for temporary rows
    uint8_t tempi[ncols];
    uint8_t tempj[ncols];
    uint8_t M[nrows][ncols];

    // copy Hprime to M, want to keep Hprime
    for (int i = 0; i < nrows; i++){
        for (int j = 0; j < ncols; j++){
            M[i][j] = Hprime[i][j];
        }
    }

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
                //printf("no non-zero pivot in column %i\n", i);
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
    }
    // M now in systematic upper triangular
    return true;
}

// Restricts coefficients of parity-check over GF(2^m) to a binary matrix,
// replacing entries with columns big-endian to little-endian
// GF 100101 -> column (1,0,1,0,0,1)^t
void restrict_pcheck(GF H[CODE_R][CODE_N], uint8_t Hprime[FLD_DIM * CODE_R][CODE_N]) {
    GF temp;
    for (int i = 0; i < CODE_R; i++) {
        for (int j = 0; j < CODE_N; j++) {
            temp = H[i][j];
            for (int k = 0; k < FLD_DIM; k++){
                Hprime[FLD_DIM * i + k][j] = temp & 1;
                temp >>= 1;
            }
        }
    }
}

// grabs bits FLD_DIM at a time to form GF elements
void binary_to_gf(uint8_t v[], int n, GF vprime[]){
    GF temp;
    for (int i = 0; i < n; i += FLD_DIM){
        temp = 0;
        for (int j = 0; j < FLD_DIM; j++){
            temp += (1<<j)*v[i+j];
        }
        vprime[i/FLD_DIM] = temp;
    }
}

/*********************************************************************************************
 *
 * randomized stuff
 *
 * *******************************************************************************************/
// permutes the integer array L
void random_permutation(int n, int L[]){
    int j;
    int temp;
    for (int i = n - 1; i > 0; i--) {
        j = rand() % (i + 1);
        temp = L[j];
        L[j] = L[i];
        L[i] = temp;
    }
}

void permute_units(GF units[]){
    int j;
    GF temp;
    for (int i = (1<<FLD_DIM) - 2; i > 0; i--) {
        j = rand() % (i + 1);
        temp = units[j];
        units[j] = units[i];
        units[i] = temp;
    }
}

GF rand_gf(){
    return rand() % (1 << FLD_DIM);
}

// choose t coordinates of e[] to be 1,
// others to 0
void random_binary_error(int t, int n, uint8_t e[]){
    int first_n[n];
    for (int i = 0; i < n; i++){
        first_n[i] = i;
        e[i] = 0;
    }
    random_permutation(n, first_n);
    for (int i = 0; i < t; i++){
        e[first_n[i]] = 1;
    }
}

// fills Hprime with binary parity check matrix of [CODE_N, CODE_R*FLD_DIM]
// bits are stores as uint8_t *shrug*
bool random_alternant_code(uint8_t Hprime[CODE_R*FLD_DIM][CODE_N], GF h_list[], GF alpha_list[]){
    int nonzero = 0;
    GF alpha = 1;
    bool OK = false;
    bool success = true;
    int fail_limit = 20; // in case of bad parameters making random code unlikely
    GF H[CODE_R][CODE_N]; // canonical parity-check
    // entire field to permute
    GF units[(1 << FLD_DIM) - 1];
    for (int i = 0; i < (1 << FLD_DIM) - 1; i++){
        units[i] = i+1;
    }
    // rejection sampling
    int counter = 0;
    while (!OK && (counter < fail_limit)){
        counter++;
        // create h_list
        nonzero = 0;
        while (nonzero < CODE_N){
            h_list[nonzero] = rand_gf();
            if (h_list[nonzero] != 0){
                nonzero++;
            }
        }
        // create alpha_list
        permute_units(units);
        for (int i = 0; i < CODE_N; i++){
            alpha_list[i] = units[i];
        }
        // create canonical parity-check
        for (int j = 0; j < CODE_N; j++){
            alpha = 1; // holds powers
            for (int i = 0; i < CODE_R; i++){
                H[i][j] = gf_mult(h_list[j], alpha);
                alpha = gf_mult(alpha, alpha_list[j]);
            }
        }
        // restrict to binary
        restrict_pcheck(H, Hprime);
        // check rref for control over code parameters
        OK = binary_rref(CODE_R*FLD_DIM, CODE_N, Hprime);
    }
    if (counter == fail_limit){
        success = false;
    }
    return success;
}

/*********************************************************************************************
 *
 * encode/decode
 *
 * *******************************************************************************************/

// matrix multiplication of error (or received message) and Hprime
void syndrome(uint8_t e[CODE_N], uint8_t Hprime[CODE_R*FLD_DIM][CODE_N], uint8_t S[CODE_R*FLD_DIM]){
    for (int i = 0; i < CODE_R*FLD_DIM; i++){
        S[i] = 0;
        for (int j = 0; j < CODE_N; j++){
            S[i] ^= Hprime[i][j] & e[j];
        }
    }
}

// wikipedia version from Massey's paper pg. 124
// http://crypto.stanford.edu/~mironov/cs359/massey.pdf
// perhaps compare with the constant time version in McEliece stuff
void BM_syndrome_decode(uint8_t S[], uint8_t e[], GF alpha_list[], GF inverse_table[]){

    // syndrome from GF(2) to GF(2^FLD_DIM)
    GF S_gf[CODE_R];
    binary_to_gf(S, CODE_R*FLD_DIM, S_gf);
    GF temp[CODE_T + 1]; // holds copies of whatever
    GF T[CODE_T + 1]; // holds copies of whatever

    /*
    // for "polar form" error locator
    for (int i = 0; i < CODE_R; i++){
        temp[i] = S_gf[i];
    }
    for (int i = 0; i < CODE_R; i++){
        S_gf[i] = temp[CODE_R - 1 - i];
    }
    */

    int L = 0;
    int m = 1;
    GF b = 1;
    GF d; // discrepancy
    GF scalar; // holds d/b

    GF C[CODE_T + 1]; // degree CODE_T possible errors
    for (int i = 0; i < CODE_T + 1; i++) {
        C[i] = 0;
    }
    C[0] = 1;

    GF B[CODE_T + 1];
    for (int i = 0; i < CODE_T + 1; i++) {
        B[i] = 0;
    }
    B[0] = 1;

    for (int n = 0; n < CODE_R; n++) {
        // discrepancy
        d = 0;
        for (int i = 0; i <  L + 1; i++) {
            d = gf_add(d, gf_mult(C[i], S_gf[n - i]));
        }
        if (d == 0) {
            m = m + 1;
        } else if (2 * L > n) {
            // copy of B to shift
            for (int i = 0; i < CODE_T + 1; i++){
                temp[i] = B[i];
            }
            // update C = C - (d/b) * x^m * B
            scalar = gf_mult(d, gf_inv(b, inverse_table));
            // temp = x^m*B
            for (int i = CODE_T; i >= m; i--) {
                temp[i] = temp[i-m];
            }
            for (int i = 0; i < m; i++){
                temp[i] = 0;
            }
            poly_scale(temp, scalar, temp, CODE_T);
            poly_add(C, temp, C, CODE_T);

            m = m + 1;
        } else {
            // copy of C to replace B later
            for (int i = 0; i < CODE_T + 1; i++){
                T[i] = C[i];
            }

            // copy of B to shift
            for (int i = 0; i < CODE_T + 1; i++){
                temp[i] = B[i];
            }
            // update C = C - (d/b) * x^m * B
            scalar = gf_mult(d, gf_inv(b, inverse_table));
            // temp = x^m*B
            for (int i = CODE_T; i >= m; i--) {
                temp[i] = temp[i-m];
            }
            for (int i = 0; i < m; i++){
                temp[i] = 0;
            }
            poly_scale(temp, scalar, temp, CODE_T);
            poly_add(C, temp, C, CODE_T);
            L = n + 1 - L;

            // update B = old C
            for (int i = 0; i < CODE_T + 1; i++){
                B[i] = T[i];
            }
            b = d;
            m = 1;
        }
    }

    GF value;
    for (int i = 0; i < CODE_N; i++) {
        value = horner_eval(gf_inv(alpha_list[i], inverse_table), C, CODE_T);
        if (value == 0) {
            e[i] = 1;
        } else {
            e[i] = 0;
        }
    }
    /*
    // for polar form
    GF value;
    for (int i = 0; i < CODE_N; i++) {
        value = horner_eval(alpha_list[i], C, CODE_T);
        if (value == 0) {
            e[i] = 1;
        } else {
            e[i] = 0;
        }
    }
    */
}

void euclid_syndrome_decode(uint8_t S[], uint8_t e[], GF alpha_list[], GF inverse_table[]){
    // syndrome from GF(2) to GF(2^FLD_DIM)
    GF S_gf[CODE_R];
    binary_to_gf(S, CODE_R*FLD_DIM, S_gf);

    // working with common degree upper bound of common - 1
    // so extend syndrome
    int common = 2*CODE_R+2;
    GF monomial[common]; // holds x^r
    GF S_ext[common]; // copy of S_gf with zero appended
    for (int i = 0; i < common; i++){
        monomial[i] = 0;
        if (i < CODE_R){
            S_ext[i] = S_gf[i];
        } else {
            S_ext[i] = 0;
        }
    }
    monomial[CODE_R] = 1;

    GF sigma[common];
    GF omega[common];
    GF u[common];

    EEA(monomial, S_ext, u, sigma, omega, common, CODE_T, inverse_table); // might need CODE_R even...
    GF value;
    for (int i = 0; i < CODE_N; i++) {
        value = horner_eval(gf_inv(alpha_list[i], inverse_table), sigma, common - 1);
        if (value == 0) {
            e[i] = 1;
        } else {
            e[i] = 0;
        }
    }
}

// for testing, computes actual error locator (inverse roots version)
void actual_sigma(uint8_t e[], int t, GF sigma[], GF alpha_list[]){
    GF temp[t + 1];
    for (int i = 0; i <= t; i++){
        sigma[i] = 0;
    }
    sigma[0] = 1;
    for (int i = 0; i < CODE_N; i++){
        if (e[i] == 1){
            for (int j = 0; j < t+1; j++){
                temp[j] = sigma[j];
            }
            poly_scale(sigma, alpha_list[i], sigma, t);
            for (int j = t; j > 0 ; j--){
                sigma[j] = sigma[j-1];
            }
            sigma[0] = 0;
            poly_add(temp, sigma, sigma, t);
        }
    }
}
/*********************************************************************************************
 *
 * main(), tests out the above
 *
 * *******************************************************************************************/

int main(){

    srand(time(NULL));

    uint16_t inverse_table[1<<FLD_DIM];
    create_inverse_table(inverse_table);

    int num_codes = 5;
    int num_trials = 5;

    uint8_t Hprime[CODE_R*FLD_DIM][CODE_N]; // stores parity-check matrix
    uint8_t S[CODE_R*FLD_DIM]; // syndrome over GF(2)
    uint8_t e[CODE_N]; // for random errors
    uint8_t e_found_BM[CODE_N]; // for detected errors
    uint8_t e_found_euclid[CODE_N]; // for detected errors

    int t; // for number of errors
    bool success_BM; // if e = e_found_BM
    bool success_euclid; // if e = e_found_euclid

    GF h_list[CODE_N]; // holds code parameters
    GF alpha_list[CODE_N]; // holds code parameters

    bool code_found; // prevent infinite or long loops in random code generation

    GF sigma[CODE_T + 1]; // for testing actual error locator, degree at most CODE_T

    printf("**************************************************************************\n");
    printf("* Creating %i random alternant codes GF(2^%i)/GF(2) with\n", num_codes, FLD_DIM);
    printf("* parameters [n, k, d] = [%i, %i, >= %i],\n",CODE_N, CODE_N - CODE_R*FLD_DIM, CODE_R + 1);
    printf("* i.e. a %i-dimensional subspace of an %i-dimensional subspace over GF(2)\n", CODE_N - CODE_R*FLD_DIM, CODE_N);
    printf("* with minimum weight codewords of weight d >= %i.\n", CODE_R+1);
    printf("*\n");
    printf("* Then correcting %i random error vectors with each\n", num_trials);
    printf("* using two decoding algorithms (Berlekamp-Massy and Euclid).\n");
    printf("*\n");
    printf("* The code can correct at least %i errors, i.e. %.2f%%.\n", CODE_T, (double)100*CODE_T/CODE_N);
    printf("* The number of errors to be corrected is randomly chosen from [0, %i].\n", CODE_T);
    printf("*\n");
    for (int i = 0; i < num_codes; i++){
        code_found = random_alternant_code(Hprime, h_list, alpha_list); // fills Hprime, h_list, and alpha_list
        if (code_found){
            printf("*\tRandom code %i constructed successfully.\n", i);
            /*
            for (int k = 0; k < CODE_N; k++){
                if (alpha_list[k] == 0){
                    printf("*\t(Zero in support)\n");
                }
            }
            */
            for (int j = 0; j < num_trials; j++){
                // generate random error and compute syndrome
                t = rand() % (CODE_T + 1); // random number of errors below max
                random_binary_error(t, CODE_N, e); // random error e
                /*
                for (int k = 0; k < CODE_N; k++){
                    if (e[k] == 1 && alpha_list[k] == 0){
                        printf("*\t\t(Error at alpha = 0)\n");
                    }
                }
                */
                syndrome(e, Hprime, S); // S is syndrome of e via Hprime
                // decode
                BM_syndrome_decode(S, e_found_BM, alpha_list, inverse_table);
                euclid_syndrome_decode(S, e_found_euclid, alpha_list, inverse_table);
                success_BM = true;
                success_euclid = true;
                for (int k = 0; k < CODE_N; k++){
                    if (e[k] != e_found_BM[k]){
                        success_BM = false;
                        //printf("*\t\t\tFailure at index %i h = %i alpha = %i : (%i, %i)\n", k, h_list[k], alpha_list[k], e[k], e_found_BM[k]);
                    }
                    if (e[k] != e_found_euclid[k]){
                        success_euclid = false;
                        //printf("*\t\t\tFailure at index %i h = %i alpha = %i : (%i, %i)\n", k, h_list[k], alpha_list[k], e[k], e_found_euclid[k]);
                    }
                }
                // report success or failure
                if (success_BM){
                    printf("*\t\tSuccess BM\t(corrected %i errors)\n", t);
                } else {
                    printf("*\t\tFailure BM\t(%i actual errors)\n", t);
                }
                // report success or failure
                if (success_euclid){
                    printf("*\t\tSuccess Euclid\t(corrected %i errors)\n", t);
                } else {
                    printf("*\t\tFailure Euclid\t(%i actual errors)\n", t);
                }
                /*
                // print the actual error locator
                actual_sigma(e, CODE_T, sigma, alpha_list);
                printf("*\t\tsigma_actual:\t");
                for (int i = 0; i < CODE_T + 1; i++){
                    printf("%i ", sigma[i]);
                }
                printf("\n*\n");
                */
            }
        } else {
            printf("*\tFailed to generate random code %i; check parameters.\n", i);
        }
    }
    printf("**************************************************************************\n");

    return 0;
}
