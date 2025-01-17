/********************************************************************************************************
 *
 * Simple version of Kyber
 *
 * ******************************************************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

const int Q = 3329; // modulus
const int N = 256; // poly dimension
const int ZETA = 17; // primitive 256th root of unity mod Q

typedef int16_t poly[256]; // elements of Z[x]/(Q, x^256+1)

/* Alter params for security *******************/
const int K = 4; // dimension, 2/3/4
const int d_u = 11; // ciphertext compression, 10/10/11
const int d_v = 5; // ciphertext compression, 4/4/5
const int eta1 = 2; // CBD parameter, 3/2/2
const int eta2 = 2; // CBD parameter, 2/2/2

typedef poly vec[4]; // vectors over poly, K = 2/3/4
/***********************************************/




// powers of ZETA = 17 mod Q, won't let me use N?
const int16_t ZETA_POWERS[257] =
{1, 17, 289, 1584, 296, 1703, 2319, 2804, 1062, 1409, 650, 1063, 1426, 939, 2647, 1722, 2642, 1637, 1197, 375, 3046, 1847, 1438, 1143, 2786, 756, 2865, 2099, 2393, 733, 2474, 2110, 2580, 583, 3253, 2037, 1339, 2789, 807, 403, 193, 3281, 2513, 2773, 535, 2437, 1481, 1874, 1897, 2288, 2277, 2090, 2240, 1461, 1534, 2775, 569, 3015, 1320, 2466, 1974, 268, 1227, 885, 1729, 2761, 331, 2298, 2447, 1651, 1435, 1092, 1919, 2662, 1977, 319, 2094, 2308, 2617, 1212, 630, 723, 2304, 2549, 56, 952, 2868, 2150, 3260, 2156, 33, 561, 2879, 2337, 3110, 2935, 3289, 2649, 1756, 3220, 1476, 1789, 452, 1026, 797, 233, 632, 757, 2882, 2388, 648, 1029, 848, 1100, 2055, 1645, 1333, 2687, 2402, 886, 1746, 3050, 1915, 2594, 821, 641, 910, 2154, 3328, 3312, 3040, 1745, 3033, 1626, 1010, 525, 2267, 1920, 2679, 2266, 1903, 2390, 682, 1607, 687, 1692, 2132, 2954, 283, 1482, 1891, 2186, 543, 2573, 464, 1230, 936, 2596, 855, 1219, 749, 2746, 76, 1292, 1990, 540, 2522, 2926, 3136, 48, 816, 556, 2794, 892, 1848, 1455, 1432, 1041, 1052, 1239, 1089, 1868, 1795, 554, 2760, 314, 2009, 863, 1355, 3061, 2102, 2444, 1600, 568, 2998, 1031, 882, 1678, 1894, 2237, 1410, 667, 1352, 3010, 1235, 1021, 712, 2117, 2699, 2606, 1025, 780, 3273, 2377, 461, 1179, 69, 1173, 3296, 2768, 450, 992, 219, 394, 40, 680, 1573, 109, 1853, 1540, 2877, 2303, 2532, 3096, 2697, 2572, 447, 941, 2681, 2300, 2481, 2229, 1274, 1684, 1996, 642, 927, 2443, 1583, 279, 1414, 735, 2508, 2688, 2419, 1175, 1};

// bit reversal of 7-bit integers 0 to 127
const int reversed[128] = {0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120, 4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124, 2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122, 6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126, 1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121, 5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125, 3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123, 7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127};
/********************************************************************************************************
 *
 * Algebra
 *
 * ******************************************************************************************************/

/*******************************************
 * Barrett reduction, 1/Q approx. 5039/2^24
 * Intermediate steps fit into 32 bits
 * Verified for |a| < (Q-1)^2,
 * i.e. handles reducing multiplication.
 * *****************************************/
int16_t barrett_reduce(int32_t a){
    int16_t b = a - ((int32_t)(5039*(a>>10))>>14)*Q;
    if (b >= Q){
        return b - Q;
    } else if (b >= 0){
        return b;
    } else {
        return b + Q;
    }
}

int16_t mod_Q_mult(int16_t a, int16_t b){
    return barrett_reduce((int32_t) a * (int32_t) b);
}

void poly_add(poly a, poly b, poly c){
    for (int i = 0; i < N; i++){
        c[i] = barrett_reduce(a[i] + b[i]);
    }
}

void poly_sub(poly a, poly b, poly c){
    for (int i = 0; i < N; i++){
        c[i] = barrett_reduce(a[i] - b[i]);
    }
}

void vec_add(vec a, vec b, vec c){
    for (int i = 0; i < K; i++){
        poly_add(a[i], b[i], c[i]);
    }
}

void vec_sub(vec a, vec b, vec c){
    for (int i = 0; i < K; i++){
        poly_sub(a[i], b[i], c[i]);
    }
}

/********************************************************************************************************
 *
 * NTT stuff
 *
 * ******************************************************************************************************/

/****************************************
 * size 128 FFT
 * in-place decimation in frequency
 * input in standard order
 * output in bit reversed order
 * **************************************/
void fft(int16_t f[128]){
    //printf("in ftt\n");
    int16_t u, v;
    int16_t zeta = ZETA;
    int16_t omega = 1;

    for (int m = 128; m > 1; m >>= 1){
        zeta = mod_Q_mult(zeta, zeta);
        omega = 1;
        for (int j = 0; j < m/2; j++){
            for (int r = 0; r < 128; r += m){
                u = f[r + j];
                v = f[r + j + m/2];
                f[r+j] = barrett_reduce(u + v);
                f[r + j + m/2] = mod_Q_mult(barrett_reduce(u - v), omega);
            }
            omega = mod_Q_mult(omega, zeta);
        }
    }
}

/*******************************************
 * size 128 inverse FFT
 * in-place decimation in time
 * input in bit reversed order
 * output in standard order
 * *****************************************/
void inv_fft(int16_t fhat[128]){
    //printf("in inv_ftt\n");
    int16_t a, b;
    int16_t omega = 1;

    for (int m = 1; m < 128; m <<= 1){
        for (int j = 0; j < 128; j += 2*m){
            omega = 1;
            for (int k = 0; k < m; k++){
                //printf("ftt loops %i %i %i\n", m, j, k);
                a = fhat[k + j];
                b = mod_Q_mult(omega, fhat[k + j + m]);
                fhat[k + j] = barrett_reduce(a + b);
                fhat[k + j + m] = barrett_reduce(a - b);
                omega = mod_Q_mult(ZETA_POWERS[256 - 128/m], omega);
            }
        }
    }
    // scale by 1/128, could be done in the loop above
    for (int i = 0; i < 128; i++){
        fhat[i] = mod_Q_mult(3303, fhat[i]);
    }
}

// NTT in place, tested and works
void ntt(poly f){
    int16_t f_even[128];
    int16_t f_odd[128];
    for (int i = 0; i < 128; i++){
        // input scaled
        f_even[i] = mod_Q_mult(ZETA_POWERS[i], f[2*i]);
        f_odd[i] = mod_Q_mult(ZETA_POWERS[i], f[2*i+1]);
    }
    fft(f_even);
    fft(f_odd);
    for (int i = 0; i < 128; i++){
        f[2*i] = f_even[i];
        f[2*i+1] = f_odd[i];
    }
}

// inverse NTT in place, tested and works
void inv_ntt(poly fhat){
    int16_t fhat_even[128];
    int16_t fhat_odd[128];
    for (int i = 0; i < 128; i++){
        fhat_even[i] = fhat[2*i];
        fhat_odd[i] = fhat[2*i+1];
    }
    inv_fft(fhat_even);
    inv_fft(fhat_odd);
    for (int i = 0; i < 128; i++){
        // output scaled
        fhat[2*i] = mod_Q_mult(ZETA_POWERS[N-i], fhat_even[i]);
        fhat[2*i+1] = mod_Q_mult(ZETA_POWERS[N-i], fhat_odd[i]);
    }
}

void vec_ntt(vec v){
    for (int i = 0; i < K; i++){
        //printf("vec_ntt\n");
        ntt(v[i]);
    }
}

void vec_inv_ntt(vec vhat){
    for (int i = 0; i < K; i++){
        //printf("vec_inv_ntt\n");
        inv_ntt(vhat[i]);
    }
}

// bit-reversal of 7-bit i
int bit_reversal(int i){
    int j = 0;
    for (int k = 0; k < 7; k++){
        j += (1<<(7-1-k))*(i&1);
        i >>= 1;
    }
    return j;
}

int iexp(int x, int e){
    int y = 1;
    for (int i = 0; i < e; i++){
        y *= x;
        y %= Q;
    }
    return y;
}

// a*b = c in NTT domain.......
void ntt_domain_mult(poly fhat, poly ghat, poly hhat){

    int16_t x; // temp product
    int16_t y; // temp product

    for (int i = 0; i < N/2; i++){
        x = mod_Q_mult(fhat[2*i], ghat[2*i]);
        y = mod_Q_mult(fhat[2*i+1], ghat[2*i+1]);

        hhat[2*i] = barrett_reduce(x + mod_Q_mult(ZETA_POWERS[2*reversed[i]+1], y));
        hhat[2*i+1] = barrett_reduce(mod_Q_mult(fhat[2*i]+fhat[2*i+1], ghat[2*i]+ghat[2*i+1]) - x - y);
    }
}

void ntt_dot_prod(vec uhat, vec vhat, poly phat){
    // initialize phat
    for (int i =0; i < N; i++){
        phat[i] = 0;
    }
    // multiply components of uhat and vhat, store in temp, add to phat
    poly temp;
    for (int i =0; i < K; i++){
        ntt_domain_mult(uhat[i], vhat[i], temp);
        poly_add(phat, temp, phat);
    }
}

/*********************************************************
 *
 * Compression/decompression
 *
 * *******************************************************/

// [0, Q) -> [0, 2^d)
int16_t mod_Q_compress(int16_t a, int d){
    int num = a*(1<<d);
    int den = Q;
    int rem = num % den;
    int quot = (num - rem)/den;
    if (2*rem < den) {
        return quot % (1<<d);
    } else {
        return (quot + 1) % (1<<d);
    }
}

// [0, 2^d) -> [0, Q)
int16_t mod_Q_decompress(int16_t a, int d){
    int num = a*Q;
    int den = (1<<d);
    int rem = num % den;
    int quot = (num - rem)/den;
    if (2*rem < den) {
        return quot % Q;
    } else {
        return (quot + 1) % Q;
    }
}

void poly_compress(poly a, int d){
    for (int i = 0; i < N; i++){
        a[i] = mod_Q_compress(a[i], d);
    }
}

void poly_decompress(poly a, int d){
    for (int i = 0; i < N; i++){
        a[i] = mod_Q_decompress(a[i], d);
    }
}

void vec_compress(vec a, int d){
    for (int i = 0; i < K; i++){
        poly_compress(a[i], d);
    }
}

void vec_decompress(vec a, int d){
    for (int i = 0; i < K; i++){
        poly_decompress(a[i], d);
    }
}
/********************************************************************************************************
 *
 * Randomness
 *
 * ******************************************************************************************************/
// waste of randomness I suppose
int16_t flip(){
    return rand() % 2;
}

int16_t CBD(int eta){
    int16_t total = 0;
    for (int i = 0; i < eta; i++){
        total += flip();
        total -= flip();
    }
    return total;
}

void poly_CBD(poly a, int eta){
    for (int i = 0; i < N; i++){
        a[i] = CBD(eta);
    }
}

void vec_CBD(vec a, int eta){
    for (int i = 0; i < K; i++){
        poly_CBD(a[i], eta);
    }
}

int16_t mod_Q_uniform(){
    return rand() % Q;
}

void poly_uniform(poly a){
    for (int i = 0; i < N; i++){
        a[i] = mod_Q_uniform();
    }
}

void vec_uniform(vec a){
    for (int i = 0; i < K; i++){
        poly_uniform(a[i]);
    }
}

// fill m with {0, 1}, probability 1/2
void random_message(poly m){
    for (int i =0; i < N; i++){
        m[i] = (flip()+1)/2;
    }
}

/********************************************************************************************************
 *
 * key_gen(), encrypt(), decrypt()
 *
 * ******************************************************************************************************/
// testing prototypes, functions below
void print_poly(poly f);
void print_vec(vec v);

void key_gen(poly Ahat[K][K], vec bhat, vec shat){

    // fill Ahat
    for (int i = 0; i < K; i++){
        for (int j = 0; j < K; j++){
            poly_uniform(Ahat[i][j]);
        }
    }

    // fill s/hat
    vec_CBD(shat, eta1);

    // create and fill error vector e/hat
    vec ehat;
    vec_CBD(ehat, eta1);

    // put s/hat and e/hat into NTT domain
    vec_ntt(shat);
    vec_ntt(ehat);

    // form bhat = Ahat*shat+ehat
    vec row;
    for (int i = 0; i < K; i++){
        // fill row
        for (int j = 0; j < K; j++){
            for (int k = 0; k < N; k++){
                row[j][k] = Ahat[i][j][k];
            }
        }
        // dot row with shat store in bhat[i]
        ntt_dot_prod(row, shat, bhat[i]);
    }
    // add ehat
    vec_add(ehat, bhat, bhat);
}

// encrypt m with (Ahat, bhat) to ciphertext (ucomp, vcomp)
// odd indices of polys in ciphertext are wrong??? even indices are fine???
void encrypt(poly m, poly Ahat[K][K], vec bhat, vec ucomp, poly vcomp){
    // random components
    vec rhat, e1;
    poly e2;
    vec_CBD(rhat, eta1);

    vec_ntt(rhat); // to NTT domain
    vec_CBD(e1, eta2);
    poly_CBD(e2, eta2);

    // form u = A^t*r + e1
    vec col;
    for (int i = 0; i < K; i++){
        // fill col
        for (int j = 0; j < K; j++){
            for (int k = 0; k < N; k++){
                col[j][k] = Ahat[j][i][k];
            }
        }
        // dot col with rhat store in u[i]
        ntt_dot_prod(col, rhat, ucomp[i]);
    }
    // move out of NTT domain
    vec_inv_ntt(ucomp);
    // add e1
    vec_add(e1, ucomp, ucomp);
    // form v = <b,r> + e2 + mdecomp
    ntt_dot_prod(bhat, rhat, vcomp);
    inv_ntt(vcomp);
    poly_add(vcomp, e2, vcomp);
    poly_decompress(m, 1);
    poly_add(m, vcomp, vcomp);

    // compress u, v -> ucomp, vcomp
    vec_compress(ucomp, d_u);
    poly_compress(vcomp, d_v);
}

void decrypt(vec ucomp, poly vcomp, vec shat, poly m){
    // decompress ciphertext
    vec_decompress(ucomp, d_u);
    poly_decompress(vcomp, d_v);
    // form v - <s,u>, storing intermediates in m
    vec_ntt(ucomp);
    ntt_dot_prod(shat, ucomp, m);
    inv_ntt(m);
    poly_sub(vcomp, m, m);
    // compress to output message bits
    poly_compress(m, 1);
}

/********************************************************************************************************
 *
 * Testing bullshit
 *
 * ******************************************************************************************************/
// fancy print binary messages
void print_message(poly f){
    printf("* \n");
    for (int i = 0; i < 4; i++){
        printf("* \t\t");
        for (int j = 0; j < 64; j++){
            printf("%i", f[64*i+j]);
        }
        printf("\n");
    }
    printf("* \n");
}

void print_poly(poly f){
    printf("\n");
    for (int i = 0; i < 16; i++){
        for (int j = 0; j < 16; j++){
            printf("%i, ", f[16*i+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vec(vec v){
    for (int i = 0; i < K; i++){
        print_poly(v[i]);
    }
}

// copy a to b
void poly_copy(poly a, poly b){
    for (int i = 0; i<N; i++){
        b[i] = a[i];
    }
}

/********************************************************************************************************
 *
 * MAIN: encrypt/decrypt random messages to see if this works
 *
 * ******************************************************************************************************/

int main(){

    // RNG
    srand(time(NULL));
    // public key
    poly Ahat[K][K];
    vec bhat;
    // private key
    vec shat;
    // ciphertext
    vec ucomp;
    poly vcomp;
    // message sent/decrypted
    poly msent;
    poly mdecrypt;
    // counts errors in decryption
    int count = 0;
    // number of keys and messages
    int trials = 10;
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("* Encrypting %i random 256-bit messages each from %i key pairs.\n", trials, trials);
    printf("* Recording success or failure (with messages in case of failure).\n");
    printf("* \n");
    for (int i = 0; i<trials; i++){
        printf("* \tKey pair %i:\n", i);

        key_gen(Ahat, bhat, shat);

        for (int j = 0; j<trials; j++){
            printf("* \t\tMessage %i:  ", j);
            random_message(msent);

            //printf("m");
            //print_poly(msent);

            // copy message
            for (int k = 0; k < N; k++){
                mdecrypt[k] = msent[k];
            }

            // encrypt/decrypt
            encrypt(mdecrypt, Ahat, bhat, ucomp, vcomp);
            decrypt(ucomp, vcomp, shat, mdecrypt);

            // check for successful decryption
            count = 0;
            for (int k = 0; k < N; k++){
                if (mdecrypt[k] != msent[k]) {
                    count++;
                }
            }
            if (count == 0){
                printf("Success\n");
            } else {
                printf("Failure %i/256\n", count);
                print_message(msent);
                print_message(mdecrypt);
            }
        }
    }
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");

    return 0;
}
