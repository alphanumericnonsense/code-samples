/* ****************************************************************************************
 *
 * A simple implementation of Saber.
 *
 * Doesn't use keccak for randomness or pack data as in the specs.
 *
 * Uses karatsuba multiplication for polynomials (deg 255),
 * not fully recursive, switches over to schoolbook at degree 7.
 *
 * main() generates keys and random messages, encrypts/decrypts and checks for success.
 *
 * One nice thing about Saber is that modular arithmetic is modulo powers of 2,
 * so overflow is fine and we can use masking x & ((1<<d) - 1).
 *
 * ****************************************************************************************/

#include <stdio.h>
#include <stdint.h>
//#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

const int N = 256; // ring dimension over base ring
const int L = 4; // module dimension, 2/3/4

const unsigned int E_P = 10; // fixed
const unsigned int E_Q = 13; // fixed
const unsigned int E_T = 6; // 3/4/6
const unsigned int P = 1<<10; // fixed
const unsigned int Q = 1<<13; // fixed
const unsigned int T = 1<<6; // 8/16/64
const int MU = 6; // CBD parameter, 10/8/6

typedef int16_t poly[256];
typedef poly vec[4]; // L = 2/3/4

// constants for deterministic error
const int16_t h1 = 1 << (E_Q - E_P - 1);
const int16_t h2 = ( 1<< (E_P - 2) ) - (1 << (E_P - E_T - 1)) + (1 << (E_Q - E_P - 1));
// some unnecessary constants built from h1, h2
poly H1 = {h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1};
// const vec H = {H1, H1, H1, H1}; // L = 2/3/4 times
poly H2 = {h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2,
                h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2, h2};

/* ***********************************
 *
 * testing stuff
 *
 * ***********************************/

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
    for (int i = 0; i < L; i++){
        print_poly(v[i]);
    }
}

/**************************************************************
 *
 * rounding (right shifts, reductions mod 2^d)
 *
 * ************************************************************/

// coefficients modulo 2^d
void poly_mod(poly f, unsigned int d){
    for (int i = 0; i < N; i++){
        f[i] = f[i] & ((1<<d) - 1);
    }
}

// coefficients modulo 2^d
void vec_mod(vec v, unsigned int d){
    for (int i = 0; i < L; i++){
        poly_mod(v[i], d);
    }
}

// coefficients d most significant bits
void poly_rshift(poly f, unsigned int d){
    for (int i = 0; i < N; i++){
        f[i] = f[i] >> d;
    }
}

// coefficients d most significant bits
void vec_rshift(vec v, unsigned int d){
    for (int i = 0; i < L; i++){
        poly_rshift(v[i], d);
    }
}

/**************************************************************
 *
 * randomness
 *
 * ************************************************************/

// uniform mod 2^d
int16_t uniform(unsigned int d){
    return rand() % (1<<d);
}

//uniform mod 2^d
void poly_uniform(poly f, unsigned int d){
    for (int i = 0; i < N; i++){
        f[i] = uniform(d);
    }
}

int flip(){
    return rand() % 2;
}

// centered binomial distribution supported in [-mu/2, mu/2]
int CBD(int mu){
    int16_t total = 0;
    for (int i = 0; i < mu/2; i++){
        total += flip();
        total -= flip();
    }
    return total;
}

void poly_CBD(poly f, int mu){
    for (int i = 0; i < N; i++){
        f[i] = CBD(mu);
    }
}

void vec_CBD(vec v, int mu){
    for (int i = 0; i < L; i++){
        poly_CBD(v[i], mu);
    }
}

// random uniform 256-bit message
void random_message(poly m){
    poly_uniform(m, 1);
}

/*************************************************************
 *
 * algebra
 *
 * ***********************************************************/

// add mod 2^d
void poly_add(poly a, poly b, poly c, unsigned int d){
    for (int i = 0; i < N; i++){
        c[i] = (a[i] + b[i]);
    }
    poly_mod(c, d);
}

// subtract mod 2^d
void poly_sub(poly a, poly b, poly c, unsigned int d){
    for (int i = 0; i < N; i++){
        c[i] = (a[i] - b[i]);
    }
    poly_mod(c, d);
}

void vec_add(vec a, vec b, vec c, unsigned int d){
    for (int i = 0; i < L; i++){
        poly_add(a[i], b[i], c[i], d);
    }
}

void vec_sub(vec a, vec b, vec c, unsigned int d){
    for (int i = 0; i < L; i++){
        poly_sub(a[i], b[i], c[i], d);
    }
}


// schoolbook polynomial multiplication,
// used as base in KM()
void SM(int16_t a[], int16_t b[], int16_t c[], int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            c[i + j] = a[i]*b[j];
        }
    }
}

// karatsuba for integer polynomials
// assumes n = degree + 1 is a power of 2
// runs fastest with switch to schoolbook around n = 8 or 16
void KM(int16_t a[], int16_t b[], int16_t c[], int n){
    for (int i = 0; i < 2*n; i++){
        c[i] = 0;
    }
    // can change base case between n = 256 and n = 1
    // n = 8 or 16 best
    if (n == 8){
        SM(a, b, c, n);
        //c[0] = a[0]*b[0];
    } else {
        int16_t a0[n/2], a1[n/2], b0[n/2], b1[n/2], a0_a1[n/2], b0_b1[n/2], x[n], y[n], z[n];
        for (int i = 0; i < n/2; i++){
            a0[i] = a[i];
            a1[i] = a[i + n/2];
            b0[i] = b[i];
            b1[i] = b[i + n/2];
            a0_a1[i] = a0[i] + a1[i];
            b0_b1[i] = b0[i] + b1[i];
        }
        KM(a1, b1, x, n/2);
        KM(a0, b0, y, n/2);
        KM(a0_a1, b0_b1, z, n/2);
        for (int i = 0; i < n; i++){
            c[i] += y[i];
            c[i + n/2] += (z[i] - x[i] - y[i]);
            c[i + n] += x[i];
        }
    }
}

void poly_mult(poly a, poly b, poly c, unsigned int d){
    // initialize output and temp to zero
    int16_t temp[2*N];
    for (int i = 0; i < N; i++){
        c[i] = 0;
        temp[i] = 0;
        temp[i + N] = 0;

    }
    // schoolbook for testing
    // seems faster than karatsuba?
    /*
    int16_t eps, index;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            index = i + j;
            if (index >= N) {
                eps = -1;
            } else {
                eps = 1;
            }
            index %= N;
            c[index] += a[i]*b[j]*eps;
            //c[index] &= ((1<<d) - 1);
        }
    }
    poly_mod(c, d);
    */
    KM(a, b, temp, N);
    // fill c, include negacyclic shift of deg 256+ from temp
    for (int i = 0; i < N; i++){
        c[i] = temp[i] - temp[i + N];
    }
    // reduction mod 2^d
    poly_mod(c, d);

}

void dot_prod(vec u, vec v, poly p, unsigned int d){
    // initialize p = 0
    for (int i = 0; i < N; i++){
        p[i] = 0;
    }
    // accumulate
    poly temp;
    for (int i = 0; i < L; i++){
        poly_mult(u[i], v[i], temp, d);
        poly_add(temp, p, p, d);
    }
}

// scale by 2^d
void poly_scale(poly f, unsigned int d){
    for (int i = 0; i < N; i++){
        f[i] = f[i] << d;
    }
}

/*************************************************************
 *
 * key generation, encryption, decryption
 *
 * ***********************************************************/

// generate uniform A, CBD s, and and b = round(A^ts+h)
void key_gen(poly A[L][L], vec b, vec s){
    // fill A
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            poly_uniform(A[i][j], E_Q);
        }
    }
    // fill s
    vec_CBD(s, MU);
    // fill b
    vec col;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < N; k++){
                col[j][k] = A[j][i][k];
            }
        }
        dot_prod(col, s, b[i], E_Q);
    }
    // build a constant I couldn't declare globally
    // and don't really need
    vec H;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < N; j++){
            H[i][j] = h1;
        }
    }
    vec_add(b, H, b, E_Q);
    vec_rshift(b, E_Q - E_P);
}

// encrypt m with (A, b) and store ciphertext (c_m, b')
void encrypt(poly m, poly A[L][L], vec b, poly c_m, vec bprime){
    vec temp;

    // nonce s'
    vec sprime;
    vec_CBD(sprime, MU);

    // fill bprime
    vec row;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            for (int k = 0; k < N; k++){
                row[j][k] = A[i][j][k];
            }
        }
        dot_prod(row, sprime, bprime[i], E_Q);
    }
    // build a constant I couldn't declare globally
    // and don't really need
    vec H;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < N; j++){
            H[i][j] = h1;
        }
    }
    vec_add(bprime, H, bprime, E_Q);
    vec_rshift(bprime, E_Q - E_P);

    // fill c_m
    poly vprime;
    vec_mod(sprime, E_P);
    dot_prod(b, sprime, vprime, E_P);
    poly_scale(m, E_P - 1); // changing m in-place
    poly_sub(H1, m, c_m, E_P);
    poly_add(vprime, c_m, c_m, E_P);
    poly_rshift(c_m, E_P - E_T);
}

// decrypt (c_m, bprime) with s and store in mprime
void decrypt(poly c_m, vec bprime, vec s, poly mprime){
    //printf("in decrypt()\n");
    poly v;
    vec_mod(s, E_P);
    dot_prod(bprime, s, v, E_P);
    poly_scale(c_m, E_P - E_T); // scaling ciphertext in-place
    poly_sub(v, c_m, mprime, E_P);
    poly_add(mprime, H2, mprime, E_P);
    poly_rshift(mprime, E_P - 1);
}

/********************************************************************************************************
 *
 * MAIN: encrypt/decrypt random messages to see if this works
 *
 * ******************************************************************************************************/

int main() {

    // RNG
    srand(time(NULL));
    // timer stuff
    //time_t start, end;
    //double total_time = 0;
    // public key
    poly A[L][L];
    vec b;
    // private key
    vec s;
    // ciphertext
    vec bprime;
    poly c_m;
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
    //start = clock();
    for (int i = 0; i < trials; i++){
        printf("* \tKey pair %i:\n", i);

        key_gen(A, b, s);
        /*
        for (int j = 0; j < L; j++){
            for (int k = 0; k < L; k++){
                printf("A[%i][%i]:\n", j, k);
                print_poly(A[j][k]);
            }
        }
        printf("b:\n");
        print_vec(b);
        printf("s:\n");
        print_vec(s);
        */
        for (int j = 0; j < trials; j++){
            printf("* \t\tMessage %i:  ", j);
            random_message(msent);

            // copy message
            for (int k = 0; k < N; k++){
                mdecrypt[k] = msent[k];
            }
            //printf("msent:\n");
            //print_poly(msent);
            // encrypt/decrypt
            encrypt(mdecrypt, A, b, c_m, bprime);
            //printf("c_m:\n");
            //print_poly(c_m);
            //printf("bprime:\n");
            //print_vec(bprime);
            decrypt(c_m, bprime, s, mdecrypt);

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
    //end = clock();
    //printf("* clock() elapsed: %f\n", (double)(end-start));
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");

    return 0;
}
