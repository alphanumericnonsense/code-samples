#ifndef NODOUBLES
#define NODOUBLES
#include <stdint.h>
#include <stdbool.h>

// Field elements stored as integers, binary coefficients
// are coefficients modulo the minimal polynomial
typedef uint16_t GF;

/******************************
* Field constants definitions
******************************/

//MIN_POLY = "1000000001001" degree high downto low
//EXTRA_MIN_POLY = [(64,1),(3,1),(1,1),(0,2)]

extern const int GOPPA_T;
extern const int CODE_N;

extern const int GF_M;
extern const GF MASK;
extern const GF REDUCER;
extern GF EXTRA_MIN_POLY[65];

// compiler complained about 1<<GF_M in array declaration...
extern const GF INVERSE_TABLE[1 << 12];

// valid secret key for testing, other testing stuff...
extern GF G_TEST[65]; // Goppa polynomial
extern GF L_TEST[3488]; //  support
extern uint8_t PLAIN_TEST[3488];
extern uint8_t CIPHER_TEST[768];
extern GF S_TEST[128]; // double syndrome unrestricted
extern uint8_t C2_TEST[1536]; // double syndrome
extern GF SIGMA_TEST[65]; // error locator

/*********************************************
 * Prototypes
 * for all functions from each *.c
 *********************************************/

// GF.c
GF gf_pow(GF a, int e);
GF gf_add(GF a, GF b);
GF gf_mult(GF a, GF b);
GF gf_div(GF a, GF b);
GF gf_fermat_inv(GF a);
GF gf_inv(GF a);
void create_inverse_table(GF* inverse_table);

// pk_gen.c
GF horner_eval(GF x, GF f[], int deg);
bool pk_gen(GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T]);
void canonical_pcheck(GF g[GOPPA_T + 1], GF support[CODE_N], GF H[GOPPA_T][CODE_N]);
void restrict_pcheck(GF H[GOPPA_T][CODE_N], uint8_t K[GF_M * GOPPA_T][CODE_N]);
bool binary_rref(int nrows, int ncols, uint8_t M[nrows][ncols]);

// sk_gen.c
void permute_field(GF L[1 << GF_M]);
void sk_gen(GF g[GOPPA_T + 1], GF support[CODE_N]);
void poly_mod_prod(GF a[], int d_a, GF b[], int d_b, GF f[], int d_f, GF result[]);
void poly_mult(GF a[], int d_a, GF b[], int d_b, GF c[]);
bool gf_rref(int nrows, int ncols, GF M[nrows][ncols]);

// encrypt.c
void encrypt(uint8_t plain[CODE_N], uint8_t K[GF_M*GOPPA_T][CODE_N - GF_M*GOPPA_T], uint8_t cipher[GF_M*GOPPA_T]);
void random_message(uint8_t plain[CODE_N]);
void random_permutation(int N, uint16_t L[N]);

// decrypt.c
void decrypt(uint8_t cipher[GF_M*GOPPA_T], GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t plain[CODE_N]);
void canonical_pcheck2(GF g2[2*GOPPA_T + 1], GF support[CODE_N], GF H2[2*GOPPA_T][CODE_N]);
void g_to_g2(GF g[GOPPA_T + 1], GF g2[2*GOPPA_T + 1]);
void restrict_pcheck2(GF H2[2*GOPPA_T][CODE_N], uint8_t H2prime[2*GF_M * GOPPA_T][CODE_N]);
void decode(GF S[2*GOPPA_T], GF sigma[GOPPA_T + 1]);
#endif
