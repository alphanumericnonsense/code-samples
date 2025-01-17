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

// MIN_POLY =  "11011000000001"degree high down to low
// EXTRA_MIN_POLY = [(128,1),(7,1),(2,1),(1,1),(0,1)], (degree, GF)

extern const int GOPPA_T;
extern const int CODE_N;

extern const int GF_M;
extern const GF MASK;
extern const GF REDUCER;
extern GF EXTRA_MIN_POLY[129];

// compiler complained about 1<<GF_M in array declaration...
extern const GF INVERSE_TABLE[1 << 13];

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
bool pk_gen(GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t (*K)[CODE_N - GF_M*GOPPA_T]);
void canonical_pcheck(GF g[GOPPA_T + 1], GF support[CODE_N], GF (*H)[CODE_N]);
void restrict_pcheck(GF (*H)[CODE_N], uint8_t (*Hprime)[CODE_N]);
bool binary_rref(int nrows, int ncols, uint8_t (*M)[ncols]);

// sk_gen.c
void permute_field(GF L[1 << GF_M]);
void sk_gen(GF g[GOPPA_T + 1], GF support[CODE_N]);
void poly_mod_prod(GF a[], int d_a, GF b[], int d_b, GF f[], int d_f, GF result[]);
void poly_mult(GF a[], int d_a, GF b[], int d_b, GF c[]);
bool gf_rref(int nrows, int ncols, GF M[nrows][ncols]);

// encrypt.c
void encrypt(uint8_t plain[CODE_N], uint8_t (*K)[CODE_N - GF_M*GOPPA_T], uint8_t cipher[GF_M*GOPPA_T]);
void random_message(uint8_t plain[CODE_N]);
void random_permutation(int N, uint16_t L[N]);

// decrypt.c
void decrypt(uint8_t cipher[GF_M*GOPPA_T], GF g[GOPPA_T + 1], GF support[CODE_N], uint8_t plain[CODE_N]);
void canonical_pcheck2(GF g2[2*GOPPA_T + 1], GF support[CODE_N], GF (*H2)[]);
void g_to_g2(GF g[GOPPA_T + 1], GF g2[2*GOPPA_T + 1]);
void restrict_pcheck2(GF (*H2)[CODE_N], uint8_t (*H2prime)[CODE_N]);
void decode(GF S[2*GOPPA_T], GF sigma[GOPPA_T + 1]);
#endif
