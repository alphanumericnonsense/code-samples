#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
// types, constants
typedef uint8_t WORD[4];

const uint8_t SBOX[256] =
{
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

const uint8_t SBOXINV[256] =
{
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d
};

const WORD RCON[11] = // x^(i-1) mod x^8 + x^4 + x^3 + x + 1, first word doesn't matter.
{
    {0,0,0,0},
    {1,0,0,0},
    {2,0,0,0},
    {4,0,0,0},
    {8,0,0,0},
    {16,0,0,0},
    {32,0,0,0},
    {64,0,0,0},
    {128,0,0,0},
    {27,0,0,0},
    {54,0,0,0}
};

// finite field arithmetic GF(2^8)

uint8_t GF8Mult (uint8_t a, uint8_t b) {
    int DIM = 8;
    uint8_t MASK = (1 << DIM) - 1; // 2^DIM - 1
    uint8_t REDUCER = 27; // x^8 = x^4 + x^3 + x + 1 = 16+8+2+1

    uint8_t s = 0;
    int eps = 0;
    for (int i = DIM - 1; i > -1; i--) {
        eps = (s >> (DIM - 1)) & 1;
        s = ((s << 1) & MASK) ^ (eps * REDUCER) ^ (((b >> i) & 1) * a);
    }
    return s & MASK;
}

uint8_t GF8Inv (uint8_t a);
uint8_t GF8Pow (uint8_t a, int e) {
    //printf("in gf8pow: %d %d\n", a, e);
    if (e == 0) {
        return 1;
    } else if (e == 1) {
        return a;
    } else if (e < 0) {
        return GF8Pow(GF8Inv(a), -e);
    } else {
        if (e % 2 == 0) {
            return GF8Pow(GF8Mult(a, a), e/2);
        } else {
            return GF8Mult(a, GF8Pow(GF8Mult(a, a), e/2));
        }
    }
}

uint8_t GF8Inv (uint8_t a) {
    return GF8Pow(a, 254);
}

// stupid "hexchar array <-> byte array" conversion shit

uint8_t hchar2byte (char x) {
    switch (x) {
        case '0':
            return 0;
            break;
        case '1':
            return 1;
            break;
        case '2':
            return 2;
            break;
        case '3':
            return 3;
            break;
        case '4':
            return 4;
            break;
        case '5':
            return 5;
            break;
        case '6':
            return 6;
            break;
        case '7':
            return 7;
            break;
        case '8':
            return 8;
            break;
        case '9':
            return 9;
            break;
        case 'a':
            return 10;
            break;
        case 'b':
            return 11;
            break;
        case 'c':
            return 12;
            break;
        case 'd':
            return 13;
            break;
        case 'e':
            return 14;
            break;
        case 'f':
            return 15;
            break;
        default:
            printf("error in hchar2byte\n");
            return 0;
            break;
    }
}

char halfbyte2h (uint8_t x) {
    switch (x) {
        case 0:
            return '0';
            break;
        case 1:
            return '1';
            break;
        case 2:
            return '2';
            break;
        case 3:
            return '3';
            break;
        case 4:
            return '4';
            break;
        case 5:
            return '5';
            break;
        case 6:
            return '6';
            break;
        case 7:
            return '7';
            break;
        case 8:
            return '8';
            break;
        case 9:
            return '9';
            break;
        case 10:
            return 'a';
            break;
        case 11:
            return 'b';
            break;
        case 12:
            return 'c';
            break;
        case 13:
            return 'd';
            break;
        case 14:
            return 'e';
            break;
        case 15:
            return 'f';
            break;
        default:
            printf("error in halfbyte2h\n");
            return '0';
            break;
    }
}

uint8_t hh2byte (char a, char b) {
    return (hchar2byte(a) << 4) | (hchar2byte(b));
}

void hex_to_bytes (char hstring[], uint8_t bytes[], int L) {
    for (int i = 0; i < L; i++) {
        bytes[i] = hh2byte(hstring[2*i], hstring[2*i+1]);
    }
}

void bytes_to_hex (char hstring[], uint8_t bytes[], int L) {
    for (int i = 0; i < L; i++) {
        hstring[2*i] = halfbyte2h(bytes[i] >> 4);
        hstring[2*i+1] = halfbyte2h(bytes[i] & 15);
    }
}

// testing stuff

void print_state(uint8_t state[4][4]) {
    char hstring4[9];
    hstring4[9] = '\0';
    uint8_t row[4];
    printf("\n");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            row[j] = state[i][j];
        }
        bytes_to_hex(hstring4, row, 4);
        printf("%s\n", hstring4);
    }
    printf("\n");
}

/***********************************************************************
 * One-block AES-128 encryptin/decryption.
 * main() has some test vectors from NIST FIPS 197
 * Have some stupid bytes <-> hexstrings functions in AES.h
 * to make testing/printing easier.
 **********************************************************************/

void AddRoundKey (uint8_t state_in[4][4], uint8_t state_out[4][4], WORD schedule[], int round) {
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            state_out[row][col] = state_in[row][col] ^ schedule[4*round+col][row];
        }
    }
}

uint8_t Sbox (uint8_t x);
void KeySchedule (uint8_t key[], WORD schedule[]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <4; j++) {
            schedule[i][j] = key[4*i+j];
        }
    }
    WORD temp;
    for (int i = 4; i < 44; i++) {
        for (int j = 0; j <4; j++) {
            temp[j] = schedule[i-1][j];
        }
        if (i % 4 == 0) {
            for (int j = 0; j < 4; j++) {
                temp[j] = Sbox(schedule[i-1][(j+1) % 4]) ^ RCON[i/4][j];
            }
        }
        for (int j = 0; j < 4; j++) {
            schedule[i][j] = schedule[i - 4][j] ^ temp[j];
        }
    }
}

/**********************************************************************/

// uint8_t Sbox (uint8_t x) {
//     return SBOX[x];
// }

// more complicated Sbox, inversion + affine
uint8_t Sbox (uint8_t x) {
    uint8_t bits[8];
    uint8_t outbits[8];
    uint8_t row0[8] = {1,0,0,0,1,1,1,1};// first row of circulant matrix
    uint8_t temp;
    uint8_t xinv = GF8Inv(x);
    int k;

    for (int i = 0; i < 8; i++) {
        bits[i] = xinv & 1;
        xinv = xinv >> 1;
        //printf("bits[%d]=%d\n", i, bits[i]);
    }
    for (int i = 0; i < 8; i++) {
        temp = 0;
        for (int j = 0; j < 8; j++) {
            k = (j-i) % 8;
            if (k < 0) {
                k = k + 8;
            }
            temp = temp ^ (bits[j] & row0[k]);
            //printf("%d", row0[k]);
        }
        outbits[i] = temp;
        //printf("\noutbits[%d]=%d\n", i, outbits[i]);
    }
    temp = 0;
    for (int i = 0; i < 8; i++) {
        temp = temp + (outbits[i] * (1 << i));
        //printf("%d ", temp);
    }
    //printf("\n");
    return temp ^ 0x63;
}


void SubBytes (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            state_out[i][j] = Sbox(state_in[i][j]);
        }
    }
}

void ShiftRows (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            state_out[i][j] = state_in[i][(j+i) % 4];
        }
    }
}

void MixCols (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    for (int col = 0; col < 4; col++) {
        state_out[0][col] = GF8Mult(2, state_in[0][col]) ^ GF8Mult(3, state_in[1][col]) ^ state_in[2][col] ^ state_in[3][col];
        state_out[1][col] = state_in[0][col] ^ GF8Mult(2, state_in[1][col]) ^ GF8Mult(3, state_in[2][col]) ^ state_in[3][col];
        state_out[2][col] = state_in[0][col] ^ state_in[1][col] ^ GF8Mult(2, state_in[2][col]) ^ GF8Mult(3, state_in[3][col]);
        state_out[3][col] = GF8Mult(3, state_in[0][col]) ^ state_in[1][col] ^ state_in[2][col] ^ GF8Mult(2, state_in[3][col]);
    }
}

void AESRoundEnc (uint8_t state_in[4][4], uint8_t state_out[4][4], WORD schedule[], int round) {
    uint8_t temp1[4][4];
    uint8_t temp2[4][4];
    SubBytes(state_in, temp1);

    //printf("Round %d after subbytes:\n", round);
    //print_state(temp1);

    ShiftRows(temp1, temp2);

    //printf("Round %d after shiftrows:\n", round);
    //print_state(temp2);

    MixCols(temp2, temp1);

    //printf("Round %d after mixcols:\n", round);
    //print_state(temp1);

    AddRoundKey(temp1, state_out, schedule, round);

    //printf("Round %d after addroundkey:\n", round);
    //print_state(state_out);
}

void InitialRoundEnc (uint8_t plain[], uint8_t state_out[4][4], WORD schedule[]) {
    // arrange plaintext into 4x4 byte array (column major)
    uint8_t state_in[4][4];
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            state_in[row][col] = plain[4*col+row];
        }
    }
    AddRoundKey(state_in, state_out, schedule, 0);
}

void FinalRoundEnc (uint8_t state_in[4][4], uint8_t cipher[], WORD schedule[]) {
    uint8_t temp1[4][4];
    uint8_t temp2[4][4];

    SubBytes(state_in, temp1);

    //printf("Final round after subbytes:\n");
    //print_state(temp1);

    ShiftRows(temp1, temp2);

    //printf("Final round after shiftrows:\n");
    //print_state(temp2);

    //MixColumns(state_in, state_out); NO MIXCOLS

    AddRoundKey(temp2, temp1, schedule, 10);

    // serialize ciphertext (column major)
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            cipher[4*col+row] = temp1[row][col];
        }
    }
}

/**********************************************************************/

// uint8_t SboxInv (uint8_t x) {
//     return SBOXINV[x];
// }

// more complicated Sbox, y =ax^(-1)+b, x = [a^(-1)(y+b)]^(-1)
uint8_t SboxInv (uint8_t y) {
    uint8_t bits[8];
    uint8_t outbits[8];
    uint8_t row0[8] = {0,0,1,0,0,1,0,1};// first row of circulant matrix
    uint8_t temp;
    //uint8_t xinv = GF8Inv(x);
    int k;

    uint8_t x = y ^ 0x63;

    for (int i = 0; i < 8; i++) {
        bits[i] = x & 1;
        x = x >> 1;
        //printf("bits[%d]=%d\n", i, bits[i]);
    }
    for (int i = 0; i < 8; i++) {
        temp = 0;
        for (int j = 0; j < 8; j++) {
            k = (j-i) % 8;
            if (k < 0) {
                k = k + 8;
            }
            temp = temp ^ (bits[j] & row0[k]);
            //printf("%d", row0[k]);
        }
        outbits[i] = temp;
        //printf("\noutbits[%d]=%d\n", i, outbits[i]);
    }
    temp = 0;
    for (int i = 0; i < 8; i++) {
        temp = temp + (outbits[i] * (1 << i));
        //printf("%d ", temp);
    }
    //printf("\n");
    return GF8Inv(temp);
}

void SubBytesInv (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            state_out[i][j] = SboxInv(state_in[i][j]);
        }
    }
}

void ShiftRowsInv (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    int k;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            k = (j - i) % 4;
            if (k < 0) {
                k = k + 4;
            }
            state_out[i][j] = state_in[i][k];
        }
    }
}

void MixColsInv (uint8_t state_in[4][4], uint8_t state_out[4][4]) {
    for (int col = 0; col < 4; col++) {
        state_out[0][col] = GF8Mult(14, state_in[0][col]) ^ GF8Mult(11, state_in[1][col]) ^ GF8Mult(13, state_in[2][col]) ^ GF8Mult(9, state_in[3][col]);
        state_out[1][col] = GF8Mult(9, state_in[0][col]) ^ GF8Mult(14, state_in[1][col]) ^ GF8Mult(11, state_in[2][col]) ^ GF8Mult(13, state_in[3][col]);
        state_out[2][col] = GF8Mult(13, state_in[0][col]) ^ GF8Mult(9, state_in[1][col]) ^ GF8Mult(14, state_in[2][col]) ^ GF8Mult(11, state_in[3][col]);
        state_out[3][col] = GF8Mult(11, state_in[0][col]) ^ GF8Mult(13, state_in[1][col]) ^ GF8Mult(9, state_in[2][col]) ^ GF8Mult(14, state_in[3][col]);
    }
}

void AESRoundDec (uint8_t state_in[4][4], uint8_t state_out[4][4], WORD schedule[], int round) {
    uint8_t temp1[4][4];
    uint8_t temp2[4][4];

    ShiftRowsInv(state_in, temp1);

    //printf("Round %d after shiftrowsinv:\n", round);
    //print_state(temp1);

    SubBytesInv(temp1, temp2);

    //printf("Round %d after subbytesinv:\n", round);
    //print_state(temp2);

    AddRoundKey(temp2, temp1, schedule, 10 - round);

    //printf("Round %d after addroundkeyinv:\n", round);
    //print_state(temp1);

    MixColsInv(temp1, state_out);

    //printf("Round %d after mixcolsinv:\n", round);
    //print_state(state_out);
}

void InitialRoundDec (uint8_t cipher[], uint8_t state_out[4][4], WORD schedule[]) {
    // arrange ciphertext into 4x4 byte array (column major)
    uint8_t state_in[4][4];
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            state_in[row][col] = cipher[4*col+row];
        }
    }
    AddRoundKey(state_in, state_out, schedule, 10);
}

void FinalRoundDec (uint8_t state_in[4][4], uint8_t plain[], WORD schedule[]) {
    uint8_t temp1[4][4];
    uint8_t temp2[4][4];

    ShiftRowsInv(state_in, temp1);

    //printf("Final round after shiftrowsinv:\n");
    //print_state(temp1);

    SubBytesInv(temp1, temp2);

    //printf("Final round after shiftrowsinv:\n");
    //print_state(temp2);

    AddRoundKey(temp2, temp1, schedule, 0);

    // serialize plaintext (column major)
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            plain[4*col+row] = temp1[row][col];
        }
    }
}

/**********************************************************************/

void AESEncrypt (uint8_t plain[], uint8_t cipher[], uint8_t key[]) {
    uint8_t state_in[4][4];
    uint8_t state_out[4][4];
    WORD schedule[44]; // 4*(10+1)

    KeySchedule(key, schedule);

    InitialRoundEnc(plain, state_in, schedule);

    //printf("Round 0 state:\n");
    //print_state(state_in);

    for (int i = 1; i < 10; i++) {
        AESRoundEnc(state_in, state_out, schedule, i);
        // switch in/out
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                state_in[i][j] = state_out[i][j];
            }
        }
    }
    FinalRoundEnc(state_out, cipher, schedule);
}

void AESDecrypt (uint8_t plain[], uint8_t cipher[], uint8_t key[]) {
    uint8_t state_in[4][4];
    uint8_t state_out[4][4];
    WORD schedule[44]; // 4*(10+1)

    KeySchedule(key, schedule);

    InitialRoundDec(cipher, state_in, schedule);

    //printf("Round 0 state:\n");
    //print_state(state_in);

    for (int i = 1; i < 10; i++) {
        AESRoundDec(state_in, state_out, schedule, i);
        // switch in/out
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                state_in[i][j] = state_out[i][j];
            }
        }
    }
    FinalRoundDec(state_out, plain, schedule);
}
