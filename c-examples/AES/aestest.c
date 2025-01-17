#include "AES.h"

int main () {
    uint8_t plain[16];
    uint8_t cipher[16];;
    WORD schedule [44];
    uint8_t key[16];
    char hstring16[33];//32 hex chars, 16 bytes
    hstring16[32] = '\0';// NULL char for C-string
    char hstring4[9];//8 hex chars, 4 bytes
    hstring4[8] = '\0';//NULL char for C-string

    /*//test GF(2^8) multiplication
    uint8_t gf;
    gf = GF8Mult(0x57, 0x02);
    printf("%x\n", gf);
    gf = GF8Mult(0x57, 0x04);
    printf("%x\n", gf);
    gf = GF8Mult(0x57, 0x08);
    printf("%x\n", gf);
    gf = GF8Mult(0x57, 0x10);
    printf("%x\n", gf);
    gf = GF8Mult(0x57, 0x13);
    printf("%x\n", gf);
    */

    /*//test GF8Pow, GF8Inv
    uint8_t x = 1;
    for (int i =  0; i < 256; i++) {
        printf("%d %x ?= %x\n", i, x, GF8Pow(2,i));
        x = GF8Mult(2, x);
    }
    */

    /*// test complicated Sbox
    for (int i = 0; i < 16; i++) {
        printf("%x\n", Sbox((uint8_t) i));
    }
    */

    // test encrypt 1
    hex_to_bytes("3243f6a8885a308d313198a2e0370734", plain, 16);
    hex_to_bytes("2b7e151628aed2a6abf7158809cf4f3c", key, 16);
    hex_to_bytes("3925841d02dc09fbdc118597196a0b32", cipher, 16);

    printf("Test encryption 1\n");

    printf("\tKey:\n");
    bytes_to_hex(hstring16, key, 16);
    printf("\t\t%s\n", hstring16);
    /*
    printf("\tKey schedule:\n");
    KeySchedule(key, schedule);
    for (int i = 0; i < 44; i++) {
        bytes_to_hex(hstring4, schedule[i], 4);
        printf("\t\t%s\n", hstring4);
    }
    */
    printf("\tPlaintext in:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    printf("\tCiphertext should be:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    AESEncrypt(plain, cipher, key);

    printf("\tEncryption is:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    // test decrypt 1
    printf("Test decryption 1\n");

    printf("\tKey:\n");
    bytes_to_hex(hstring16, key, 16);
    printf("\t\t%s\n", hstring16);
    /*
     printf("\tKey schedule:\n"*);
     KeySchedule(key, schedule);
     for (int i = 0; i < 44; i++) {
         bytes_to_hex(hstring4, schedule[i], 4);
         printf("\t\t%s\n", hstring4);
    }
    */
    printf("\tCiphertext in:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    printf("\tPlaintext should be:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    AESDecrypt(plain, cipher, key);

    printf("\tDecryption is:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    // test encrypt 2
    hex_to_bytes("00112233445566778899aabbccddeeff", plain, 16);
    hex_to_bytes("000102030405060708090a0b0c0d0e0f", key, 16);
    hex_to_bytes("69c4e0d86a7b0430d8cdb78070b4c55a", cipher, 16);

    printf("Test encryption 2\n");

    printf("\tKey:\n");
    bytes_to_hex(hstring16, key, 16);
    printf("\t\t%s\n", hstring16);
    /*
     printf("\tKey schedule:\n");
     KeySchedule(key, schedule);
     for (int i = 0; i < 44; i++) {
         bytes_to_hex(hstring4, schedule[i], 4);
         printf("\t\t%s\n", hstring4);
    }
    */
    printf("\tPlaintext in:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    printf("\tCiphertext should be:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    AESEncrypt(plain, cipher, key);

    printf("\tEncryption is:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    // test decrypt 2
    printf("Test decryption 2\n");

    printf("\tKey:\n");
    bytes_to_hex(hstring16, key, 16);
    printf("\t\t%s\n", hstring16);
    /*
        printf("\tKey schedule:\n");
    KeySchedule(key, schedule);
    for (int i = 0; i < 44; i++) {
        bytes_to_hex(hstring4, schedule[i], 4);
        printf("\t\t%s\n", hstring4);
    }
    */
    printf("\tCiphertext in:\n");
    bytes_to_hex(hstring16, cipher, 16);
    printf("\t\t%s\n", hstring16);

    printf("\tPlaintext should be:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    AESDecrypt(plain, cipher, key);

    printf("\tDecryption is:\n");
    bytes_to_hex(hstring16, plain, 16);
    printf("\t\t%s\n", hstring16);

    return 0;
}
