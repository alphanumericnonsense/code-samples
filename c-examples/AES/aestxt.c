#include "AES.h"

void en_de_crypt_file(char keyfile[], char plainfile[], char cipherfile[], char mode[]) {
    FILE* key_fp;
    FILE* plain_fp;
    FILE* cipher_fp;
    uint8_t key[16];
    char keystring[16];
    uint8_t blockin[16];
    uint8_t blockout[16];
    int num_char = 0;// number of characters in file
    int num_blocks = 0;// number of blocks to endecrypt
    int rembytes = 0;//partial block size for encryption
    int c;

    //get key
    key_fp = fopen(keyfile, "r");
    for (int i = 0; i < 32; i++) {
        keystring[i] = fgetc(key_fp);
    }
    fclose(key_fp);
    hex_to_bytes (keystring, key, 16); // key is now uint8_t array


    printf("*\tkeyfile: %s\n*\tplainfile: %s\n*\tcipherfile: %s\n*\tmode: %s\n", keyfile, plainfile, cipherfile, mode);
    //start encrypting or decrypting
    if (mode[0] == 'e') {
        // get file size
        plain_fp = fopen(plainfile, "r");
        do {
            c = fgetc(plain_fp);
            if (c != EOF) {
                num_char += 1;
            }
        } while (c != EOF);
        fclose(plain_fp);
        num_blocks = num_char*8/128;
        rembytes = num_char - num_blocks*16;

        printf("*\tnumber of characters encrypted: %d\n", num_char);
        printf("*\tnumber of full blocks encrypted: %d\n", num_blocks);
        printf("*\tnumber of characters in partial block: %d\n", rembytes);

        cipher_fp = fopen(cipherfile, "w");
        plain_fp = fopen(plainfile, "r");
        for (int i = 0; i < num_blocks; i++) {
            //read/load block to encrypt
            for (int j = 0; j < 16; j++) {
                c = fgetc(plain_fp);
                blockin[j] = (uint8_t)c;
            }
            //encrypt a block
            AESEncrypt(blockin, blockout, key);
            //write encrypted block
            for (int j = 0; j < 16; j++) {
                c = (int)blockout[j];
                fputc(c, cipher_fp);

                // write to hex
                //int c1 = halfbyte2h(c >> 4);
                //int c2 = halfbyte2h(c & 15);
                //fputc(c1, cipher_fp);
                //fputc(c2, cipher_fp);
            }
        }
        //padding if partial final block
        if (rembytes > 0) {
            for (int j = 0; j < 16; j++) {
                if (j < rembytes) {
                    c = fgetc(plain_fp);
                } else {
                    c = 35;//pad with #
                }
                blockin[j] = (uint8_t)c;
            }
            //encrypt a block
            AESEncrypt(blockin, blockout, key);
            //write encrypted block
            for (int j = 0; j < 16; j++) {
                c = (int)blockout[j];
                fputc(c, cipher_fp);
            }
        }
        fclose(cipher_fp);
        fclose(plain_fp);
    } else if (mode[0] == 'd') {
        // get file size
        cipher_fp = fopen(cipherfile, "r");
        do {
            c = fgetc(cipher_fp);
            if (c != EOF) {
                num_char += 1;
            }
        } while (c != EOF);
        fclose(cipher_fp);
        num_blocks = num_char*8/128;
        printf("*\tnumber of characters decrypted: %d\n", num_char);
        printf("*\tnumber of full blocks decrypted: %d\n", num_blocks);
        cipher_fp = fopen(cipherfile, "r");
        plain_fp = fopen(plainfile, "w");
        for (int i = 0; i < num_blocks; i++) {
            //read/load block to decrypt
            for (int j = 0; j < 16; j++) {
                c = fgetc(cipher_fp);
                blockin[j] = (uint8_t)c;
            }
            //decrypt a block
            AESDecrypt(blockout, blockin, key);
            //write decrypted block
            for (int j = 0; j < 16; j++) {
                c = (int)blockout[j];
                fputc(c, plain_fp);
            }
        }
        fclose(cipher_fp);
        fclose(plain_fp);
    } else {
        printf("Incorrect mode (choices are encrypt `e` or decrypt `d`)\n");
    }
    return;
}

int main (int argc, char* argv[]) {
    if (argc == 5) {
        en_de_crypt_file(argv[1], argv[2], argv[3], argv[4]);
    } else {
        printf("Too few arguments: keyfile, plainfile, cipherfile, mode (encrypt `e` or decrypt `d`).\n");
        return 0;
    }
    return 0;
}
