# AES notes

To encrypt or decrypt (AES128-ECB), run

> `./aestxt <keyfile> <plainfile> <cipherfile> <mode>`

from the command line e.g.

> `./aestxt key.txt plain.txt cipher e`

or

>  `./aestxt key.txt plain.txt cipher d`

to encrypt or decrypt respectively.

- The keyfile should be plain text with the key in hex (128 bits or 16 bytes, i.e. 32 hex characters), e.g.

  > `00112233445566778899aabbccddeeff`

- The file being written to will be created or overwritten (cipherfile for encryption, plainfile for decryption).

- The block cipher mode is ECB, i.e. the plain/cipher is broken into blocks and each block is (en/de)crypted separately.

- The plaintext is padded with `#` if not a multiple of the block size.

- `aestest` is just a functional test, encrypting/decrypting some test vectors.