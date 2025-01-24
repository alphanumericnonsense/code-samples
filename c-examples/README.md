# Assorted C snippets

Some random stuff in C as evidence that I did things in C.  None of the crypto should be used for cryptographic purposes; it was written for learning purposes only.

## Contents

- `huffman/`.  Huffman coding, does the bible as an example.  Just produces the codes, doesn't do (de)compression.
- `AES/`.  AES in ECB mode, nothing fancy, a little program to encrypt/decrypt text files.
- `ALTCODES/`.  Encoding/decoding alternant codes (subfield subcodes of generalized Reed--Solomon codes).  Can futz with parameters ($n$, $k$, $q$, lower bounds for $d$).  Includes some randomized testing.
- `KARATSUBA/`.  A little test of integer coefficient polynomial Karatsuba multiplication.  $O(n^{\log_23})$ as opposed to $O(n^2)$ schoolbook multiplication.  [Wikipedia](https://en.wikipedia.org/wiki/Karatsuba_algorithm).
- The three PQC KEM algorithms below aren't following specs, just a loose version of the algorithms, although they could be upgraded with some Keccak, bit/byte-packing, etc.
  - `KYBER/`.  Edit parameters for different security levels.  Generates keys and runs some test encryption/decryption.
  - `MCELIECE/`.  Two versions for small and large parameter sets.  Key generation is slow, need to speed up row reduction.  Generates keys and runs some test encryption/decryption.
  - `SABER/`.   Edit parameters for different security levels.  Generates keys and runs some test encryption/decryption.