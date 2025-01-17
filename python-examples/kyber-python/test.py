#-----------------------------------------------------------------------------------#
# KYBER VERSION 3.02
#
# SPECS:  https://pq-crystals.org/kyber/resources.shtml
#
# Reference implementation:
# https://github.com/pq-crystals/kyber
#
# Follows the specs and passes KAT from the reference implementation;
# Beware sample order for ref and specs differs in one spot which affects KAT.
#
# Another python implementation:
# https://github.com/jack4818/kyber-py
#
# A simple test of encryption/decryption, encapsulate/decapsulate
# export keys commented out for git
#-----------------------------------------------------------------------------------#

#
# a simple test of encryption/decryption, encapsulate/decapsulate
# export keys commented out for git
#

import pke
import kem

num_tests = 5
keylength = 32
print("pke test")
for level in [1,3,5]:
    print(f"level {level}")
    pk, sk = pke.pke_key_gen(level=level)
    #pke.pke_export_keys(pk,sk,"test-keys/test",level)
    print(f"|pk| = {len(pk)}, |sk| = {len(sk)}")
    print(f"begin {num_tests} tests...")
    for i in range(num_tests):
        m = pke.random_message()
        c = pke.pke_encrypt(pk, m, level=level)
        m_prime = pke.pke_decrypt(sk, c, level)
        if m == m_prime:
            print("\tsuccessful encrypt/decrypt")
    print(f"|ciphertext| = {len(c)}, |message| = {len(m)}\n")

print("kem test")
print(f"chosen keylength, |K| = {keylength}")
for level in [1,3,5]:
    print(f"level {level}")
    pk, sk = kem.kem_key_gen(level=level)
    #kem.kem_export_keys(pk,sk,"test-keys/test",level)
    print(f"|pk| = {len(pk)}, |sk| = {len(sk)}")
    print(f"begin {num_tests} tests...")
    for i in range(num_tests):
        c, K = kem.kem_encaps(pk, keylength, level=level)
        K_prime = kem.kem_decaps(sk, c, keylength, level)
        if K == K_prime:
            print("\tsuccessful encaps/decaps")
    print(f"|ciphertext| = {len(c)}, |K| = {len(K)}\n")
