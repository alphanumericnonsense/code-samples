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
# Below is the additional functionality for the CCAKEM.
#-----------------------------------------------------------------------------------#

import pke
from random import random
import os

# N = 256
# Q = 3329

def kem_key_gen(d = None, z=None, level=5):
    """
    return (pk,sk) as a pair of bytearrays.
    the secret key contains the public key in the KEM specs,
    along with a hash of pk and a random 32-byte message for use
    in implicite rejection in the ciphertext comparison in decapsulation.
    """
    if z is None:
        z = os.urandom(32) # bytearray([int(256*random()) for i in range(32)]) # 32 random bytes
    pk, sk0 = pke.pke_key_gen(d, level) # CPA key generation
    sk = sk0 + pk + pke.H_hash(pk) + z # add randomness and hashes to secret key
    return (pk, sk)

def kem_encaps(pk, keylength, m=None, level=5):
    """
    return (c,K) ciphertext c and derived key K of length keylength.
    c and K are bytearrays.
    """
    if m is None:
        m = os.urandom(32) # bytearray([int(256*random()) for i in range(32)]) # 32 random bytes
    m = pke.H_hash(m)
    Kbar_r = pke.G_hash(m + pke.H_hash(pk))
    Kbar = Kbar_r[:32]
    r = Kbar_r[32:]
    c = pke.pke_encrypt(pk, m, r, level)
    K = pke.KDF(Kbar + pke.H_hash(c), keylength)
    return (c, K)

def kem_decaps(sk, c, keylength, level=5):
    """
    return derived key K, a bytearray of length keylength.
    """
    params = pke.param_dict(level)
    l = 12*256*params['k']//8
    sk0 = sk[:l]
    pk = sk[l:2*l+32]
    h = sk[2*l+32:2*l+64]
    z = sk[2*l+64:]
    m_prime = pke.pke_decrypt(sk0, c, level)
    Kbarprime_rprime = pke.G_hash(m_prime + h)
    Kbar_prime = Kbarprime_rprime[:32]
    r_prime = Kbarprime_rprime[32:]
    c_prime = pke.pke_encrypt(pk, m_prime, r_prime, level)
    if c == c_prime:
        K = pke.KDF(Kbar_prime + pke.H_hash(c), keylength)
    else:
        K = pke.KDF(z + pke.H_hash(c), keylength)
    return K

def kem_export_keys(pk, sk, name, level = 5):
    """
    Export keys to "name-pke-sec" and "name-pke-pub"
    """
    with open(name + f"-kem{level}-sec", "wb") as f:
        f.write(sk)
    with open(name + f"-kem{level}-pub", "wb") as f:
        f.write(pk)
