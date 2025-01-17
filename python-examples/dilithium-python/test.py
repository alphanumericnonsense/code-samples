#--------------------------------------------------------------------------------------------------#
# DILITHIUM VERSION 3.1
#
# SPECS: https://pq-crystals.org/dilithium/resources.shtml
#
# Reference implementation:
# https://github.com/pq-crystals/dilithium
#
# Random file used for testing, can do timing evaluation for keygen/sign/verify.
#
# Another python implementation:
# https://github.com/jack4818/dilithium-py
#--------------------------------------------------------------------------------------------------#

from dilithium import *
from random import random
from math import floor
import aes256_ctr_drbg
import CompactFIPS202 as FIPS
#------------------------------#
# testing
#------------------------------#

def expand_KAT_seed(seed):
    """
    expand the 48 byte seed from NIST submission KAT.
    uses AES256 in counter mode to generate all necessary randomness from seed.

    """
    X = aes256_ctr_drbg.AES256_CTR_DRBG(seed=seed)
    zeta = X.random_bytes(32) # used in KeyGen
    return zeta

def RandomMessage(length):
    """
    random binary message of given byte length
    """
    M = []
    for l in range(length):
        M.append(math.floor(256*random()))
    return bytearray(M)

def SignRandomLoop(num_iters, message_size, level = 5):
    """
    keygen, sign, verify some random messages
    """
    genavg = 0
    signavg = 0
    veravg = 0
    params = param_dict(level)
    print(f"Security level {level}")
    print(f"Message size = {message_size} bytes")
    print(f"Running {num_iters} iterations of keygen/sign/verify.")
    for i in range(num_iters):
        #print(f"\tMessage {i} keygen/sign/verify")
        M = RandomMessage(message_size)
        #print(M)
        t1 = time()
        pk, sk = KeyGen(level=level)
        t2 = time()
        #print("KeyGen time", t2-t1)
        sigma = Sign(sk, M, level = level)
        t3 = time()
        #print("Sign time", t3-t2)

        #-------------------------------------------#
        # perturb sig or message to see if it fails
        #-------------------------------------------#
        #sigma[2][0][0] = int(1)
        #sigma[1][0][0] += 1
        #x = sigma[0] + 1
        #sigma = (x,sigma[1],sigma[2])
        #M[0] = 1 - M[0]
        #-------------------------------------------#
        accept_sig  = Verify(pk, M, sigma, level = level)
        t4 = time()
        #print("Verify time", t4-t3)
        genavg += t2-t1
        signavg += t3-t2
        veravg += t4-t3
        #print("sig, size", sigma[0], math.log2(sigma[0]))
        #print("\tNumber of hints", VecCount1s(sigma[2], KK))
        #print("\tAccepted?", accept_sig)
        if accept_sig == True:
            print("\tSUCCESS")
        else:
            print("\tFAILURE")
    print("keygen, sign, verify time avg.", genavg/num_iters, signavg/num_iters, veravg/num_iters)
    print("")

if __name__ == "__main__":
    #ZETA  = 1753

    # # x = [0 for i in range(256)]
    # # y = FFT256(x, ZETA)
    # # z = IFFT256(y, ZETA)
    # # print(x)
    # # print(y)
    # # print(z)
    #
    # x = [1] + [0 for i in range(255)]
    # y = FFT256(x, ZETA)
    # z = IFFT256(y, ZETA)
    # print(x)
    # print(y)
    # print(z)
    #
    # x = [0] + [1] + [0 for i in range(254)]
    # y = FFT256(x, ZETA)
    # z = IFFT256(y, ZETA)
    # print(x)
    # print(y)
    # print(z)
    #
    # x = [floor(QQ*random()) for i in range(256)]
    # y = FFT256(x, ZETA)
    # z = IFFT256(y, ZETA)
    # #print(x)
    # #print(y)
    # #print(z)
    # print(x == z)
    #
    # d = 13

    # x = [1] + [0 for i in range(255)]
    # y = NTT(x)
    # z = INTT(y)
    # print(x)
    # print(y)
    # print(z)
    # print(x == z)
    # x = [0] + [1] + [0 for i in range(254)]
    # y = NTT(x)
    # z = INTT(y)
    # print(x)
    # print(y)
    # print(z)
    # print(x == z)
    # x = [floor(QQ*random()) for i in range(256)]# RqScale(2**d, RqOne()) #
    # y = NTT(x)
    # z = INTT(y)
    # print(x)
    # print(y)
    # print(z)
    # print(x == z)

    num_iters = 10
    message_size = 128
    level = 2
    SignRandomLoop(num_iters, message_size, level = level)
    num_iters = 10
    message_size = 128
    level = 3
    SignRandomLoop(num_iters, message_size, level = level)
    num_iters = 10
    message_size = 128
    level = 5
    SignRandomLoop(num_iters, message_size, level = level)

    # seed = bytes.fromhex('061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1')
    # pk, sk = KeyGen(zeta = expand_KAT_seed(seed), level = 2)
    # M = bytes.fromhex("0123456789abcdef")
    # sigma = Sign(sk, M, randomized = False, coins = None, level = 2)
    #print(ba2h(sigma))

    #XOF = FIPS.gen_Keccak(1344, 256, seed, 0x1F)
    #X = FIPS.gen_SHAKE256(seed)
    #for i in range(10000):
    #    print(next(X))

