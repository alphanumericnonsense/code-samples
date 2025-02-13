import os
import math

def witness(a, n, bitlen):
    """
    returns True if n is composite
    false if no witness detected by:
        - fermat or
        - non-trivial sqrt of 1 in squaring chain of modular exponentiation
    """
    d = 1
    for i in range(bitlen, -1, -1):
        b = ((n - 1) >> i) & 1
        x = d
        d = (d**2) % n
        if (d == 1) and (x != 1) and (x != n - 1):
            #print(n, a, "nontrivsqrt1")
            return True
        if b == 1:
            d = (d*a) % n
    if d != 1:
        #print(n, a, "fermatfail")
        return True
    return False

def miller_rabin(n, s):
    """
    miller-rabin primality test, s iterations testing primality of n
    """
    bitlen = math.ceil(math.log2(n-1))
    i = 0
    while i < s:
        a = int.from_bytes(os.urandom((math.ceil(math.log2(n)) + 7)//8), "big")
        if 0 < a < n:
            i += 1
            #print(n, a)
            if witness(a, n, bitlen):
                return False
    return True

def random_prime(bytewidth, s):
    found = False
    while True:
        n = int.from_bytes(os.urandom(bytewidth), "big")
        isprime = miller_rabin(n, s)
        if isprime:
            return n

if __name__ == "__main__":
    s = 4
    pmax = 100
    print(f"probable primes in range(2,{pmax}), {s} MR iterations:")
    for n in range(2,pmax):
        isprime = miller_rabin(n, s)
        if isprime:
            print('\t', n, "is probably prime")

    bytewidth = 64
    print("")
    print(random_prime(bytewidth, s), f"is probably a {bytewidth}-byte prime, {s} MR iterations")
