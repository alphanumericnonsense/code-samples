#------------------------------------------------------------------------------------------------------#
# Saber MLWR public-key scheme
#
# Uses arithmetic in R_Q = Z[x]/(Q, x^N - 1) mostly, with coefficient-wise rounding modulo Q, P, T, and 2
#
# N = 256, Q = 2^13, P = 2^10, T = 2^? varies based on security
#
# Module is small-rank (2/3/4) over uniform_R_Q
#
# Errors are from centered binomial distributions (CBD)
#
# Speed boosted by efficient bit-shifting >> and extracting LSB (say % 2^k with a bitmask &"11...11").
#
# See https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf for specifications,
# https://www.esat.kuleuven.be/cosic/pqcrypto/saber/index.html for other resources
#
# Multiplication is slowest part, currently uses 4-way karatsuba. Times for keygen/encrypt/decrypt:
#
# 0.86 seconds for schoolbook
# 0.21 seconds for karatsuba2 (one recursion degrees 256 -> 128 -> 64 then schoolbook)
#------------------------------------------------------------------------------------------------------#

# imports
from random import random
from math import floor
from time import time
from copy import copy

# Security levels: increasing dimension and T, decreasing MU

# LightSaber
#N = 256
#L = 2
#P = 2**10
#EPS_P = 10
#Q = 2**13
#EPS_Q = 13
#T = 2**3
#EPS_T = 3
#MU = 10

#Saber
#N = 256
#L = 3
#P = 2**10
#EPS_P = 10
#Q = 2**13
#EPS_Q = 13
#T = 2**4
#EPS_T = 4
#MU = 8

#FireSaber
N = 256
L = 4
P = 2**10
EPS_P = 10
Q = 2**13
EPS_Q = 13
T = 2**6
EPS_T = 6
MU = 6

# some constants
H1 = [2**(EPS_Q-EPS_P-1) for i in range(N)]
H = [H1 for i in range(L)]
H2 = [2**(EPS_P-2) - 2**(EPS_P-EPS_T - 1) + 2**(EPS_Q - EPS_P - 1) for i in range(N)]

#------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
# algebra, etc.
#-------------------------------------------------------------------------------------------#
def ring_add(a,b,modulus):
    return [(a[i] + b[i]) % modulus for i in range(N)]

def ring_sub(a,b,modulus):
    return [(a[i] - b[i]) % modulus for i in range(N)]

def polymult64to256(a, b, modulus):
    """
    for multiplying pieces in karatsuba2_helper()
    schoolbook multiplication
    """
    c = [0 for i in range(256)]
    for i in range(64):
        for j in range(64):
            k = i + j
            c[k] += a[i]*b[j]
            c[k] = c[k] % modulus
    return c

def karatsuba2_helper(a,b,modulus):
    """
    Basically a recursive karatsuba2 call, but degrees are low so x^N = -1 doesn't matter.
    """
    a1 = copy(a[64:])
    a0 = copy(a[:64])
    b1 = copy(b[64:])
    b0 = copy(b[:64])

    a0b0 = polymult64to256(a0,b0,modulus)
    a1b1 = polymult64to256(a1,b1,modulus)
    a0plusa1 = [(a0[i] + a1[i]) % modulus for i in range(64)]
    b0plusb1 = [(b0[i] + b1[i]) % modulus for i in range(64)]
    plusprod = polymult64to256(a0plusa1,b0plusb1,modulus)
    a1b1shifted = [a1b1[(i - 128) % 256] for i in range(256)]
    middle = ring_sub(plusprod, ring_add(a0b0, a1b1, modulus), modulus)
    middleshifted = [middle[(i - 64) % 256] for i in range(256)]
    return ring_add(a0b0, ring_add(a1b1shifted, middleshifted, modulus), modulus)

def karatsuba2(a, b, modulus):
    """
    2-way karatsuba (fewer mults, more adds)
    (a1*y+a2)*(b1*y+b2) = a1*b1*y^2 +((a1 + a2)*(b1 + b2) - a1*b1 - a2*b2)*y + a2*b2

    Here y = x^128, y^2 = x^256 = -1
    """
    a1 = copy(a[128:])
    a2 = copy(a[:128])
    b1 = copy(b[128:])
    b2 = copy(b[:128])

    a1b1 = karatsuba2_helper(a1, b1, modulus)
    a2b2 = karatsuba2_helper(a2, b2, modulus)
    a1plusa2 = [(a1[i] + a2[i]) % modulus for i in range(128)]
    b1plusb2 = [(b1[i] + b2[i]) % modulus for i in range(128)]
    plusprod = karatsuba2_helper(a1plusa2, b1plusb2, modulus)

    x = ring_sub(a2b2, a1b1, modulus) # y^2 and y^0 terms

    # middle term, needs negacyclic x^128 shift
    y = ring_sub(plusprod, ring_add(a1b1, a2b2, modulus), modulus)

    z = [0 for i in range(256)]
    for i in range(256):
        if i < 128:
            z[i] = modulus - y[i + 128]
        else:
            z[i] = y[i - 128]
    return ring_add(x, z, modulus)


def ring_mult(a, b, modulus):
    return karatsuba2(a, b, modulus)

# old ring multiplication, karatsuba much faster
# def ring_mult(a,b,modulus):
#     """
#     basic N^2 schoolbook multiplication with reduction modulo x^N = -1 and modulus = 0
#     """
#     c = [0 for i in range(N)]
#     #print(a)
#     #print(b)
#     for i in range(N):
#         for j in range(N):
#             r = (i + j) % N
#             q = (i + j)//N
#             c[r] += a[i]*b[j]*(-1)**(q)
#             c[r] = c[r] % modulus
#     return c

def mat_vec_mult(M, x, modulus):
    y = [[0 for i in range(N)] for l in range(L)] # zero vector
    for l in range(L):
        y[l] = dot(M[l], x, modulus)
    return y

def transpose(M):
    """
    return transpose of the L-by-L matrix M as list of rows
    """
    Mt = []
    for l in range(L):
        Mt.append([M[i][l] for i in range(L)])
    return Mt

def dot(u,v, modulus):
    S = [0 for i in range(N)]
    w = [ring_mult(u[l], v[l], modulus) for l in range(L)]
    for l in range(L):
        S = ring_add(S, w[l], modulus)
    return S

def vec_add(x, y, modulus):
    return [ring_add(x[l], y[l], modulus) for l in range(L)]

def vec_sub(x, y, modulus):
    return [ring_sub(x[l], y[l], modulus) for l in range(L)]

def ring_shift(a, d):
    """
    params: ring element a in R_Z, integer d
    return: ring element in R_(Z-d) with components right-shifted d places, i.e. xi//2^d
    """
    return [ai >> d for ai in a]

def vec_shift(x, d):
    """
    params: vector x, integer d
    return: vector with components coefficients right-shifted d places, i.e. xi//2^d
    """
    return [ring_shift(xi, d) for xi in x]

def ring_mod(a, modulus):
    """
    return a with coeffs reduced mod modulus
    """
    return [ai % modulus for ai in a]

def vec_mod(x, modulus):
    """
    return x with entries' coefficients reduced mod modulus
    """
    return [ring_mod(xi,modulus) for xi in x]

#-------------------------------------------------------------------------------------------#
# random value generation
#-------------------------------------------------------------------------------------------#
def uniform_R_Q():
    """
    return uniformly random list of 256 integers in [0,Q)
    """
    return [floor(Q*random()) for i in range(N)]

def uniform_matrix():
    """
    return L-by-L matrix with uniform_R_Q() entries
    matrix is a list of rows
    """
    M = []
    for l1 in range(L):
        M.append([uniform_R_Q() for l2 in range (L)])
    return M

def flip():
    """
    return 0 or 1 with probability 1/2
    """
    return floor(2*random())

def CBD():
    """
    returns an integer in [-MU/2, MU/2] as a difference of sums of MU {0,1}-Bernoulli trials
    """
    return sum([flip() for i in range(MU)]) - sum([flip() for i in range(MU)])

def ring_CBD():
    """
    return ring element with CBD() entries
    """
    return [CBD() for i in range(N)]

def vec_CBD():
    """
    return vector with CBD entry coefficients
    """
    return [ring_CBD() for l in range(L)]
def random_message():
    """
    return a uniformly random 256-bit binary string
    """
    return [flip() for i in range(N)]

#-------------------------------------------------------------------------------------------#
# key generation
#-------------------------------------------------------------------------------------------#
def key_gen():
    A = uniform_matrix()

    s = vec_CBD()

    #print(s)

    b = vec_add(mat_vec_mult(transpose(A), s, Q), H, Q)
    b = vec_shift(b, EPS_Q - EPS_P)

    return (A,b,s)

#-------------------------------------------------------------------------------------------#
# encryption
#-------------------------------------------------------------------------------------------#
def encrypt(A,b,m):
    s1 = vec_CBD()

    b1 = vec_add(mat_vec_mult(A,s1, Q), H, Q)
    b1 = vec_shift(b1, EPS_Q - EPS_P)

    v1 = dot(b, vec_mod(s1, P), P)
    c_m = ring_sub(ring_add(v1, H1, P), [(2**(EPS_P - 1) * m[i]) % P for i in range(N)], P)
    c_m = ring_shift(c_m, EPS_P - EPS_T)

    return (c_m, b1)

#-------------------------------------------------------------------------------------------#
# decryption
#-------------------------------------------------------------------------------------------#
def decrypt(c_m, b1, s):
    v = dot(b1, vec_mod(s, P), P)

    m1 = ring_add(ring_sub(v, [2**(EPS_P-EPS_T) * c_mi for c_mi in c_m], P), ring_mod(H2, P), P)
    m1 = ring_shift(m1, EPS_P - 1)

    return m1

#-------------------------------------------------------------------------------------------#
# testing
#-------------------------------------------------------------------------------------------#

def random_trials(n_trials = 5):
    for i in range(n_trials):
        print("="*25)
        print(f"Trial {i+1}/{n_trials}:")
        m = random_message()
        #print("m:\n",m)
        A, b, s = key_gen()
        #print("A:\n", A)
        #print("b:\n", b)
        #print("s:\n", s)
        c_m, b1 = encrypt(A,b,m)
        #print("c_m:\n", c_m)
        #print("b1:\n", b1)
        m1 = decrypt(c_m,b1,s)
        #print("m1:\n", m1)
        if m1 == m:
            print("SUCCESS")
        else:
            fails = sum([(m1[i]+m[i]) % 2 for i in range(N)])
            print(f"FAILURE {fails}")

#-------------------------------------------------------------------------------------------#
# some janky 4-way karatsuba multiplication without any recursion,
# slower than depth 2 karatsuba 2-way.
# keeping it around just in case...
#-------------------------------------------------------------------------------------------#
# def karatsuba4(a,b,modulus):
#     """
#     4-way karatsuba multiplication in Z[x]/(modulus, x^526 + 1)
#     """
#     # lower degree pieces
#     a3 = copy(a[192:])
#     a2 = copy(a[128:192])
#     a1 = copy(a[64:128])
#     a0 = copy(a[:64])
#     b3 = copy(b[192:])
#     b2 = copy(b[128:192])
#     b1 = copy(b[64:128])
#     b0 = copy(b[:64])

#     # diagonal terms
#     a3b3 = polymult64to256(a3,b3,modulus)
#     a2b2 = polymult64to256(a2,b2,modulus)
#     a1b1 = polymult64to256(a1,b1,modulus)
#     a0b0 = polymult64to256(a0,b0,modulus)

#     # additions before multiplications
#     a0plusa1 = [(a0[i] + a1[i]) % modulus for i in range(64)]
#     b0plusb1 = [(b0[i] + b1[i]) % modulus for i in range(64)]

#     a0plusa2 = [(a0[i] + a2[i]) % modulus for i in range(64)]
#     b0plusb2 = [(b0[i] + b2[i]) % modulus for i in range(64)]

#     a0plusa3 = [(a0[i] + a3[i]) % modulus for i in range(64)]
#     b0plusb3 = [(b0[i] + b3[i]) % modulus for i in range(64)]

#     a1plusa2 = [(a1[i] + a2[i]) % modulus for i in range(64)]
#     b1plusb2 = [(b1[i] + b2[i]) % modulus for i in range(64)]

#     a1plusa3 = [(a1[i] + a3[i]) % modulus for i in range(64)]
#     b1plusb3 = [(b1[i] + b3[i]) % modulus for i in range(64)]

#     a2plusa3 = [(a2[i] + a3[i]) % modulus for i in range(64)]
#     b2plusb3 = [(b2[i] + b3[i]) % modulus for i in range(64)]

#     # more complicated terms
#     y6 = a3b3
#     y5 = ring_sub(polymult64to256(a2plusa3, b2plusb3, modulus), ring_add(a2b2, a3b3, modulus), modulus)
#     y4 = ring_sub(ring_add(polymult64to256(a1plusa3, b1plusb3, modulus), a2b2, modulus), ring_add(a1b1, a3b3, modulus), modulus)

#     # pieces for y3
#     y3_1 = ring_sub(polymult64to256(a0plusa3, b0plusb3, modulus), ring_add(a0b0, a3b3, modulus), modulus)
#     y3_2 = ring_sub(polymult64to256(a1plusa2, b1plusb2, modulus), ring_add(a1b1, a2b2, modulus), modulus)
#     y3 = ring_add(y3_1, y3_2, modulus)

#     y2 = ring_sub(ring_add(polymult64to256(a0plusa2, b0plusb2, modulus),a1b1, modulus), ring_add(a0b0, a2b2, modulus), modulus)
#     y1 = ring_sub(polymult64to256(a0plusa1, b0plusb1, modulus),ring_add(a0b0, a1b1, modulus), modulus)
#     y0 = a0b0

#     # putting it together, need some negacyclic shifts
#     z = [0 for i in range(256)]
#     z6 = [-y6[(i-128) % N] for i in range(N)]
#     z = ring_add(z, z6, modulus)
#     z5 = [-y5[(i-64) % N] for i in range(N)]
#     z = ring_add(z, z5, modulus)
#     z4 = [-y4[i] for i in range(N)]
#     z = ring_add(z, z4, modulus)
#     z3 = [-y3[(i-192) % N] if i < 192 else y3[(i-192) % N] for i in range(N)]
#     z = ring_add(z, z3, modulus)

#     z2 = [y2[(i-128) % N] for i in range(N)]
#     z = ring_add(z, z2, modulus)
#     z1 = [y1[(i-64) % N] for i in range(N)]
#     z = ring_add(z, z1, modulus)
#     z0 = y0
#     z = ring_add(z, z0, modulus)
#     #print(z0,"\n")
#     #print(z1,"\n")
#     #print(z2,"\n")
#     #print(z3,"\n")
#     #print(z4,"\n")
#     #print(z5,"\n")
#     #print(z6,"\n")
#     #print(z,"\n")
#     return z
