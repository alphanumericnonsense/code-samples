#-------------------------------------------------------------------------------------------------#
# Frodo LWE (Don't wear the ring!)
# Convservative LWE scheme working with linear algebra over Z/QZ for Q = 2^15 or 2^16
# and not some big fancy ring with optimized multiplication (Saber, Kyber, ~NTRU)

# https://frodokem.org/files/FrodoKEM-specification-20210604.pdf
#-------------------------------------------------------------------------------------------------#

import numpy as np
from time import time
from random import random
from math import floor

#------------------------------------#
# Security levels
#
# CHI_TABLE[abs(i)] is 2^16 times the probability
# that the discrete Gaussian takes the value i.
# Used to build the CDF for inverse sampling.
#------------------------------------#

# low
N = 640
D = 15
SUPCHI = 12
Bconst = 2
CHI_TABLE = {0:9288, 1:8720, 2:7216, 3:5264, 4:3384, 5:1918, 6:958, 7:422, 8:164, 9:56, 10:17, 11:4, 12:1}

# medium
#N = 976
#D = 16
#SUPCHI = 10
#Bconst = 3
#CHI_TABLE = {0:11278, 1:10277, 2:7774, 3:4882, 4:2545, 5:1101, 6:396, 7:118, 8:29, 9:6, 10:1}

# high
# N = 1344
# D = 16
# SUPCHI = 6
# Bconst = 4
# CHI_TABLE = {0:18286, 1:14320, 2:6876, 3:2023, 4:364, 5:40, 6:2}
#------------------------------------#

Q = 2**D # working in Z/QZ

# other matrix dimensions
MBAR = 8
NBAR = 8

# CDF for discrete Gaussian
def CDF_BUILD():
    table = {}
    S = 0
    for i in range(-SUPCHI, SUPCHI + 1):
        S += CHI_TABLE[abs(i)] * (2**(-16))
        table[i] = S
    return table

CDF_TABLE = CDF_BUILD()

L = Bconst*MBAR*NBAR # length of messages, i.e. 128/192/256
#-------------------------------------------------------------------------------------------------#

#-----------------------------------#
# algebra
#-----------------------------------#
def matrix_add(X, Y):
    """
    add matrices X+Y mod Q
    """
    h = len(X)
    w = len(X[0])
    Z = [[0 for i in range(w)] for j in range(h)]
    for i in range(h):
        for j in range(w):
            Z[i][j] = (X[i][j] + Y[i][j]) % Q
    return Z

def matrix_sub(X, Y):
    """
    subtract matrices X-Y mod Q
    """
    h = len(X)
    w = len(X[0])
    Z = [[0 for i in range(w)] for j in range(h)]
    for i in range(h):
        for j in range(w):
            Z[i][j] = (X[i][j] - Y[i][j]) % Q
    return Z

def matrix_mult(X, Y):
    """
    multiply matrices X*Y mod Q
    """
    hx = len(X)
    wx = len(X[0])
    hy = len(Y)
    wy = len(Y[0])
    Z = [[0 for i in range(wy)] for j in range(hx)]
    for i in range(hx):
        for j in range(wy):
            for k in range(wx):
                Z[i][j] += (X[i][k]*Y[k][j]) % Q
                Z[i][j] = Z[i][j] % Q
    return Z

#-----------------------------------#
# encode/decode binary message
# to/from matrix over Z/QZ
#-----------------------------------#
def encode(mu):
    """
    takes bitstring mu in blocks of size B to elements of Z/QZ and fills
    MBAR-by-NBAR matrix, row-major
    """
    M = [[0 for i in range(NBAR)] for j in range(MBAR)]
    for i in range(0,L,Bconst):
        w = mu[i : i+Bconst]
        #print(w)
        k = bit_list_to_int(w)
        #print(k)
        k1 = k*(2**(D-Bconst))
        #print(k1)
        ind = i//Bconst
        M[ind // NBAR][ind % NBAR] = k1
    return M

def decode(M):
    mu = []
    for i in range(MBAR):
        for j in range(NBAR):
            k = M[i][j]
            kbits = int_to_bit_list(k)
            for b in kbits:
                mu.append(b)
    return mu

def bit_list_to_int(bits):
    """
    little-endian, 1101 -> 1+2+0+8
    """
    i = 0
    result = 0
    for b in bits:
        result += b*(2**i)
        i += 1
    return result

def int_to_bit_list(n):
    x = round(n*(2**Bconst/Q)) % 2**Bconst
    y = [int(b) for b in bin(x)[2:]]
    y.reverse()
    z = y + [0 for i in range(Bconst - len(y))]
    return z

#-----------------------------------#
# random generation
#-----------------------------------#
def flip():
    return round(random())

def uniform_mod_Q():
    return floor(Q*random())

def chi_sample():
    """
    inverse sampling for the discrete Gaussian

    In practice, should be timing-resistant
    """
    u = random()
    for i in range(-SUPCHI, SUPCHI + 1):
        if u < CDF_TABLE[i]:
            return i

def random_uniform_matrix(d1,d2):
    X = [[0 for i in range(d2)] for j in range(d1)]
    for i in range(d1):
        for j in range(d2):
            X[i][j] = uniform_mod_Q()
    return X

def random_chi_matrix(d1, d2):
    X = [[0 for i in range(d2)] for j in range(d1)]
    for i in range(d1):
        for j in range(d2):
            X[i][j] = chi_sample()
    return X

#-----------------------------------#
# key generation
#-----------------------------------#

def key_gen():
    A = random_uniform_matrix(N, N)
    S = random_chi_matrix(N, NBAR)
    E = random_chi_matrix(N, NBAR)
    B = matrix_add(matrix_mult(A, S), E)
    return (A, B, S)

#-----------------------------------#
# encryption
#-----------------------------------#
def encrypt(A, B, mu):
    Sp = random_chi_matrix(MBAR, N)
    Ep = random_chi_matrix(MBAR, N)
    Epp = random_chi_matrix(MBAR, NBAR)
    Bp = matrix_add(matrix_mult(Sp, A), Ep)
    V = matrix_add(matrix_mult(Sp, B), Epp)
    C1 = Bp
    C2 = matrix_add(V, encode(mu))
    return (C1, C2)

#-----------------------------------#
# decryption
#-----------------------------------#
def decrypt(C1, C2, S):
    M = matrix_sub(C2, matrix_mult(C1, S))
    mup = decode(M)
    return mup

#-----------------------------------#
# testing
#-----------------------------------#

def random_message():
    return [flip() for i in range(L)]

def random_trials(n_trials = 5, verbose = False):
    for i in range(n_trials):
        print("="*25)
        print(f"Trial {i+1}/{n_trials}")
        A, B, S = key_gen()
        mu = random_message()
        C1, C2 = encrypt(A, B, mu)
        mup = decrypt(C1, C2, S)
        if verbose:
            #print("A\n", A)
            #print("B\n", B)
            #print("S\n", S)
            print("message mu\n", mu)
            print("\nciphertext C1\n", C1)
            print("\nciphertext C2\n", C2)
            print("-"*25)
        if mup == mu:
            print("SUCCESS")
        else:
            x = sum([(mup[i]+mu[i]) % 2 for i in range(L)]) # XOR
            print(f"FAILURE {x}/{L}")
