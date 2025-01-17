import numpy as np
from random import random
from time import time
import copy
from math import floor, ceil, log2

def get_bit(n,i):
    """
    returns the ith bit of the integer n
    """
    return (n >> i) & 1

#------------------------------------------------------------------------------------------#
# GF(2^DIM) class
#
# Finite field elements are integers in [0, 2^DIM) and bitwise operations are used,
# e.g. 0 = 0, 1 = 1, x = 2, x+1 = 3, x^2 = 4, ..., x^{DIM-1}+x^{DIM-2}+...+x+1 = 2^{DIM}-1
#------------------------------------------------------------------------------------------#

# some choices for defining polynomial
# monic degree DIM to degree 0, left-to-right
MIN_POLY_LIST = ["1", "10", "111", "1011", "10011", "100101", "1000011", "10000011", "100011011", "1100000001", "10000001001", "100000000101", "1000000001001", "10000000011011", "10110011100001"]

class GF:

    def __init__(self, coeffs, min_poly):
        """
        initialize with a non-negative integer
        """
        self.min_poly = min_poly
        self.dim = len(self.min_poly) - 1
        self.mask = 2**self.dim - 1
        self.coeffs = coeffs & self.mask
        self.reducer = int(self.min_poly[1:], 2)

    def __add__(self, other):
        assert self.min_poly == other.min_poly
        return GF((self.coeffs ^ other.coeffs) & self.mask)

    def __sub__(self, other):
        """
        addition and subtraction are the same in char = 2
        """
        return self + other

    def __neg__(self):
        """
        identity in characteristic two
        """
        return self

    def __mul__(self, other):
        """
        MSE algorithm, could use something else (e.g. LUT)
        """
        assert self.min_poly == other.min_poly
        s = 0
        a = self.coeffs
        b = other.coeffs
        for i in range(self.dim - 1, -1, -1):
            eps = get_bit(s, self.dim - 1)
            s = ((s << 1) & self.mask) ^ (eps * self.reducer) ^ (get_bit(b, i) * a)
            #s = ((s << 1) & MASK) ^ ((eps * REDUCER) & MASK) ^ ((get_bit(b, i) * a) & MASK)
        return GF(s & self.mask, self.min_poly)

    def __pow__(self, e):
        """
        square and multiply recursively
        """
        if e == 0:
            return GF(1, self.min_poly)
        elif e == 1:
            return self
        elif e < 0:
            return self.inv()**(-e)
        else:
            if abs(e) % 2 == 0:
                return (self*self)**(e//2)
            else:
                return self*(self*self)**(e//2)

    def __eq__(self, other):
        assert self.min_poly == other.min_poly
        if self.coeffs ^ other.coeffs == 0:
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return format(self.coeffs, f'0{self.dim}b')

    def fermat_inv(self):
        """
        using Fermat x^q = x, fixes 0
        """
        assert self != GF(0, self.min_poly)
        return self**(2**self.dim - 2)

    def euclid_inv(self):
        """
        inversion using Euclidean algorithm
        """
        assert self != GF(0, self.min_poly)
        r = self.coeffs
        u = 1
        s = int(self.min_poly, 2)
        v = 0
        d = 0
        m = self.dim
        for i in range(1,2*self.dim + 1):
            if (r >> m) & 1 == 0:
                r <<= 1
                u <<= 1
                d += 1
            else:
                if (s >> m) & 1 == 1:
                    s = s ^ r
                    v = v ^ u
                s <<= 1
                if d == 0:
                    r, s = s, r
                    u, v = v << 1, u
                    d = 1
                else:
                    u >>= 1
                    d -= 1

        return GF(u, self.min_poly)

def random_GF(min_poly):
    """
    return a uniformly random field element
    """
    r = floor(random() * 2**(len(min_poly)-1))
    return GF(r, min_poly)

def random_nz_GF(min_poly):
    """
    return a uniformly random field element
    """
    r = floor(random() * (2**(len(min_poly)-1)-1))
    return GF(r+1, min_poly)

def fermat_inv_test(num_trials):
    times = []
    for i in range(2, len(MIN_POLY_LIST)):
        min_poly = MIN_POLY_LIST[i]
        t0 = time()
        for j in range(num_trials):
            x = random_nz_GF(min_poly)
            x_inv = x.fermat_inv()
            # if x*x_inv != GF(1, min_poly):
            #     print("FUCK")
        t1 = time()
        times.append(t1 - t0)
    return times

def euclid_inv_test(num_trials):
    times = []
    for i in range(2, len(MIN_POLY_LIST)):
        min_poly = MIN_POLY_LIST[i]
        t0 = time()
        for j in range(num_trials):
            x = random_nz_GF(min_poly)
            x_inv = x.euclid_inv()
            # if x*x_inv != GF(1, min_poly):
            #     print(x, x_inv)
        t1 = time()
        times.append(t1 - t0)
    return times

# num_trials = 100000
# F = fermat_inv_test(num_trials)
# E = euclid_inv_test(num_trials)
#
# import matplotlib.pyplot as plt
# plt.figure(figsize=(16,9))
# plt.title("random inversion\n100k trials per field")
# plt.xlabel("dimension m")
# plt.ylabel("total time (s)")
# plt.plot([i for i in range(2, len(MIN_POLY_LIST))], F, label="fermat")
# plt.plot([i for i in range(2, len(MIN_POLY_LIST))], E, label="euclid")
# plt.legend(loc="upper left")
# plt.savefig("fermat-vs-euclid.png")
# plt.show()
