#------------------------------------------------------------------------------------------#
# Imports, constants
#------------------------------------------------------------------------------------------#

import numpy as np
from random import random
from time import time
import copy
from math import floor
#import sympy # for binary rref maybe

#---------------------------#
# some cryptographic levels
#---------------------------#

# key gen ~ 24 secs avg., encryption 0.8 secs, decryption 11 secs
MIN_POLY = "1000000001001" # monic degree 0 to degree DIM
DIM = 12
DEGREE = 64
CODE_DIM = 3488
EXTRA_MIN_POLY = [(64,1),(3,1),(1,1),(0,2)] # non-zero coefficient GF(j) in degree i for (i,j)

# key gen ~ 64 secs avg., encryption 1.5 secs, decryption 24 secs
# MIN_POLY = "10000000011011" # monic degree DIM to degree 0, left-to-right
# DIM = 13
# DEGREE = 96
# CODE_DIM = 4608
# EXTRA_MIN_POLY = [(96,1),(10,1),(9,1),(6,1),(0,1)] # non-zero coefficient GF(j) in degree i for (i,j)

# MIN_POLY = "10000000011011" # monic degree DIM to degree 0, left-to-right
# DIM = 13
# DEGREE = 128
# CODE_DIM = 6688
# EXTRA_MIN_POLY = [(128,1),(7,1),(2,1),(1,1),(0,1)] # non-zero coefficient GF(j) in degree i for (i,j)

# key gen in ~5 mins?, encryption 3 secs, decryption 45 secs
# MIN_POLY = "10000000011011" # monic degree DIM to degree 0, left-to-right
# DIM = 13
# DEGREE = 119
# CODE_DIM = 6960
# EXTRA_MIN_POLY = [(119,1),(8,1),(0,1)] # non-zero coefficient GF(j) in degree i for (i,j)

# MIN_POLY = "10000000011011" # monic degree DIM to degree 0, left-to-right
# DIM = 13
# DEGREE = 128
# CODE_DIM = 8192
# EXTRA_MIN_POLY = [(128,1),(7,1),(2,1),(1,1),(0,1)] # non-zero coefficient GF(j) in degree i for (i,j)

#---------------------------#
# some fields for testing
#---------------------------#

# GF(4)
# for super testing
# MIN_POLY = "111"
# DIM = 2
# DEGREE = 2
# CODE_DIM = 0
# EXTRA_MIN_POLY = [(2,1),(1,1),(0,2)] # non-zero coefficient GF(j) in degree i for (i,j)

# GF(2**7)
# MIN_POLY = "10000011"
# DIM = 7

# GF(2**9)
# MIN_POLY = "1000000011"
# DIM = 9

# GF(2**10)
# MIN_POLY = "10000001001"
# DIM = 10

# GF(2**11)
# MIN_POLY = "100000000101"
# DIM = 11

# GF(2**12)
# MIN_POLY = "1000000001001"
# DIM = 12

MASK = 2**DIM - 1 # used as a bitmask
REDUCER = int(MIN_POLY[1:], 2) # used in finite field multiplication

#------------------------------------------------------------------------------------------#
# GF(2^DIM) class
#
# Finite field elements are integers in [0, 2^DIM) and bitwise operations are used
#------------------------------------------------------------------------------------------#

class GF:

    def __init__(self, coeffs):
        """
        initialize with a non-negative integer
        """
        self.coeffs = coeffs & MASK

    def __add__(self, other):
        return GF((self.coeffs ^ other.coeffs) & MASK)

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
        s = 0
        a = self.coeffs
        b = other.coeffs
        for i in range(DIM - 1, -1, -1):
            eps = get_bit(s, DIM - 1)
            s = ((s << 1) & MASK) ^ (eps * REDUCER) ^ (get_bit(b, i) * a)
            #s = ((s << 1) & MASK) ^ ((eps * REDUCER) & MASK) ^ ((get_bit(b, i) * a) & MASK)
        return GF(s & MASK)

    def __pow__(self, e):
        """
        square and multiply recursively
        """
        if e == 0:
            return FLD_ONE()
        elif e == 1:
            return self
        elif e < 0:
            return self.inv()**(-e)
        else:
            if abs(e) % 2 == 0:
                return (self*self)**(e//2)
            else:
                return self*(self*self)**(e//2)

#     def __pow__(self, e):
#         """
#         square and multiply iteratively
#         """
#         if e == 0:
#             return FLD_ONE()
#         elif e == 1:
#             return self
#         elif e < 0:
#             return self.inv()**(-e)
#         else:
#             prod = GF(1)
#             bitstring = bin(e)[2:]
#             for c in bitstring:
#                 prod = prod * prod
#                 if c == '1':
#                     prod = prod*self
#             return prod

    def __eq__(self, other):
        if self.coeffs ^ other.coeffs == 0:
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return bin(self.coeffs)

    def inv(self):
        """
        return inverse of self using LUT
        """
        if self == FLD_ZERO():
            print("zero not invertible, returning None")
            return None
        else:
            return GF(INVERSES[self.coeffs])

    def fermat_inv(self):
        """
        using Fermat x^q = x, fixes 0
        """
        return self**(2**DIM - 2)

# useful functions using or used by GF
# constants, inversion LUT

def get_bit(n,i):
    """
    returns the ith bit of the integer n
    """
    return (n >> i) & 1

def random_GF():
    """
    return a uniformly random field element
    """
    r = floor(random() * (MASK + 1))
    return GF(r)

def  FLD_ZERO():
    return GF(0)

def FLD_ONE():
    return GF(1)

def inv_dict():
    """
    precompute inverses for the finite field

    returns: dictionary {n : ninv} with integer key/value pairs
    """
    return {i : GF(i).fermat_inv().coeffs for i in range(2**DIM)}

INVERSES = inv_dict() # generates LUT for inverses

#------------------------------------------------------------------------------------------#
# Polynomials over GF(2^DIM) class
#
# Polynomials are lists of field elements with trailing zero coefficients trimmed,
# i.e. 0 = [] has degree -1, 1 + x = [1,1] degree 1
#
# degree increases with index, [deg0, deg1, etc.]
#------------------------------------------------------------------------------------------#

class PolyFF(object):

    def __init__(self, coeffs):
        """
        coeffs should be a list of GF
        """
        # trim trailing zero coefficients
        i = len(coeffs)
        for c in coeffs[::-1]:
            if c != FLD_ZERO():
                break
            i -= 1
        self.coeffs = coeffs[:i]
        self.degree = len(self.coeffs) - 1 # zero poly [] has degree -1

    def __mul__(self, other):
        m = len(self.coeffs)
        n = len(other.coeffs)
        C = [FLD_ZERO() for k in range(m+n)]
        for i in range(m):
            for j in range(n):
                C[i+j] = C[i+j] + self.coeffs[i]*other.coeffs[j]
        return PolyFF(C)

    def __add__(self, other):
        m = len(self.coeffs)
        n = len(other.coeffs)
        C = [FLD_ZERO() for i in range(max(m,n))]
        if m >= n:
            for i in range(n):
                C[i] = self.coeffs[i] + other.coeffs[i]
            for i in range(n,m):
                C[i] = self.coeffs[i]
        else:
            for i in range(m):
                C[i] = self.coeffs[i] + other.coeffs[i]
            for i in range(m,n):
                C[i] = other.coeffs[i]
        return PolyFF(C)

    def __sub__(self, other):
        """
        same as addition in char = 2
        """
        return self + other

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __ne__(self, other):
        return not self.coeffs == other.coeffs

    def __neg__(self):
        """
        identity in char = 2
        """
        return self

    def __repr__(self):
        if self.coeffs == []:
            return "zero_poly"
        string = ""
        for i in range(self.degree+1):
            string += "deg "
            string += str(i)
            string += "\n"
            string += str(self.coeffs[i])
            string +="\n"
        return string

    def quot_rem(self, other):
        """
        usual polynomial division, probably something faster???
        This is really slow

        return: quotient, remainder such that self = quotient*other + remainder
        """
        quot = PolyFF([])
        rem = copy.deepcopy(self)
        if other.degree == -1:
            print("division by zero, returning None")
            return None
        while rem.degree >= other.degree:
            remlead = rem.coeffs[rem.degree]
            otherlead = other.coeffs[other.degree]
            scalar = remlead*(otherlead.inv())
            scalar0 = PolyFF([scalar])
            monomial = PolyFF([FLD_ZERO() for i in range(rem.degree-other.degree)] + [scalar])
            rem = rem - monomial*other
            quot = quot + monomial
        return (quot, rem)

# useful functions using/used by PolyFF

def mod_prod(a, b, f):
    """
    return: a*b mod f assuming degree of a, b less than deg(f)
            and monic f
    """
    s = PolyFF([]) # polynomial zero
    r = PolyFF(f.coeffs[:-1]) # for reducing x^deg(f) = sum of lower degree terms
    for i in range(b.degree+1):
        s = PolyFF([FLD_ZERO()] + s.coeffs) # multiply by x, i.e. shift
        if s.degree >= f.degree:
            s = PolyFF(s.coeffs[:-1]) + PolyFF([s.coeffs[-1]])*r
        s = s + PolyFF([b.coeffs[-i-1]])*a
    return s

def poly_EEA(f, g, early_stop = 0):
    """
    extended Euclidean algorithm for polynomials over finite fields
    NEED TO MAKE FASTER

    params: polynomials f and g with GF coefficients
    return: (a,b,d) such that af+bg = d = gcd(f,g)
            Note: if stopped early, stops when degree of d < early_stop
            so d != gcd(f,g) necessarily
    """
    if f == PolyFF([]) and g == PolyFF([]):
        print("0/0, return None")
        return None
    rold, rnew = f, g
    aold, anew = PolyFF([]), PolyFF([FLD_ONE()])
    bold, bnew = PolyFF([FLD_ONE()]), PolyFF([])
    i = 0
    while rnew != PolyFF([]):
        q, r = rold.quot_rem(rnew)
        #print("q",q)
        #print("r",r)
        aold, anew = anew, q*anew + aold
        bold, bnew = bnew, q*bnew + bold
        rold, rnew = rnew, r
        i += 1
        if rold.degree < early_stop:
            break
    if i % 2 == 0:
        #print("poly_EEA test")
        #print(bold*f-aold*g)
        #print(rold)
        return (bold, -aold, rold)
    else:
        #print("poly_EEA test")
        #print(-bold*f+aold*g)
        #print(rold)
        return (-bold, aold, rold)

def poly_gcd(f,g):
    """
    iterative binary algorithm (recursive version was too deep)
    """
    a = f
    b = g
    d = PolyFF([FLD_ONE()])

    while a != PolyFF([]) and b != PolyFF([]):
        a0 = a.coeffs[0]
        b0 = b.coeffs[0]
        #print(a.degree, b.degree)
        if a0 == FLD_ZERO():
            if b0 == FLD_ZERO():
                d = d*PolyFF([FLD_ZERO(), FLD_ONE()])
                a = PolyFF(a.coeffs[1:])
                b = PolyFF(b.coeffs[1:])
            else:
                a = PolyFF(a.coeffs[1:])
        elif b0 == FLD_ZERO():
            b = PolyFF(b.coeffs[1:])
        else:
            if a.degree >= b.degree:
                a = a - PolyFF([a.coeffs[0]*b.coeffs[0].inv()])*b
                a = PolyFF(a.coeffs[1:])
            else:
                b = b - PolyFF([b.coeffs[0]*a.coeffs[0].inv()])*a
                b = PolyFF(b.coeffs[1:])
    if a == PolyFF([]):
        if b == PolyFF([]):
            print("gcd(0,0) undefined, returning None")
            return None
        else:
            return d*b*PolyFF([b.coeffs[-1].inv()])
    elif b == PolyFF([]):
        return d*a*PolyFF([a.coeffs[-1].inv()])

def coprime(f,g):
    """
    params: polynomials f, g
    return: True if f, g coprime, False otherwise
    """
    d = poly_gcd(f,g)
    if d.degree == 0:
        return True
    else:
        return False

def eval_poly(f, x):
    """
    evaluate polynomial f at point x in GF(q) in an iterative (Horner) fashion
    """
    s = FLD_ZERO()
    d = f.degree
    for i in range(d + 1):
        s = s*x + f.coeffs[d - i]
    return s

def is_separable(f):
    """
    return True if gcd(f,f')=1, else False
    """
    fprime = PolyFF([f.coeffs[i] if i % 2 == 1 else FLD_ZERO() for i in range(f.degree+1)])
    return coprime(f, fprime)

#------------------------------------------------------------------------------------------#
# Key generation
# Encryption
# Decryption
# Testing
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
# key generation
#------------------------------------------------------------------------------------------#

def random_monic_poly(degree = DEGREE):
    """
    random monic polynomial of degree 'degree' over GF(q)
    """
    return PolyFF([random_GF() for i in range(degree)] + [FLD_ONE()])

def random_poly(degree):
    """
    random polynomial of degree 'degree' over GF(q)
    """
    return PolyFF([random_GF() for i in range(degree+1)])

def new_separable_goppa(L):
    """
    random monic separable goppa polynomial
    with zeros disjoint from the support L through rejection sampling.
    """
    disjoint = False
    while True:
        g = random_monic_poly()
        disjoint = True
        for l in L:
            if eval_poly(g,l) == FLD_ZERO():
                disjoint = False
                break
        if disjoint:
            if is_separable(g):
                return g

def new_irred_goppa():
    """
    bad implementation because I didn't want to generalize the GF class for other fields

    return: An irreducible monic Goppa polynomial as the minimal polynomial
            of a random element of GF(2^(DIM*DEG)).
            Uses fast_right_rref_GF to solve the relevant system of equations
            and check that the powers of the random element are linearly independent
    """
    # f is the defining polynomial for a degree DEGREE extenstion of GF(2^DIM)
    coeffs = [FLD_ZERO() for i in range(DEGREE + 1)]
    for pair in EXTRA_MIN_POLY:
        i, j = pair
        coeffs[i] = GF(j)
    f = PolyFF(coeffs)

    irred = False
    while not irred:
        # r is a stand-in for a random element of GF(2^DIM*DEGREE)
        r = random_poly(DEGREE - 1)
        #print("random\n",r)
        # generate powers of r modulo f
        power_list = [PolyFF([FLD_ONE()])]
        for i in range(DEGREE):
            power_list.append(mod_prod(r, power_list[i], f))

        # building a GF(2^m) matrix from coefficients of powers of r modulo f
        M = []
        for i in range(DEGREE):
            row = []
            for j in range(DEGREE + 1):
                k = DEGREE - 1 - i
                if power_list[j].degree >= k:
                    row.append(power_list[j].coeffs[k])
                else:
                    row.append(FLD_ZERO())
            M.append(row)

        # row reduce M to get the goppa polynomial as the last column if M full rank
        good, g = fast_right_rref_GF(M)
        if good:
            g.append(FLD_ONE())
            irred = True
        else:
            print("goppa poly reducible...")
    #print(g)
    return PolyFF(g)

def new_support():
    """
    return distinct list of CODE_DIM elements of GF(q)
    uses a random permutation of GF(q)
    """
    perm = np.random.permutation(2**DIM)
    return [GF(int(perm[i])) for i in range(CODE_DIM)]

def cannonical_pcheck(g, L):
    """
    generate cannonical parity check matrix for the goppa code
    """
    H = np.zeros((g.degree, len(L)))
    for j in range(len(L)):
        #print(f"column {j} of pcheck")
        x = eval_poly(g, L[j])
        if x == FLD_ZERO():
            print("L not disjoint from g, shouldn't happen!\nAbout to explode...")
            print("g\n",g)
            print("alpha\n",L[j])
        x = x.inv()
        y = FLD_ONE()
        for i in range(g.degree):
            # storing integers instead of GF objects
            H[i,j] = (y * x).coeffs
            y = L[j] * y
    return H

def restrict(H):
    """
    params: H over GF(2^DIM) (numpy integer array)
    return: H restricted to GF(2) = {0,1}, entries to columns
    """
    t = H.shape[0]
    R = np.zeros((t*DIM,H.shape[1]))
    for i in range(t):
        for j in range(CODE_DIM):
            col = int(H[i,j]) # this is an integer in [0, 2^DIM)
            for k in range(DIM):
                R[i*DIM + k, j] = get_bit(col, DIM - k -1)
    return R

# def fast_right_rref(M):
#     """
#     "fast" in that it exits early if not systematic-full-rank,
#     i.e. not in form (I|K)

#     params: matrix over GF(2) = {0, 1} (numpy array)
#     return: the non-identity part of rref(M) if the matrix is wider than it is tall
#             and has full rank, preceded by success Boolean,
#             either (False, None) or (True, public key)
#     """
#     m = M.shape[0] #len(M)
#     n = M.shape[1] #len(M[0])
#     # first do upper triangular to speed up non-systematic early abort
#     for i in range(m):
#         if M[i][i] == 0: # need a non-zero pivot
#             zero_pivot = True
#             for j in range(i+1,m): # look below for pivot
#                 if M[j][i] == 1: # switch rows
#                     tempi = np.array(M[i])
#                     tempj = np.array(M[j])
#                     M[i] = tempj
#                     M[j] = tempi
#                     zero_pivot = False
#                     break
#             if zero_pivot == True: # no non-zero pivot, not full-rank
#                 #print(f"no pivot in column {i}")
#                 return (False, None)
#         for j in range(i+1, m): # clear out below
#             if M[j][i] == 1: # add correct multiple of row i to row j
#                 tempj = np.array(M[j])
#                 M[j] = np.remainder(tempj + M[i], 2)

#     # now in systematic upper triangular
#     # clear out above
#     for i in range(m):
#         for j in range(i):
#             if M[j,i] == 1:
#                 tempj = np.array(M[j])
#                 M[j] = np.remainder(tempj + M[i], 2)

#     pubkey = np.array(M[:,m:])
#     #print(M[:, :m]) # should be identity
#     return (True, pubkey)

def fast_right_rref(M):
    """
    "fast" in that it exits early if not systematic-full-rank,
    i.e. not in form (I|K)

    Differs from the commented-out version above by using Boolean arrays instead of integer arrays,
    2-3 times faster.

    params: matrix over GF(2) = {0, 1} (numpy array)
    return: the non-identity part of rref(M) if the matrix is wider than it is tall
            and has full rank, preceded by success Boolean,
            either (False, None) or (True, public key)
    """
    m = M.shape[0] #len(M)
    n = M.shape[1] #len(M[0])
    M = np.array(M,bool)
    # first do upper triangular to speed up non-systematic early abort
    for i in range(m):
        if not M[i][i]: # need a non-zero pivot
            zero_pivot = True
            for j in range(i+1,m): # look below for pivot
                if M[j][i]: # switch rows
                    tempi = np.array(M[i])
                    tempj = np.array(M[j])
                    M[i] = tempj
                    M[j] = tempi
                    zero_pivot = False
                    break
            if zero_pivot: # no non-zero pivot, not full-rank
                #print(f"no pivot in column {i}")
                return (False, None)
        for j in range(i+1, m): # clear out below
            if M[j][i]: # add correct multiple of row i to row j
                tempj = np.array(M[j])
                M[j] = np.logical_xor(tempj, M[i])

    # now in systematic upper triangular
    # clear out above
    for i in range(m):
        for j in range(i):
            if M[j,i]:
                tempj = np.array(M[j])
                M[j] = np.logical_xor(tempj, M[i])

    pubkey = np.array(M[:,m:], int)
    #print(np.array(M[:, :m],int)) # should be identity
    return (True, pubkey)

def fast_right_rref_GF(M):
    """
    "fast" in that it aborts early; probably pretty slow.

    params: M a matrix (list of row lists) with GF entries
            intended dimensions DEGREE x (DEGREE + 1)
    return: (True, col) if rref(M)=(I|col), col a list of GF,
            else return (False, None)
    """
    m = len(M)
    n = len(M[0])
    # first do upper triangular to speed up non-systematic early abort
    #print("M to be reduced\n",M)
    for i in range(m):
        if M[i][i] == FLD_ZERO(): # need a non-zero pivot
            zero_pivot = True
            for j in range(i+1,m): # look below for pivot
                if M[j][i] != FLD_ZERO(): # switch rows
                    tempi = copy.deepcopy(M[i])
                    tempj = copy.deepcopy(M[j])
                    M[i] = tempj
                    M[j] = tempi
                    zero_pivot = False
                    break
            if zero_pivot == True: # no non-zero pivot, not full-rank
                return (False, None)
            #print("swap\n", M)

        # scale row i
        scalar = copy.deepcopy(M[i][i])
        for j in range(i,n):
            M[i][j] = scalar.inv()*M[i][j]
        #print("scale\n", M)

        # clear out below
        for j in range(i+1, m):
            tempj = copy.deepcopy(M[j])
            # add correct multiple of row i to row j
            if M[j][i] != FLD_ZERO():
                for k in range(i,n):
                    M[j][k] = tempj[k] + M[i][k]*tempj[i]

    #print("upper triangular\n", M)
    # now in systematic upper triangular with 1s diagonal
    # clear out above
    for i in range(m):
        tempi = copy.deepcopy(M[i])
        for j in range(i):
            tempj = copy.deepcopy(M[j])
            if M[j][i] != FLD_ZERO():
                for k in range(i,n):
                    M[j][k] = tempj[k] + tempi[k]*tempj[i]
        #print("clear above\n", M)

    # return last column
    #print("reduced\n",M)
    col = [M[i][n-1] for i in range(m)]
    return (True, col)


def key_gen(verbose = True, poly = "irred"):
    """
    return: goppa poly g
            support L
            public key K
    """
    if verbose:
        print("-"*25)
        print("key generation:\n")
    s1 = time()
    while True:
        t1 = time()
        L = new_support()
        t2 = time()
        if verbose:
            print(f"support in {t2-t1} seconds")
        t3 = time()
        if poly == "irred":
            g = new_irred_goppa()
        elif poly == "separable":
            g = new_separable_goppa(L)
        t4 = time()
        if verbose:
            print(f"goppa poly in {t4-t3} seconds")
        t5 = time()
        H = cannonical_pcheck(g, L)
        t6 = time()
        if verbose:
            print(f"cannon. pcheck in {t6-t5} seconds")
        t7 = time()
        M = restrict(H)
        good, K = fast_right_rref(M)
        t8 = time()
        if verbose:
            print(f"rref in {t8-t7} seconds")
        if good:
            s2 = time()
            if verbose:
                print(f"SUCCESS! keys created in {s2-s1} seconds")
                print("-"*25)
            return (g, L, K)
        else:
            if verbose:
                print("rejected...\n")

#------------------------------------------------------------------------------------------#
# encrypt
#------------------------------------------------------------------------------------------#

def encrypt(m, K):
    """
    params: plaintext binary message m of length CODE_DIM and weight DEGREE
            public key K
    return: ciphertext syndrome (I|K)*m
    """
    xcols = K.shape[1] # len(K[0]) # CODE_DIM - DIM*DEGREE
    rows = K.shape[0] # len(K) # DIM*DEGREE
    cipher = [0 for i in range(rows)]
    for i in range(rows):
        cipher[i] += m[i]
        for j in range(xcols):
            cipher[i] += ( m[rows + j] * K[i][j] )
    return [c % 2 for c in cipher]

#------------------------------------------------------------------------------------------#
# decrypt
#------------------------------------------------------------------------------------------#

def decrypt(c, g, L, decode = "alg6"):
    """
    double syndrome decoding, two algorithms available
    option = "alg6", "euclid"
    Patterson decoding not implemented (yet).
    """
    gsq = g*g
    H2 = cannonical_pcheck(gsq, L)
    H2r = restrict(H2)
    zeros = [0 for i in range(CODE_DIM - DIM*DEGREE)]
    c0 = c + zeros
    c2 = [0 for i in range(2*DIM*DEGREE)]
    for i in range(2*DIM*DEGREE):
        for j in range(CODE_DIM):
            c2[i] += int(H2r[i][j])*c0[j]
    c2 = [ int(c % 2) for c in c2]
    #print("cipher\n", c)
    #print("double cipher\n", c2)
    S = PolyFF(unrestrict_vector(c2))
    #print("syndrome polynomial\n", S)

    if decode == "alg6":
        sigma1 = alg6decode(S)
        #print("sigma1\n", sigma1)
        m1 = [0 for i in range(CODE_DIM)]
        for i in range(CODE_DIM): # evaluate sigma along the support L
            if L[i] == FLD_ZERO():
                if sigma1.degree < DEGREE:
                    m1[i] = 1
            else:
                ev = eval_poly(sigma1, L[i].inv())
                if ev == FLD_ZERO():
                    m1[i] = 1
        return m1

    if decode == "euclid":
        sigma2 = dbl_syn_decode(S)
        #print("sigma2\n", sigma2)
        m2 = [0 for i in range(CODE_DIM)]
        for i in range(CODE_DIM): # evaluate sigma along the support L
            if L[i] == FLD_ZERO():
                if sigma2.degree < DEGREE:
                    m2[i] = 1
            else:
                ev = eval_poly(sigma2, L[i].inv())
                if ev == FLD_ZERO():
                    m2[i] = 1
        return m2

def dbl_syn_decode(S):
    """
    NOTE: there are other ways to decode, this is the "double syndrome" decoding

    params: syndrome polynomial S
    return: error locator polynomial sigma
    """
    monomial = PolyFF([FLD_ZERO() for i in range(2*DEGREE)] + [FLD_ONE()])
    a,b,d = poly_EEA(monomial, S, early_stop = DEGREE)
    sigma = b # error locator polynomial
    #omega = d # error evaluation polynomial, not needed

    #print(a*monomial + b*S == d) # the euclidean algorithm works
    #print(f"sigma, degree {sigma.degree}\n", sigma) # should be degree DEGREE
    #print(f"omega, degree {omega.degree}\n", omega) # should be degree DEGREE -1

    return sigma


def alg6decode(S):
    """
    algorithm from https://eprint.iacr.org/2017/1180.pdf

    params: syndrome polynomial S
    return: error locator polynomial sigma
    """
    # initializations
    sigma = PolyFF([FLD_ONE()])
    beta = PolyFF([FLD_ZERO(),FLD_ONE()])
    l = 0
    delta = FLD_ONE()

    # supposed to construct sigma
    for k in range(2*DEGREE):
        # build d
        d = FLD_ZERO()
        for i in range(DEGREE+1):
            if k - i >= 0 and i <= sigma.degree:
                d = d + sigma.coeffs[i]*S.coeffs[k - i]
        # updates
        if d == FLD_ZERO() or k < 2*l:
            sigma = sigma - PolyFF([d*delta.inv()])*beta
            beta = PolyFF([FLD_ZERO(),FLD_ONE()])*beta
            #l = l
            #delta = delta
        else:
            sigma, beta = sigma - PolyFF([d*delta.inv()])*beta, PolyFF([FLD_ZERO(),FLD_ONE()])*sigma
            l = k - l + 1
            delta = d
        #print("partial sigma1\n", sigma)
    return sigma

def unrestrict_vector(x):
    """
    takes a vector with {0,1} coefficients and return a list of GF
    reading DIM entries at a time (i.e. using explicit linear isomorphism GF(2)^mn GF(2^m)^n, m = DIM)

    params: x a binary list of length DIM*n
    return: y, a list of GF of length n
    """
    y = []
    n = int(len(x)/DIM)
    for i in range(n):
        z = 0
        for j in range(DIM):
            z += x[i*DIM + j] * 2**(DIM - j - 1)
        y.append(GF(z))
    return y

# for debugging
def actual_sigma(L, m):
    """
    given message m and support L, return the error locator polynomial
    """
    prod = PolyFF([FLD_ONE()])
    for i in range(len(L)):
        if m[i] == 1:
            prod = prod * PolyFF([FLD_ONE(), L[i]])
    return prod

#------------------------------------------------------------------------------------------#
# testing, exporting keys to file
#------------------------------------------------------------------------------------------#

def random_message(weight = DEGREE):
    """
    return: random weight DEGREE binary message of length CODE_DIM
    """
    m = [0 for i in range(CODE_DIM)]
    perm = np.random.permutation(CODE_DIM)
    for i in range(weight):
        m[perm[i]] = 1
    return m

def random_trials(num_trials = 10, decode="alg6", poly = "irred"):
    g, L, K = key_gen(poly=poly)
    enc_tot = 0
    dec_tot = 0
    for i in range(num_trials):
        print(f"Trial {i+1}/{num_trials}")
        m = random_message()
        t1 = time()
        c = encrypt(m, K)
        t2 = time()
        print(f"encryption time: {t2-t1} seconds")
        enc_tot += t2 - t1
        t3 = time()
        m0 = decrypt(c, g, L, decode=decode)
        t4 = time()
        print(f"decryption time: {t4-t3} seconds")
        dec_tot += t4 - t3
        #print("sigma should be\n", actual_sigma(L, m))
        #print("message:\n", m)
        #print("cipher:\n", c)
        #print("decrypted:\n", m0)
        if m == m0:
            print("decryption SUCCESS\n")
        else:
            print("decryption FAILURE\n")
    print("Average times:")
    print("encryption: ", enc_tot/num_trials)
    print("decryption: ", dec_tot/num_trials)

def key_gen_to_file(name="generic"):
    """
    writes g, L, K to file.  just for show, can be compressed and standardized.

    params: name is the text file name; if generic records parameter info
    return: txt file with g, L, K
    """
    g,L,K = key_gen()
    if name == "generic":
        newname = f"m{DIM}n{CODE_DIM}t{DEGREE}.txt"
    else:
        newname = name + ".txt"
    with open(newname, 'w', encoding = 'utf-8') as f:
        f.write("Working over GF(2^m) with degree t Goppa polynomial and length n support.\n")
        f.write("Finite field elements are integers [0,2^m) whose binary expansion gives GF(2) coefficients via the field polynomial.\n")
        f.write(f"\nPARAMETERS:\nm = {DIM}\nn = {CODE_DIM}\nt = {DEGREE}\nfield polynomial (decreasing degree) = {MIN_POLY}\n")
        f.write(f"\nGOPPA POLYNOMIAL (degree {DEGREE}, increasing degree coefficients):\n")
        f.writelines([str(x.coeffs) + "\n" for x in g.coeffs])
        f.write(f"\nSUPPORT (length {CODE_DIM}):\n")
        f.writelines([str(l.coeffs) + "\n" for l in L])
        f.write(f"\nPUBLIC KEY ({DIM*DEGREE}-by-{CODE_DIM-DIM*DEGREE} binary matrix):\n")
        klines = []
        for i in range(K.shape[0]):
            row = ""
            for j in range(K.shape[1]):
                row += str(int(K[i,j]))
            klines.append(row)
        f.writelines([r + "\n" for r in klines])

def key_gen_timer():
    """
    return: goppa poly g
            support L
            public key K
    """
    s1 = time()
    while True:
        L = new_support()

        t3 = time()
        g = new_irred_goppa()
        t4 = time()

        t5 = time()
        H = cannonical_pcheck(g, L)
        t6 = time()

        t7 = time()
        M = restrict(H)
        good, K = fast_right_rref(M)
        t8 = time()
        if good:
            s2 = time()
            return (g, L, K, t4-t3, t6-t5, t8-t7, s2-s1)
