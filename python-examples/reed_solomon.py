#------------------------------------------------------------------------------------------#
# Imports, constants
#------------------------------------------------------------------------------------------#
import numpy as np
from random import random
from time import time
import copy
import os
from math import floor, ceil, log2
#---------------------------#
# Some minimal polynomials
#---------------------------#
# first number is degree
# following numbers are nonzero degrees
# except for lead (monic) and
# constant term 1
# 2,1 : x^2 + x + 1
# 3,1 : x^3 + x + 1
# 4,1
# 5,2
# 6,1
# 7,1
# 8,4,3,1 : x^8 + x^4 + x^3 + x + 1, primitive element 3 = x+1
# 9,1, primitive element 9 = x^3+1
# 10,3
# 11,2
# 12,3
# 13,4,3,1
# 14,5
# 15,1
# 16,5,3,1
# 17,3
# 18,3
# 19,5,2,1
# 20,3
# 21,2
# 22,1
# 23,5
# 24,4,3,1
# 25,3
# 26,4,3,1
# 27,5,2,1
# 28,1
# 29,2
# 30,1
#---------------------------#
# some parameters to set
#---------------------------#
MIN_POLY = "100011011" # monic degree DIM to degree 0
DIM = 8
CODE_DIM_OUT = 2**DIM-1
CODE_DIM_IN = 32
MASK = 2**DIM - 1 # used as a bitmask
REDUCER = int(MIN_POLY[1:], 2) # used in finite field multiplication

#------------------------------------------------------------------------------------------#
# GF(2^DIM) class
#
# Finite field elements are integers in [0, 2^DIM) and bitwise operations are used,
# e.g. 0 = 0, 1 = 1, x = 2, x+1 = 3, x^2 = 4, ..., x^{DIM-1}+x^{DIM-2}+...+x+1 = 2^{DIM}-1
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

##################################################
# Inverse look-up table!!!
##################################################
INVERSES = inv_dict() # generates LUT for inverses
##################################################

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
            string += " : "
            string += str(self.coeffs[i])
            string +="\n"
        return string

    def __copy__(self):
        """
        an attempt at copying?
        """
        new_coeffs = []
        for gf in self.coeffs:
            new_gf_coeffs = copy.copy(gf.coeffs)
            new_coeffs.append(GF(new_gf_coeffs))
        return PolyFF(new_coeffs)

    def scale(self, alpha):
        """
        scale poly by alpha in GF
        """
        return PolyFF([alpha*c for c in self.coeffs])

    def differentiate(self):
        """
        derivative of polynomial
        """
        return PolyFF([GF(i%2)*self.coeffs[i] for i in range(self.degree+1)][1:])

# useful functions using/used by PolyFF

def mod_prod(a, b, f):
    """
    used in new_irred_goppa()

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

def eval_poly(f, x):
    """
    evaluate polynomial f at point x in GF(q) in an iterative (Horner) fashion

    deprecated in favor of AFFT(f) which evaluates f over the whole field with AFFT
    """
    s = FLD_ZERO()
    d = f.degree
    for i in range(d + 1):
        s = s*x + f.coeffs[d - i]
    return s

def monomial(k):
    """
    return x^k = PolyFF([GF(0),...,GF(0),GF(1)])
    """
    zeros = [GF(0) for i in range(k)]
    zeros.append(GF(1))
    return PolyFF(zeros)

#------------------------------------------------------------------------------------------#
# Additive FFT evaluating a polynomial on all of GF(2^m)
# following http://www.math.clemson.edu/~sgao/papers/GM10.pdf
#------------------------------------------------------------------------------------------#

def FIELD_LIST():
    """
    a list of all the field elements
    """
    return [GF(i) for i in range(2**DIM)]

# precomputed list of all field elements
THE_FIELD = FIELD_LIST()

# tested and seems to work
def radix_conversion(f, n):
    """
    expanding f in base x^2+x recursively, f = g0 + g1*(x^2+x)^(2^k) below

    params: f a polynomial over GF(2^m) with deg(f) < n, 2^(k+1) < n <= 2^(k+2)
    """

    if n <= 2: # degree <= 1 base case
        d = f.degree
        if d == 1:
            a = f.coeffs[0]
            b = f.coeffs[1]
        elif d == 0:
            a = f.coeffs[0]
            b = GF(0)
        else:
            a = GF(0)
            b = GF(0)
        return [(a,b)]

    k = ceil(log2(n)) - 2
    f0 = PolyFF(f.coeffs[ : 2**(k+1)])
    f1 = PolyFF(f.coeffs[2**(k+1) : 2**k + 2**(k+1)])
    f2 = PolyFF(f.coeffs[2**k + 2**(k+1) : ])
    g0 = f0 + (f1 + f2)*monomial(2**k)
    g1 = f1 + f2 + f2*monomial(2**k)
    V1 = radix_conversion(g0, 2**(k+1))
    V2 = radix_conversion(g1, n-2**(k+1))
    return V1 + V2

# tested, works
def AFFT(f, m = DIM, B = THE_FIELD, basis = [GF(2**i) for i in range(DIM)]):
    """
    computes the value of f on all of GF(2^m) recursively,
    can be used on any GF(2) subspace of GH(2^m)

    params: B a subspace with basis, initialized to the whole field
            f the poly we're evaluating
            m the dimension of the subspace over GF(2)
            basis - the ordered basis for the subspace,
                    initialized to the whole field with power basis 1, x, x^2, ..., x^(m-1)
    return: list [f(GF(0)), ..., f(GF(2^m-1))]
    """
    if m == 1:
        return [eval_poly(f,GF(0)), eval_poly(f,basis[0])] # GF(2) base case

    G = [beta*(basis[m-1].inv()) for beta in B[:2**(m-1)]] # half-dimensional "twisted" subspace
    Gbasis = [beta*(basis[m-1].inv()) for beta in basis[:-1]] # basis for G
    D = [gamma**2 + gamma for gamma in G]
    Dbasis = [gamma**2 + gamma for gamma in Gbasis] # basis for D

    ftwist = PolyFF([f.coeffs[i]*basis[m-1]**i for i in range(f.degree + 1)]) # "twist" by last element of basis
    glist = radix_conversion(ftwist, ftwist.degree + 1) # [..., (gi0, gi1), ...]

    g0 = PolyFF([gi[0] for gi in glist])
    g1 = PolyFF([gi[1] for gi in glist])
    U = AFFT(g0, m - 1, D, Dbasis) # recursion
    V = AFFT(g1, m - 1, D, Dbasis) # recursion

    Wfirst = [U[i]+G[i]*V[i] for i in range(2**(m-1))] # ftwist on G
    Wlast = [Wfirst[i] + V[i] for i in range(2**(m-1))] # ftwist on affine subspace G + 1
    return Wfirst + Wlast

def field_eval(f):
    """
    evaluate poly f on all of GF(2^DIM) conventionally,
    for testing the AFFT
    """
    return [eval_poly(f, GF(i)) for i in range(2**DIM)]

#-----------------------------------------------------------------#
# some linear algebra I don't use
#-----------------------------------------------------------------#
def fast_right_rref(M):
    """
    "fast" in that it exits early if not systematic-full-rank,
    i.e. not in form (I|K)

    Uses Boolean arrays instead of integer arrays, 2-3 times faster,
    but still slowest part of entire scheme

    params: matrix over GF(2) = {0, 1} (numpy array)
    return: the non-identity part of rref(M) if the matrix is wider than it is tall
            and has full rank, preceded by success Boolean,
            either (False, None) or (True, public key)
    """
    m = M.shape[0] #len(M)
    n = M.shape[1] #len(M[0])
    M = np.array(M, bool)
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

    # now in systematic upper triangular with 1s diagonal
    # clear out above
    for i in range(m):
        tempi = copy.deepcopy(M[i])
        for j in range(i):
            tempj = copy.deepcopy(M[j])
            if M[j][i] != FLD_ZERO():
                for k in range(i,n):
                    M[j][k] = tempj[k] + tempi[k]*tempj[i]

    # return last column
    col = [M[i][n-1] for i in range(m)]
    return (True, col)

    # # modify to produce inverse
    # M = [[GF(PRIM_ELT)**(i*j) for j in range(256)]+[GF(1) if k==i else GF(0) for k in range(256)] for i in range(256)]
    # Minv = fast_right_rref_GF(M)
    # print(Minv)

#-----------------------------------------------------------------#
# decoding stuff for RS codes
#-----------------------------------------------------------------#
# Useful for RS codes
PRIM_ELT = 3 # primitive element for the field
UNITS_REV = {(GF(PRIM_ELT)**i).coeffs:i for i in range(MASK)}
UNITS_FWD = {i:(GF(PRIM_ELT)**i).coeffs for i in range(MASK)}

def BM(S):
#     """
#     Return connection poly for shortest LFSR producing sequence S,
#     sum_i C_i*S_{N-i} = 0 for L <= N < n, C_0 = 1, [1/C = S?]
#     """
#     C = PolyFF([GF(1)]) # poly over field, connection poly to be returned
#     B = PolyFF([GF(1)]) #poly over field
#     x = 1 # integer
#     L = 0 # integer, length of "connection poly"/recurrence relation/LFSR
#     b = GF(1) # field elt
#     N = 0 # integer, counts 0 to n - 1
#     n = len(S.coeffs)

#     while N < n:
#         d = GF(0) # field elt, discrepancy
#         for i in range(0,L+1):
#             d = d + C.coeffs[i]*S.coeffs[N-i]
#         if d == GF(0):
#             x += 1
#         elif 2*L > N:
#             C = C - (B*monomial(x)).scale(d*b.inv())
#             x += 1
#         else:
#             B = copy.copy(C)
#             C = C - (B*monomial(x)).scale(d*b.inv())
#             L = N + 1 - L
#             b = d
#             x = 1
#         N += 1
#     return C
    # initializations
    sigma = PolyFF([FLD_ONE()])
    beta = PolyFF([FLD_ZERO(),FLD_ONE()])
    l = 0
    delta = FLD_ONE()
    N = len(S.coeffs)
    # supposed to construct sigma
    for k in range(N):
        # build d
        d = FLD_ZERO()
        for i in range(floor((CODE_DIM_OUT-CODE_DIM_IN)/2)):
            if k - i >= 0 and i <= sigma.degree:
                d = d + sigma.coeffs[i]*S.coeffs[k - i]

        # updates
        if d == FLD_ZERO() or k < 2*l:
            sigma = sigma - PolyFF([d*delta.inv()])*beta
            beta = PolyFF([GF(0)] + beta.coeffs) # beta = PolyFF([FLD_ZERO(),FLD_ONE()])*beta
            #l = l
            #delta = delta
        else:
            sigma, beta = sigma - PolyFF([d*delta.inv()])*beta, PolyFF([GF(0)] + sigma.coeffs) # PolyFF([FLD_ZERO(),FLD_ONE()])*sigma
            l = k - l + 1
            delta = d
        #print("partial sigma1\n", sigma)
    return sigma

def encode(message_bytes):
    """
    Assuming message bytes is a multiple of CODE_DIM_IN bytes (block size for ECC).
    RS-code with parameters [n=CODE_DIM_OUT,k=CODE_DIM_IN,d=n-k+1].
    """
    n_blocks = len(message_bytes)//CODE_DIM_IN
    encoded = bytes([])
    for n in range(n_blocks):
        a = message_bytes[n*CODE_DIM_IN:(n+1)*CODE_DIM_IN]
        #a_poly = PolyFF(223*[GF(0)] + [GF(int(a[CODE_DIM_IN-1-i])) for i in range(CODE_DIM_IN)]) # a[0]*x^(n-1) + ...
        a_poly = PolyFF([GF(int(a[i])) for i in range(CODE_DIM_IN)])
        encoded += bytes([eval_poly(a_poly, GF(PRIM_ELT)**e).coeffs for e in range(CODE_DIM_OUT)])
    return encoded

def decode(noisy_bytes):
    """
    syndrome decoding using BM, then invert encoding.
    """
    n_blocks = len(noisy_bytes)//CODE_DIM_OUT
    decoded = bytes([])
    for n in range(n_blocks):
        # get error locator sigma from syndrome for this block
        c_noise = noisy_bytes[n*CODE_DIM_OUT:(n+1)*CODE_DIM_OUT]
        c_noise_poly = PolyFF([GF(int(c)) for c in c_noise])
        syn = PolyFF([eval_poly(c_noise_poly, GF(PRIM_ELT)**e) for e in range(CODE_DIM_OUT-CODE_DIM_IN)])
        sigma = BM(syn) # error locator polynomial, roots are inverses of error locations
        #print("sigma degree", sigma.degree)

        # error evaluator polynomial-ish?
        omega = syn*sigma
        #m = min(omega.degree+1, 2*sigma.degree)
        omega = PolyFF(omega.coeffs[:2*sigma.degree])
        #print("omega degree", omega.degree)

        # find error locations
        error_locs = []
        e = 0
        t = sigma.degree
        for i in range(CODE_DIM_OUT):
            pot_root = GF(PRIM_ELT)**i
            if eval_poly(sigma, pot_root) == GF(0):
                error_locs.append(UNITS_REV[pot_root.inv().coeffs])
                e += 1
            if e >= t:
                break
        #print("error_locs", error_locs)
        error_coeffs = [GF(0) for i in range(CODE_DIM_OUT)]

        # find error values
        for index in error_locs:
            x = GF(UNITS_FWD[index])
            error_coeffs[index] = -x*eval_poly(omega, x.inv())*eval_poly(sigma.differentiate(), x.inv()).inv()
#         for i in range(len(error_coeffs)):
#             if error_coeffs[i] != GF(0):
#                 print(i,error_coeffs[i])

        # remove error
        c_poly = c_noise_poly - PolyFF(error_coeffs)
        # project back to message space
        a_gf_list = [eval_poly(c_poly, GF(PRIM_ELT)**(CODE_DIM_OUT-i)) for i in range(CODE_DIM_IN)]
        # tack on bytes for this block
        decoded += bytes([a_gf_list[i].coeffs for i in range(CODE_DIM_IN)])
    return decoded

if __name__ == "__main__":
    N = 10
    for i in range(N):
        #print(f"\nerror at position {i}, {i+1}, {i+2}")
        rmsg = os.urandom(32)
        c = encode(rmsg)
        err = os.urandom(100)
        rerr = bytes([0]*i) + err + bytes([0]*(255-i-100))
        rnoise = bytes([c[j] ^ rerr[j] for j in range(255)])
        #print(bin(err[0]), bin(err[1]), bin(err[2]))
        d = decode(rnoise)
        print(rmsg.hex())
        print(d.hex())
        print(rmsg==d)
