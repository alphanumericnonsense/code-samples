#------------------------------------------------------------------------------------------#
# Additive FFT evaluating a polynomial on all of GF(2^m)
# following http://www.math.clemson.edu/~sgao/papers/GM10.pdf
#
# To be used to speed up evaluation of Goppa polynomial and error locator polynomial
#------------------------------------------------------------------------------------------#

import math

#---------------------------#
# some fields for testing
#---------------------------#

# GF(4)
# for super testing
# MIN_POLY = "111"
# DIM = 2

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

# GF(2**13)
MIN_POLY = "10000000011011" # monic degree DIM to degree 0, left-to-right
DIM = 13

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

def eval_poly(f, x):
    """
    evaluate polynomial f at point x in GF(q) in an iterative (Horner) fashion
    """
    s = FLD_ZERO()
    d = f.degree
    for i in range(d + 1):
        s = s*x + f.coeffs[d - i]
    return s

def FIELD_LIST():
    """
    a list of all the field elements
    """
    return [GF(i) for i in range(2**DIM)]

# precomputed list of all field elements
THE_FIELD = FIELD_LIST()

def monomial(n):
    return PolyFF([FLD_ZERO() for i in range(n)] + [FLD_ONE()])

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

    k = math.ceil(math.log2(n)) - 2
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

if __name__=="__main__":
    import math
    import random
    import time

    testdegree = 100
    num_trials = 5
    print(f"Testing evaluation of {num_trials} random degree {testdegree} polynomials over GF(2**{DIM})")
    for i in range(num_trials):
        testpolylist = []
        for i in range(testdegree+1):
            coeff = math.floor((2**DIM)*random.random())
            testpolylist.append(GF(coeff))
        testpoly = PolyFF(testpolylist)

        # print(f"testpoly:\n{testpoly}")
        print(f"Random testpoly {i+1}")
        print("fast eval:")
        t0 = time.time()
        X = AFFT(testpoly)
        t1 = time.time()
        print(f"time:  {t1-t0}")

        print("conventional eval:")
        t0 = time.time()
        Y = field_eval(testpoly)
        t1 = time.time()
        print(f"time:  {t1-t0}")

        if X == Y:
            print("Evaluations AGREE.\n")
        else:
            print("Evaluations DISAGREE!\n")
