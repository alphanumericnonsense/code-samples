import copy

def gcd(a,b):
    while b!=0:
        a,b=b,a%b
    return abs(a)

def lcm(args):
    """
    lcm using gcd
    """
    args = [arg for arg in args if arg != 0]
    retval = args[0]
    for a in args:
        retval = retval*a//gcd(retval, a)
    return abs(retval)

class frac(object):
    """
    simple rational number class
    """
    def __init__(self,num,den):
        assert type(num)==type(den)==int, "numerator and denominator not integers"
        assert den!=0, "division by zero"
        d=gcd(num,den)
        num=num//d
        den=den//d
        self.num=num
        self.den=den
    def __add__(self,other):
        a,b,c,d=self.num,self.den,other.num,other.den
        return frac(a*d+b*c,b*d)
    def __sub__(self,other):
        a,b,c,d=self.num,self.den,other.num,other.den
        return frac(a*d-b*c,b*d)
    def __mul__(self,other):
        a,b,c,d=self.num,self.den,other.num,other.den
        return frac(a*c,b*d)
    def __truediv__(self,other): # not in python 2.7???
    #def __div__(self,other):
        a,b,c,d=self.num,self.den,other.num,other.den
        return frac(a*d,b*c)
    def __str__(self):
        return str(self.num)+"/"+str(self.den)
    def __repr__(self):
        return str(self)
    def __eq__(self, other):
        a,b,c,d = self.num,self.den,other.num,other.den
        if a*d-b*c == 0:
            return True
        else:
            return False
    def __ne__(self, other):
        return not (self == other)

def m2frac(m):
    """
    convert input to stochastic matrix and return terminal states
    """
    new_m = []
    terminal_states = []
    for i, row in enumerate(m):
        new_row = []
        row_sum = sum(row)
        if row_sum == 0:
            new_row = [frac(1,1) if j == i else frac(0,1) for j in range(len(m))]
            terminal_states.append(i)
        else:
            for entry in row:
                new_row.append(frac(entry, row_sum))
        new_m.append(new_row)
    return (new_m, terminal_states)

def frac_mat_mult(m1, m2):
    """
    matrix multiplication over rationals, no error checking
    NxM times MxL
    """
    N = len(m1)
    M = len(m1[0])
    L = len(m2[0])
    result = []
    for i in range(N):
        new_row = []
        for j in range(L):
            eij = frac(0,1)
            for k in range(M):
                eij = eij + m1[i][k]*m2[k][j]
            new_row.append(eij)
        result.append(new_row)
    return result

def frac_mat_add(m1, m2):
    """
    matrix addition, no error checking
    """
    return [[m1[i][j]+m2[i][j] for j in range(len(m1[0]))] for i in range(len(m1))]

def frac_mat_sub(m1, m2):
    """
    matrix addition, no error checking
    """
    return [[m1[i][j]-m2[i][j] for j in range(len(m1[0]))] for i in range(len(m1))]

def frac_mat_rref(M):
    """
    rref for use in inverse
    """
    # matrix dimensions
    m = len(M)
    n = len(M[0])
    count = 0
    for j in range(n):
        found = False
        for i in range(count,m):
            if M[i][j] != frac(0,1):
                found=True
                pr = i
                pc = j
                break
        if found == True:
            M[count], M[pr] = M[pr], M[count] # switch pivot row and count row
            coeff = M[count][pc]
            for j1 in range(n):
                M[count][j1] = M[count][j1]/coeff # scale count row so first non-zero entry is 1
            for i in range(m):#clears out pivot column outside count row (elementary row operations)
                if i != count:
                    coeff1 = M[i][pc]
                    for j2 in range(n):
                        M[i][j2] = M[i][j2] - coeff1*M[count][j2]
            count += 1
    return M

def frac_mat_inv(M):
    """
    inverse of M
    """
    # create augmented (M|I) to row reduce to (I|Minv)
    n = len(M)
    augM = []
    for i in range(n):
        augM.append([])
        for j in range(2*n):
            if j < n:
                augM[i].append(M[i][j])
            elif i == j-n:
                augM[i].append(frac(1,1))
            else:
                augM[i].append(frac(0,1))
    augMinv = frac_mat_rref(augM)
    # chop off second half
    Minv = []
    for i in range(n):
        Minv.append([])
        for j in range(n):
            Minv[i].append(augMinv[i][j+n])
    return Minv

def print_mat(m):
    for row in m:
        print(row)
    print("\n")

import random
import math
def rand_mat():
    """
    random example with n > 1 and 1 < k < n terminal states
    """
    n = math.floor(1+10*random.random())
    k = math.floor(1+(n-1)*random.random())
    retmat = []
    for i in range(n-k):
        retmat.append([])
        for j in range(n):
            x = math.floor(100*random.random())
            retmat[i].append(x)
    for i in range(k):
        retmat.append([0 for j in range(n)])
    if retmat == [[0]]:
        print(n, k)
    return retmat

def solution(m):
    """
    absorbing markov chain,
    https://en.wikipedia.org/wiki/Absorbing_Markov_chain

    m~ = (Q|R)
         (0|I)

    N = 1/(1-Q), NR describes transition from transient to absorbing states.
    Want first row of NR (probability to go from state 0 to absorbing states)???

    Is Q nilpotent? No...
    """
    # "normalize" and split m
    M, tstates = m2frac(m)
    tcompl = [i for i in range(len(M)) if i not in tstates]

    if tcompl == []:
        return [1,1] # ???

    Q = [[M[i][j] for j in tcompl] for i in tcompl]
    R = [[M[i][j] for j in tstates] for i in tcompl]
    #print_mat(Q)
    #print_mat(R)
    I = [[frac(1,1) if i == j else frac(0,1) for j in range(len(tcompl))] for i in range(len(tcompl))]
    N = frac_mat_inv(frac_mat_sub(I, Q))
    probvec = frac_mat_mult(N, R)[0]
    # convert to output format...
    bigD = lcm([x.den for x in probvec])
    bigD_frac = frac(bigD,1)
    scale_output = [bigD_frac*x for x in probvec]
    integerify = [x.num for x in scale_output]
    return integerify + [bigD]
